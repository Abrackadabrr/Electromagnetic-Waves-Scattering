// Comparative matvec benchmark for five block-Toeplitz implementations.
//
// Every Eigen leaf matrix has dynamic dimensions and owns its coefficient
// allocation. BlockSize still selects its runtime square dimensions.
// Example:
//   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 \
//     ./bench_five_toeplitz_matvec_layouts 2>&1 | tee results.txt
// Iteration count and calibration time are determined by Google Benchmark.
// Each configuration uses ten independent runs to estimate mean and stddev.

#ifndef EIGEN_STACK_ALLOCATION_LIMIT
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#endif

#include <benchmark/benchmark.h>

#include "math/matrix/TeoplitzTemplateNew.hpp"
#include "math/matrix/ToeplitzFullyTemplated.hpp"
#include "math/matrix/ToeplitzMultilevelNew.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

namespace Matrix = EMW::Math::LinAgl::Matrix;

using Scalar = std::complex<double>;
using Index = Eigen::Index;
using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
using DynamicBlock = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

// Keeping one alias per registered size makes switching selected sizes back to
// fixed-size leaves a one-line change, while this benchmark uses heap-allocated
// storage for every leaf.
using Block3 = DynamicBlock;
using Block5 = DynamicBlock;
using Block16 = DynamicBlock;
using Block32 = DynamicBlock;
using Block64 = DynamicBlock;
using Block128 = DynamicBlock;
using Block256 = DynamicBlock;
using Block1024 = DynamicBlock;

constexpr std::size_t mebibyte = 1024ULL * 1024ULL;
constexpr std::size_t maximum_toeplitz_coefficients = 1200ULL * mebibyte;
constexpr std::size_t construction_headroom = 128ULL * mebibyte;
constexpr int independent_repetitions = 10;

// The original shapes below all contain 16 logical blocks. The additional
// Level3FiveSettings case contains 5^3 = 125 logical blocks.
template <int Levels> struct LevelSettings;
template <> struct LevelSettings<1> {
    inline static constexpr std::array<Index, 1> extents{16};
    inline static constexpr const char *label = "L1";
};
template <> struct LevelSettings<2> {
    inline static constexpr std::array<Index, 2> extents{4, 4};
    inline static constexpr const char *label = "L2";
};
template <> struct LevelSettings<3> {
    inline static constexpr std::array<Index, 3> extents{2, 2, 4};
    inline static constexpr const char *label = "L3";
};

struct Level3FiveSettings {
    inline static constexpr std::array<Index, 3> extents{5, 5, 5};
    inline static constexpr const char *label = "L3x5x5x5";
};

template <typename Block, int Levels> struct OldNested {
    using type = Matrix::ToeplitzStructure<Scalar, typename OldNested<Block, Levels - 1>::type>;
};
template <typename Block> struct OldNested<Block, 0> { using type = Block; };
template <typename Block, int Levels> using OldNestedT = typename OldNested<Block, Levels>::type;

template <typename Block, int Levels> struct NewNested {
    using type = Matrix::ToeplitzTemplateNew<typename NewNested<Block, Levels - 1>::type>;
};
template <typename Block> struct NewNested<Block, 0> { using type = Block; };
template <typename Block, int Levels> using NewNestedT = typename NewNested<Block, Levels>::type;

std::size_t available_memory_bytes() {
    std::ifstream meminfo("/proc/meminfo");
    std::string key, unit;
    std::size_t value = 0;
    while (meminfo >> key >> value >> unit)
        if (key == "MemAvailable:")
            return value * 1024ULL;
    return 0;
}

std::size_t initial_available_memory_bytes() {
    // Some allocators retain freed arenas after a large benchmark case. Using
    // the live MemAvailable value would therefore make later registrations
    // fail the guard solely because of their execution order.
    static const std::size_t available = available_memory_bytes();
    return available;
}

template <int Levels, typename Settings = LevelSettings<Levels>>
constexpr Index unique_leaf_blocks() {
    Index result = 1;
    for (const Index extent : Settings::extents)
        result *= 2 * extent - 1;
    return result;
}

template <int BlockSize, int Levels, typename Settings = LevelSettings<Levels>>
constexpr std::size_t coefficient_bytes() {
    return static_cast<std::size_t>(unique_leaf_blocks<Levels, Settings>()) *
           BlockSize * BlockSize * sizeof(Scalar);
}

template <int BlockSize, int Levels, typename Settings = LevelSettings<Levels>>
bool configuration_fits_memory() {
    constexpr std::size_t coefficients = coefficient_bytes<BlockSize, Levels, Settings>();
    if (coefficients > maximum_toeplitz_coefficients)
        return false;
    // Packed constructors temporarily hold source blocks and the final arena.
    const std::size_t estimated_peak = 2 * coefficients + construction_headroom;
    const std::size_t available = initial_available_memory_bytes();
    // Keep at least 20% of initially available memory outside this estimate.
    return available == 0 || estimated_peak < available * 4 / 5;
}

template <int Levels>
std::pair<Scalar, Scalar> block_coefficients(std::span<const Index, Levels> row,
                                             std::span<const Index, Levels> col) {
    double real_displacement = 0.0;
    double imaginary_displacement = 0.0;
    double weight = 1.0;
    for (int level = 0; level < Levels; ++level) {
        const double displacement = static_cast<double>(col[level] - row[level]);
        real_displacement += weight * displacement;
        imaginary_displacement += (level + 1.0) * displacement;
        weight *= 0.37;
    }
    const Scalar diagonal{1.0 + 0.013 * real_displacement,
                          0.007 * imaginary_displacement};
    const Scalar constant{0.0003 * real_displacement,
                          -0.0002 * imaginary_displacement};
    return {diagonal, constant};
}

template <typename Block, int BlockSize, int Levels>
void initialize_block(Block &result, std::span<const Index, Levels> row,
                      std::span<const Index, Levels> col) {
    const auto [diagonal, constant] = block_coefficients<Levels>(row, col);
    result.setConstant(constant);
    result.diagonal().array() += diagonal;
}

template <typename Block, int BlockSize, int Levels>
void append_block(std::vector<Block> &blocks, std::span<const Index, Levels> row,
                  std::span<const Index, Levels> col) {
    // Construct the leaf in its final vector slot, then allocate its dynamic
    // coefficient buffer once at the requested runtime dimensions.
    blocks.emplace_back(BlockSize, BlockSize);
    initialize_block<Block, BlockSize, Levels>(blocks.back(), row, col);
}

template <typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>, int Depth = 0>
OldNestedT<Block, Levels - Depth>
make_old_subtree(std::array<Index, Levels> row, std::array<Index, Levels> col) {
    constexpr Index extent = Settings::extents[Depth];
    using Node = OldNestedT<Block, Levels - Depth>;
    using Child = typename Node::block_type;
    std::vector<Child> blocks;
    blocks.reserve(static_cast<std::size_t>(2 * extent - 1));
    if constexpr (Depth + 1 == Levels) {
        const auto append = [&](Index i, Index j) {
            row[Depth] = i;
            col[Depth] = j;
            append_block<Block, BlockSize, Levels>(blocks, row, col);
        };
        for (Index col_index = 0; col_index < extent; ++col_index)
            append(0, col_index);
        for (Index row_index = 1; row_index < extent; ++row_index)
            append(row_index, 0);
    } else {
        const auto append = [&](Index i, Index j) {
            row[Depth] = i;
            col[Depth] = j;
            blocks.emplace_back(
                make_old_subtree<Block, BlockSize, Levels, Settings, Depth + 1>(row, col));
        };
        for (Index col_index = 0; col_index < extent; ++col_index)
            append(0, col_index);
        for (Index row_index = 1; row_index < extent; ++row_index)
            append(row_index, 0);
    }
    return Node(extent, extent, std::move(blocks));
}

template <typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>>
OldNestedT<Block, Levels> make_old_matrix() {
    return make_old_subtree<Block, BlockSize, Levels, Settings>({}, {});
}

template <typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>, int Depth = 0>
NewNestedT<Block, Levels - Depth>
make_new_subtree(std::array<Index, Levels> row, std::array<Index, Levels> col) {
    constexpr Index extent = Settings::extents[Depth];
    using Node = NewNestedT<Block, Levels - Depth>;
    using Child = typename Node::block_type;
    std::vector<Child> blocks;
    blocks.reserve(static_cast<std::size_t>(2 * extent - 1));
    if constexpr (Depth + 1 == Levels) {
        const auto append = [&](Index i, Index j) {
            row[Depth] = i;
            col[Depth] = j;
            append_block<Block, BlockSize, Levels>(blocks, row, col);
        };
        for (Index col_index = 0; col_index < extent; ++col_index)
            append(0, col_index);
        for (Index row_index = 1; row_index < extent; ++row_index)
            append(row_index, 0);
    } else {
        const auto append = [&](Index i, Index j) {
            row[Depth] = i;
            col[Depth] = j;
            blocks.emplace_back(
                make_new_subtree<Block, BlockSize, Levels, Settings, Depth + 1>(row, col));
        };
        for (Index col_index = 0; col_index < extent; ++col_index)
            append(0, col_index);
        for (Index row_index = 1; row_index < extent; ++row_index)
            append(row_index, 0);
    }
    return Node(extent, extent, std::move(blocks));
}

template <typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>>
NewNestedT<Block, Levels> make_new_matrix() {
    return make_new_subtree<Block, BlockSize, Levels, Settings>({}, {});
}

template <typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>>
Matrix::MultilevelToeplitz<Block> make_runtime_matrix() {
    const auto &extents = Settings::extents;
    std::vector<Index> runtime_extents(extents.begin(), extents.end());
    std::vector<Block> blocks;
    blocks.reserve(static_cast<std::size_t>(unique_leaf_blocks<Levels, Settings>()));
    std::array<Index, Levels> row{}, col{};
    for (Index leaf = 0; leaf < unique_leaf_blocks<Levels, Settings>(); ++leaf) {
        Index encoded = leaf;
        for (int reverse = Levels; reverse-- != 0;) {
            const Index unique = 2 * extents[reverse] - 1;
            const Index diagonal = encoded % unique;
            encoded /= unique;
            if (diagonal < extents[reverse]) {
                row[reverse] = 0;
                col[reverse] = diagonal;
            } else {
                row[reverse] = diagonal - extents[reverse] + 1;
                col[reverse] = 0;
            }
        }
        append_block<Block, BlockSize, Levels>(blocks, row, col);
    }
    return Matrix::MultilevelToeplitz<Block>(Levels, runtime_extents,
                                              runtime_extents, std::move(blocks));
}

template <typename T>
concept OldToeplitzNode = requires(const T &matrix, Index i) {
    matrix.rows_in_toeplitrz();
    matrix.cols_in_toeplitrz();
    matrix.rows_in_block();
    matrix.cols_in_block();
    matrix.get_block(i, i);
};

// This is the allocation-free form of the original non-batched algorithm:
// identical row/column block traversal and identical legacy storage, with the
// result supplied by the caller. ToeplitzStructure's in-place overload only
// recognizes MatrixX leaves and therefore cannot be instantiated for the
// fixed-size Eigen leaves used by this benchmark.
template <OldToeplitzNode Node>
void old_classic_multiply_add(const Node &matrix, const Scalar *rhs, Scalar *dest) {
    const Index block_rows = static_cast<Index>(matrix.rows_in_toeplitrz());
    const Index block_cols = static_cast<Index>(matrix.cols_in_toeplitrz());
    const Index rows_in_block = static_cast<Index>(matrix.rows_in_block());
    const Index cols_in_block = static_cast<Index>(matrix.cols_in_block());
    for (Index row = 0; row < block_rows; ++row)
        for (Index col = 0; col < block_cols; ++col) {
            const auto &block = matrix.get_block(row, col);
            if constexpr (OldToeplitzNode<std::remove_cvref_t<decltype(block)>>) {
                old_classic_multiply_add(block, rhs + col * cols_in_block,
                                         dest + row * rows_in_block);
            } else {
                Eigen::Map<const Vector> input(rhs + col * cols_in_block, cols_in_block);
                Eigen::Map<Vector> output(dest + row * rows_in_block, rows_in_block);
                output.noalias() += block * input;
            }
        }
}

template <int Levels, typename Settings = LevelSettings<Levels>>
void decode_logical_index(Index linear, std::array<Index, Levels> &coordinates) {
    for (int reverse = Levels; reverse-- != 0;) {
        const Index extent = Settings::extents[reverse];
        coordinates[reverse] = linear % extent;
        linear /= extent;
    }
}

template <int Levels, typename Settings = LevelSettings<Levels>>
constexpr Index logical_extent() {
    Index result = 1;
    for (const Index extent : Settings::extents)
        result *= extent;
    return result;
}

template <int BlockSize, int Levels, typename Settings = LevelSettings<Levels>>
Vector reference_matvec(const Vector &rhs) {
    constexpr Index extent = logical_extent<Levels, Settings>();
    Vector result = Vector::Zero(extent * BlockSize);
    std::array<Index, Levels> row{}, col{};
    for (Index logical_row = 0; logical_row < extent; ++logical_row) {
        decode_logical_index<Levels, Settings>(logical_row, row);
        auto output = result.segment(logical_row * BlockSize, BlockSize);
        for (Index logical_col = 0; logical_col < extent; ++logical_col) {
            decode_logical_index<Levels, Settings>(logical_col, col);
            const auto [diagonal, constant] = block_coefficients<Levels>(row, col);
            const auto input = rhs.segment(logical_col * BlockSize, BlockSize);
            output.array() += diagonal * input.array() + constant * input.sum();
        }
    }
    return result;
}

template <int BlockSize, int Levels, typename Settings = LevelSettings<Levels>>
Vector make_rhs() {
    Vector result(logical_extent<Levels, Settings>() * BlockSize);
    for (Index i = 0; i < result.size(); ++i) {
        const double value = static_cast<double>(i + 1);
        result[i] = Scalar{std::sin(0.013 * value), std::cos(0.017 * value)};
    }
    return result;
}

template <int BlockSize, int Levels, typename Settings = LevelSettings<Levels>>
bool prepare_benchmark(benchmark::State &state, Vector &rhs, Vector &output, Vector &expected) {
    if (!configuration_fits_memory<BlockSize, Levels, Settings>()) {
        state.SkipWithError("Toeplitz coefficient storage or construction peak exceeds the memory guard");
        return false;
    }
    Eigen::setNbThreads(1);
    rhs = make_rhs<BlockSize, Levels, Settings>();
    output = Vector::Zero(rhs.size());
    expected = reference_matvec<BlockSize, Levels, Settings>(rhs);
    return true;
}

void record_error(benchmark::State &state, const Vector &output, const Vector &expected) {
    const double error = (output - expected).norm() / std::max(1.0, expected.norm());
    if (error > 2e-11)
        state.SkipWithError("matvec differs from the independent reference");
}

template <typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>>
void benchmark_old_classic(benchmark::State &state) {
    Vector rhs, output, expected;
    if (!prepare_benchmark<BlockSize, Levels, Settings>(state, rhs, output, expected))
        return;
    auto matrix = make_old_matrix<Block, BlockSize, Levels, Settings>();
    old_classic_multiply_add(matrix, rhs.data(), output.data());
    record_error(state, output, expected);
    output.setZero();
    for (auto _ : state) {
        old_classic_multiply_add(matrix, rhs.data(), output.data());
        benchmark::DoNotOptimize(output.data());
        benchmark::ClobberMemory();
    }
}

template <typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>>
void benchmark_old_batched(benchmark::State &state) {
    Vector rhs, output, expected;
    if (!prepare_benchmark<BlockSize, Levels, Settings>(state, rhs, output, expected))
        return;
    auto matrix = make_old_matrix<Block, BlockSize, Levels, Settings>();
    matrix.matvec_wise(rhs.data(), static_cast<std::size_t>(rhs.size()), output.data(),
                       static_cast<std::size_t>(output.size()));
    record_error(state, output, expected);
    output.setZero();
    for (auto _ : state) {
        matrix.matvec_wise(rhs.data(), static_cast<std::size_t>(rhs.size()), output.data(),
                           static_cast<std::size_t>(output.size()));
        benchmark::DoNotOptimize(output.data());
        benchmark::ClobberMemory();
    }
}

template <typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>>
void benchmark_nested_batched(benchmark::State &state) {
    Vector rhs, output, expected;
    if (!prepare_benchmark<BlockSize, Levels, Settings>(state, rhs, output, expected))
        return;
    auto matrix = make_new_matrix<Block, BlockSize, Levels, Settings>();
    typename decltype(matrix)::workspace_type workspace;
    matrix.matvec_add(rhs, output, workspace);
    record_error(state, output, expected);
    output.setZero();
    for (auto _ : state) {
        matrix.matvec_add(rhs, output, workspace);
        benchmark::DoNotOptimize(output.data());
        benchmark::ClobberMemory();
    }
}

template <typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>>
void benchmark_runtime_batched(benchmark::State &state) {
    Vector rhs, output, expected;
    if (!prepare_benchmark<BlockSize, Levels, Settings>(state, rhs, output, expected))
        return;
    auto matrix = make_runtime_matrix<Block, BlockSize, Levels, Settings>();
    typename decltype(matrix)::workspace_type workspace;
    matrix.matvec_add(rhs, output, workspace);
    record_error(state, output, expected);
    output.setZero();
    for (auto _ : state) {
        matrix.matvec_add(rhs, output, workspace);
        benchmark::DoNotOptimize(output.data());
        benchmark::ClobberMemory();
    }
}

template <typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>>
void benchmark_nested_prepared(benchmark::State &state) {
    Vector rhs, output, expected;
    if (!prepare_benchmark<BlockSize, Levels, Settings>(state, rhs, output, expected))
        return;
    auto matrix = make_new_matrix<Block, BlockSize, Levels, Settings>();
    typename decltype(matrix)::workspace_type workspace;

    // Both the initial gather and the final conversion used for validation are
    // outside the timed loop. matvec_prepared_add itself includes GEMM,
    // reduction to logical output blocks, and repacking into the persistent
    // panel layout required by the next prepared multiplication.
    auto prepared_rhs = matrix.prepare_vector(rhs);
    auto prepared_output = matrix.make_prepared_vector();
    matrix.matvec_prepared_add(prepared_rhs, prepared_output, workspace);
    output = matrix.unpack_vector(prepared_output);
    record_error(state, output, expected);

    // Verify that the prepared result can be fed back without conversion and
    // that the alias-safe in-place form preserves the same layout.
    auto chained = prepared_rhs;
    matrix.matvec_prepared(chained, chained, workspace);
    output = matrix.unpack_vector(chained);
    record_error(state, output, expected);
    matrix.matvec_prepared(chained, chained, workspace);
    record_error(state, matrix.unpack_vector(chained),
                 reference_matvec<BlockSize, Levels, Settings>(expected));
    prepared_output = matrix.make_prepared_vector();

    for (auto _ : state) {
        matrix.matvec_prepared_add(prepared_rhs, prepared_output, workspace);
        benchmark::DoNotOptimize(prepared_output.data());
        benchmark::ClobberMemory();
    }
}

template <typename Function>
void register_benchmark(const std::string &name, Function function) {
    benchmark::RegisterBenchmark(name.c_str(), function)
        ->Unit(benchmark::kMicrosecond)
        ->UseRealTime()
        ->Repetitions(independent_repetitions)
        ->ReportAggregatesOnly(true);
}

enum class Algorithm { old_classic, old_batched, nested_batched, runtime_batched, nested_prepared };

template <Algorithm Selected, typename Block, int BlockSize, int Levels,
          typename Settings = LevelSettings<Levels>>
void register_configuration(const char *block_name) {
    const std::string suffix = "/" + std::string(Settings::label) + "/B" + block_name;
    if constexpr (Selected == Algorithm::old_classic)
        register_benchmark("OldClassic" + suffix,
                           &benchmark_old_classic<Block, BlockSize, Levels, Settings>);
    else if constexpr (Selected == Algorithm::old_batched)
        register_benchmark("OldBatched" + suffix,
                           &benchmark_old_batched<Block, BlockSize, Levels, Settings>);
    else if constexpr (Selected == Algorithm::nested_batched)
        register_benchmark("NestedBatched" + suffix,
                           &benchmark_nested_batched<Block, BlockSize, Levels, Settings>);
    else if constexpr (Selected == Algorithm::runtime_batched)
        register_benchmark("RuntimeBatched" + suffix,
                           &benchmark_runtime_batched<Block, BlockSize, Levels, Settings>);
    else
        register_benchmark("NestedPrepared" + suffix,
                           &benchmark_nested_prepared<Block, BlockSize, Levels, Settings>);
}

template <Algorithm Selected, typename Block, int BlockSize>
void register_block_size(const char *block_name) {
    register_configuration<Selected, Block, BlockSize, 1>(block_name);
    register_configuration<Selected, Block, BlockSize, 2>(block_name);
    register_configuration<Selected, Block, BlockSize, 3>(block_name);
}

template <Algorithm Selected, typename Block, int BlockSize>
void register_five_cube_block_size(const char *block_name) {
    register_configuration<Selected, Block, BlockSize, 3, Level3FiveSettings>(block_name);
}

template <Algorithm Selected> void register_algorithm() {
    register_block_size<Selected, Block3, 3>("3");
    register_block_size<Selected, Block5, 5>("5");
    register_block_size<Selected, Block16, 16>("16");
    register_block_size<Selected, Block32, 32>("32");
    register_block_size<Selected, Block64, 64>("64");
    register_block_size<Selected, Block128, 128>("128");
    register_block_size<Selected, Block256, 256>("256");
    register_block_size<Selected, Block1024, 1024>("1024");

    register_five_cube_block_size<Selected, Block3, 3>("3");
    register_five_cube_block_size<Selected, Block5, 5>("5");
    register_five_cube_block_size<Selected, Block16, 16>("16");
    register_five_cube_block_size<Selected, Block32, 32>("32");
    register_five_cube_block_size<Selected, Block64, 64>("64");
    register_five_cube_block_size<Selected, Block128, 128>("128");
    register_five_cube_block_size<Selected, Block256, 256>("256");
}

[[maybe_unused]] const bool benchmarks_registered = [] {
    register_algorithm<Algorithm::old_classic>();
    register_algorithm<Algorithm::old_batched>();
    register_algorithm<Algorithm::nested_batched>();
    register_algorithm<Algorithm::runtime_batched>();
    register_algorithm<Algorithm::nested_prepared>();
    return true;
}();

class SummaryConsoleReporter final : public benchmark::ConsoleReporter {
  public:
    using benchmark::ConsoleReporter::ConsoleReporter;

    void ReportRuns(const std::vector<Run> &reports) override {
        const bool has_statistical_aggregates =
            std::any_of(reports.begin(), reports.end(), [](const Run &run) {
                return run.run_type == Run::RT_Aggregate &&
                       (run.aggregate_name == "mean" || run.aggregate_name == "stddev");
            });
        std::vector<Run> filtered;
        filtered.reserve(reports.size());
        for (const auto &run : reports) {
            const bool selected_aggregate =
                run.run_type == Run::RT_Aggregate &&
                (run.aggregate_name == "mean" || run.aggregate_name == "stddev");
            const bool single_default_run =
                !has_statistical_aggregates && run.run_type == Run::RT_Iteration;
            if (run.skipped != benchmark::internal::NotSkipped || selected_aggregate ||
                single_default_run)
                filtered.push_back(run);
        }
        benchmark::ConsoleReporter::ReportRuns(filtered);
    }
};

} // namespace

int main(int argc, char **argv) {
    benchmark::MaybeReenterWithoutASLR(argc, argv);
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv))
        return 1;
    SummaryConsoleReporter reporter(benchmark::ConsoleReporter::OO_Tabular);
    benchmark::RunSpecifiedBenchmarks(&reporter);
    benchmark::Shutdown();
    return 0;
}
