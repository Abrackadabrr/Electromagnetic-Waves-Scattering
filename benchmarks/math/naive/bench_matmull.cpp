// bench_matmul.cpp
//
// Сравнение трех вариантов умножения комплексных матриц A (n x n) и B (n x k):
// 1) По столбцам средствами Eigen: C.col(j) = A * B.col(j)
// 2) Целиком через cblas_zgemm:      C = A * B
// 3) По столбцам через cblas_zgemv:   C.col(j) = A * B.col(j)
//
// Пример сборки:
//   c++ -O3 -DNDEBUG -march=native -std=c++20 bench_matmul.cpp \
//       -I /path/to/eigen \
//       -lbenchmark -lpthread -lopenblas
//
// Для честного сравнения лучше запускать в 1 поток:
//   export OPENBLAS_NUM_THREADS=1
//   export MKL_NUM_THREADS=1
//   export OMP_NUM_THREADS=1

#include <benchmark/benchmark.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <cblas.h>

#include <complex>
#include <cstdint>
#include <random>
#include <vector>

namespace {

using Scalar = std::complex<double>;
using MatrixXcd = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using Index = Eigen::Index;

template <class Derived>
void FillRandom(Eigen::MatrixBase<Derived>& x, std::mt19937_64& gen) {
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (Index j = 0; j < x.cols(); ++j) {
        for (Index i = 0; i < x.rows(); ++i) {
            x.derived()(i, j) = Scalar(dist(gen), dist(gen));
        }
    }
}

struct BenchmarkData {
    MatrixXcd A;
    MatrixXcd B;
    MatrixXcd C;

    BenchmarkData(Index n, Index k) : A(n, n), B(n, k), C(n, k) {
        std::mt19937_64 gen(0xC0FFEEULL + static_cast<std::uint64_t>(1315423911ULL * n + 2654435761ULL * k));
        FillRandom(A, gen);
        FillRandom(B, gen);
        C.setZero();
    }
};

void SetCounters(benchmark::State& state, Index n, Index k) {
    // Для комплексного GEMM обычно считают ~8 flops на один complex multiply-add.
    const double flops_per_call =
        8.0 * static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(k);

    state.counters["n"] = static_cast<double>(n);
    state.counters["k"] = static_cast<double>(k);
    state.counters["calls/s"] =
        benchmark::Counter(static_cast<double>(state.iterations()), benchmark::Counter::kIsRate);
    state.counters["flops/s"] =
        benchmark::Counter(state.iterations() * flops_per_call, benchmark::Counter::kIsRate);
}

static void BM_Eigen_Columnwise(benchmark::State& state) {
    openblas_set_num_threads(1);
    Eigen::setNbThreads(1);

    const Index n = static_cast<Index>(state.range(0));
    const Index k = static_cast<Index>(state.range(1));

    BenchmarkData data(n, k);

    for (auto _ : state) {
        for (Index j = 0; j < k; ++j) {
            data.C.col(j).noalias() = data.A * data.B.col(j);
        }
        benchmark::DoNotOptimize(data.C.data());
        benchmark::ClobberMemory();
    }

    SetCounters(state, n, k);
}

static void BM_CBLAS_ZGEMM(benchmark::State& state) {
    openblas_set_num_threads(1);
    Eigen::setNbThreads(1);

    const Index n = static_cast<Index>(state.range(0));
    const Index k = static_cast<Index>(state.range(1));

    BenchmarkData data(n, k);

    const Scalar alpha(1.0, 0.0);
    const Scalar beta(0.0, 0.0);

    for (auto _ : state) {
        cblas_zgemm(CblasColMajor,
                    CblasNoTrans,
                    CblasNoTrans,
                    static_cast<int>(n),
                    static_cast<int>(k),
                    static_cast<int>(n),
                    &alpha,
                    data.A.data(),
                    static_cast<int>(n),
                    data.B.data(),
                    static_cast<int>(n),
                    &beta,
                    data.C.data(),
                    static_cast<int>(n));

        benchmark::DoNotOptimize(data.C.data());
        benchmark::ClobberMemory();
    }

    SetCounters(state, n, k);
}

static void BM_CBLAS_ZGEMV_Columnwise(benchmark::State& state) {
    openblas_set_num_threads(1);
    Eigen::setNbThreads(1);

    const Index n = static_cast<Index>(state.range(0));
    const Index k = static_cast<Index>(state.range(1));

    BenchmarkData data(n, k);

    const Scalar alpha(1.0, 0.0);
    const Scalar beta(0.0, 0.0);

    for (auto _ : state) {
        for (Index j = 0; j < k; ++j) {
            cblas_zgemv(CblasColMajor,
                        CblasNoTrans,
                        static_cast<int>(n),
                        static_cast<int>(n),
                        &alpha,
                        data.A.data(),
                        static_cast<int>(n),
                        data.B.col(j).data(),
                        1,
                        &beta,
                        data.C.col(j).data(),
                        1);
        }

        benchmark::DoNotOptimize(data.C.data());
        benchmark::ClobberMemory();
    }

    SetCounters(state, n, k);
}

void CommonArguments(benchmark::Benchmark* b) {
    const std::vector<int> ns = {64, 128, 192, 256, 384, 512, 1024};
    const std::vector<int> ks = {1, 4, 8, 16, 32, 64, 100, 200, 300, 500, 1000};

    for (int n : ns) {
        for (int k : ks) {
            b->Args({n, k});
        }
    }
}

BENCHMARK(BM_Eigen_Columnwise)
    ->Apply(CommonArguments)
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_CBLAS_ZGEMM)
    ->Apply(CommonArguments)
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_CBLAS_ZGEMV_Columnwise)
    ->Apply(CommonArguments)
    ->Unit(benchmark::kMillisecond);

}  // namespace

BENCHMARK_MAIN();