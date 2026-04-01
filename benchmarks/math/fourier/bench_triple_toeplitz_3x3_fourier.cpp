//
// Created by evgen on 31.03.2026.
//

#include <benchmark/benchmark.h>

#include "math/fourier/TripleToeplitz3x3Fourier.hpp"

#include "mesh/volume_mesh/CubeMeshWithData.hpp"

#include "operators/volume/OperatorK.hpp"

#include <Utils.hpp>

using namespace EMW;

namespace {
constexpr Types::scalar SPHERE_RADUIS = 1;

using toeplitz_t = Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d>;
using skeleton_t = Math::LinAgl::Matrix::TripleToeplitzFactoredBlock<Types::complex_d>;
using fourier_t = Math::Fourier::TripleToeplitz3x3Fourier<Types::complex_d>;

class TripleToeplitz3x3FourierBenchmark : public benchmark::Fixture {
public:
    void SetUp(const benchmark::State &state) override {
        Eigen::setNbThreads(1);

        const Types::index nodes_per_axis = static_cast<Types::index>(state.range(0));

        constexpr Types::scalar cube_length = 2.1 * SPHERE_RADUIS;
        const Types::scalar mesh_one_axis_size = cube_length / (nodes_per_axis - 1);
        const Types::scalar basis_fn_norm = 1. / (mesh_one_axis_size * std::sqrt(mesh_one_axis_size));

        auto mesh_ = Mesh::VolumeMesh::CubeMeshWithData(
            Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2},
            (nodes_per_axis - 1) * mesh_one_axis_size,
            (nodes_per_axis - 1) * mesh_one_axis_size,
            (nodes_per_axis - 1) * mesh_one_axis_size,
            nodes_per_axis, nodes_per_axis, nodes_per_axis);

        auto operator_k_ = Operators::Volume::operator_K_over_cube_mesh(
            Types::complex_d{6 * M_PI, 0.}, mesh_);

        toeplitz_ = operator_k_.compute_galerkin_matrix(basis_fn_norm);
        fourier_ = fourier_t(toeplitz_);

        std::cout << "Fourier matrx: " << std::endl;
        std::cout << Utils::get_memory_usage(toeplitz_) << std::endl;

        x_ = Types::VectorXc::Random(toeplitz_.cols());
        y_ = Types::VectorXc::Zero(toeplitz_.rows());
    }

protected:
    toeplitz_t toeplitz_;
    fourier_t fourier_;

    Types::VectorXc x_;
    Types::VectorXc y_;
};

class TripleToeplitz3x3SkeletonBenchmark : public benchmark::Fixture {
public:
    void SetUp(const benchmark::State &state) override {
        Eigen::setNbThreads(1);
        openblas_set_num_threads(1);

        const Types::index nodes_per_axis = static_cast<Types::index>(state.range(0));
        const Types::index block_side = static_cast<Types::index>(state.range(1));

        constexpr Types::scalar cube_length = 2.1 * SPHERE_RADUIS;
        const Types::scalar mesh_one_axis_size = cube_length / (nodes_per_axis - 1);
        const Types::scalar basis_fn_norm = 1. / (mesh_one_axis_size * std::sqrt(mesh_one_axis_size));

        auto mesh_ = Mesh::VolumeMesh::CubeMeshWithData(
            Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2},
            (nodes_per_axis - 1) * mesh_one_axis_size,
            (nodes_per_axis - 1) * mesh_one_axis_size,
            (nodes_per_axis - 1) * mesh_one_axis_size,
            nodes_per_axis, nodes_per_axis, nodes_per_axis);

        auto operator_k_ = Operators::Volume::operator_K_over_cube_mesh(
            Types::complex_d{6 * M_PI, 0.}, mesh_);

        auto [skeleton, perm] = operator_k_.compute_galerkin_matrix_custom_blocksize_compressed(
            block_side, block_side, block_side, basis_fn_norm, 1e-5);
        skeleton_ = std::move(skeleton);
        perm_ = std::move(perm);

        x_ = Types::VectorXc::Random(skeleton_.cols());
        x_permuted_ = perm_ * x_;
        y_ = Types::VectorXc::Zero(skeleton_.rows());

        std::cout << "N = " << nodes_per_axis << std::endl;
        std::cout << Utils::get_memory_usage(skeleton_) << std::endl;
        std::cout << "M-rank = " << Utils::get_elements_for_parametrization(skeleton_) / skeleton_.cols() << std::endl;
    }

protected:
    skeleton_t skeleton_;
    Eigen::PermutationMatrix<Eigen::Dynamic> perm_;

    Types::VectorXc x_;
    Types::VectorXc x_permuted_;
    Types::VectorXc y_;
};

BENCHMARK_DEFINE_F(TripleToeplitz3x3SkeletonBenchmark, SkeletonFormatMatVec)(benchmark::State &state) {
    for (auto _ : state) {
        y_ = skeleton_.matvec(x_permuted_);
        benchmark::DoNotOptimize(y_.data());
        benchmark::ClobberMemory();
    }
    // state.counters["unknowns"] = static_cast<double>(x_.size());
}

BENCHMARK_DEFINE_F(TripleToeplitz3x3FourierBenchmark, FourierMatVec)(benchmark::State &state) {
    for (auto _ : state) {
        y_ = fourier_.matvec(x_);
        benchmark::DoNotOptimize(y_.data());
        benchmark::ClobberMemory();
    }
    // state.counters["unknowns"] = static_cast<double>(x_.size());
}

BENCHMARK_REGISTER_F(TripleToeplitz3x3SkeletonBenchmark, SkeletonFormatMatVec)
->Args({41, 4})
        ->Unit(benchmark::kMillisecond);

BENCHMARK_REGISTER_F(TripleToeplitz3x3FourierBenchmark, FourierMatVec)
->Args({41, 4})
        ->Unit(benchmark::kMillisecond);
} // namespace

BENCHMARK_MAIN();
