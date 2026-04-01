//
// Created by evgen on 31.03.2026.
//

#include <benchmark/benchmark.h>


#include <Eigen/Dense>
#include <complex>
#include <memory>
#include <random>

#include "math/matrix/Matrix.hpp"

#include <mesh/volume_mesh/CubeMeshWithData.hpp>
#include <operators/volume/OperatorK.hpp>

using Complex = std::complex<double>;
using MatrixXc = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXc = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;

using singleToep = EMW::Math::LinAgl::Matrix::ToeplitzBlock<Complex>;
using doubleToep = EMW::Math::LinAgl::Matrix::ToeplitzToeplitzBlock<Complex>;
using tripleToep = EMW::Math::LinAgl::Matrix::TripleToeplitzBlock<Complex>;

using singleToepFactor = EMW::Math::LinAgl::Matrix::ToeplitzDynFactoredBlock<Complex>;
using doubleToepFactor = EMW::Math::LinAgl::Matrix::ToeplitzToeplitzDynFactoredBlock<Complex>;
using tripleToepFactor = EMW::Math::LinAgl::Matrix::TripleToeplitzFactoredBlock<Complex>;


using namespace EMW;

class MatrixVectorBenchmark : public benchmark::Fixture {
public:
    void SetUp(const benchmark::State &state) override {
    openblas_set_num_threads(1);
        size_t Nx = state.range(0);
        size_t Ny = state.range(1);
        size_t Nz = state.range(2);
        size_t nx = state.range(3);
        size_t ny = state.range(4);
        size_t nz = state.range(5);

        // Сетка на кубе
        constexpr Types::scalar cube_length = 2.1;
        const Types::scalar mesh_one_axis_size = cube_length / (Nx - 1);
        Mesh::VolumeMesh::CubeMeshWithData mesh{Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2},
                                                (Nx - 1) * mesh_one_axis_size,
                                                (Ny - 1) * mesh_one_axis_size, (Nz - 1) * mesh_one_axis_size, Nx, Ny, Nz};
        const Types::scalar cube_measure = mesh.dx() * mesh.dy() * mesh.dz();
        const Types::scalar basis_fn_module = 1. / sqrt(cube_measure);

        // Считаем плотную и сжатую матрицы одновременно
        Types::complex_d k{2 * M_PI, 0};
        Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};
        auto mat_full = Math::LinAgl::Matrix::ZeroTripleToeplitzBlock<Types::complex_d>(1, 1, 1, 1);
        auto [mat_compressed, perm] = operator_K.
            compute_galerkin_matrix_custom_blocksize_compressed(nx, ny, nz, 1, 1e-5, mat_full);

        // Собираем матрицы разных форматов
        m_dense_3 = std::move(mat_full);
        m_compressed_3 = std::move(mat_compressed);

        m_compressed_2 = m_compressed_3.get_block(0, 0);
        m_compressed_1 = m_compressed_2.get_block(0, 0);

        m_dense_2 = m_dense_3.get_block(0, 0);
        m_dense_1 = m_dense_2.get_block(0, 0);

        // Случайные вектора, на которые умножаем
        x_3 = VectorXc::Random(m_compressed_3.cols());
        x_2 = VectorXc::Random(m_dense_2.cols());
        x_1 = VectorXc::Random(m_dense_1.cols());
        // Результаты умножения
        y_3 = VectorXc::Zero(m_compressed_3.cols());
        y_2 = VectorXc::Zero(m_dense_2.cols());
        y_1 = VectorXc::Zero(m_dense_1.cols());
    }

protected:
    VectorXc x_3;
    VectorXc x_2;
    VectorXc x_1;
    VectorXc y_3; // multiplication result;
    VectorXc y_2; // multiplication result;
    VectorXc y_1; // multiplication result;
    tripleToepFactor m_compressed_3;
    doubleToepFactor m_compressed_2;
    singleToepFactor m_compressed_1;

    tripleToep m_dense_3;
    doubleToep m_dense_2;
    singleToep m_dense_1;
};

BENCHMARK_DEFINE_F(MatrixVectorBenchmark, InplaceMatvec)(benchmark::State &state) {
    for (auto _ : state) {
        m_compressed_3.matvec(x_3, y_3);
        benchmark::DoNotOptimize(y_3.data());
        benchmark::ClobberMemory();
    }

    const auto n = static_cast<double>(state.range(0));
    state.SetComplexityN(n);
}

BENCHMARK_DEFINE_F(MatrixVectorBenchmark, RawDataMatvec)(benchmark::State &state) {
    openblas_set_num_threads(1);
    for (auto _ : state) {
        m_compressed_3.matvec_wise(x_3.data(), x_3.size(), y_3.data(), y_3.size());
        benchmark::DoNotOptimize(y_3.data());
        benchmark::ClobberMemory();
    }
    const auto n = static_cast<double>(state.range(0));
    state.SetComplexityN(n);
};

BENCHMARK_REGISTER_F(MatrixVectorBenchmark, InplaceMatvec)
->Args({161, 21, 13, 4, 4, 4})
->Unit({benchmark::kMillisecond});

BENCHMARK_REGISTER_F(MatrixVectorBenchmark, RawDataMatvec)
->Args({161, 21, 13, 4, 4, 4})
->Unit({benchmark::kMillisecond});

BENCHMARK_MAIN();
