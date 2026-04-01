//
// Created by evgen on 31.03.2026.
//

#include <benchmark/benchmark.h>


#include <Eigen/Dense>
#include <complex>
#include <cstdint>
#include <memory>
#include <random>

#include "math/matrix/Matrix.hpp"

using Complex = std::complex<double>;
using MatrixXc = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXc = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;

using singleToep = EMW::Math::LinAgl::Matrix::ToeplitzBlock<Complex>;
using doubleToep = EMW::Math::LinAgl::Matrix::ToeplitzToeplitzBlock<Complex>;
using tripleToep = EMW::Math::LinAgl::Matrix::TripleToeplitzBlock<Complex>;

using singleToepFactor = EMW::Math::LinAgl::Matrix::ToeplitzDynFactoredBlock<Complex>;
using doubleToepFactor = EMW::Math::LinAgl::Matrix::ToeplitzToeplitzDynFactoredBlock<Complex>;
using tripleToepFactor = EMW::Math::LinAgl::Matrix::TripleToeplitzFactoredBlock<Complex>;


decltype(auto) make_toeplitz_block(size_t block_size, size_t fl) {
    auto tp_mat = EMW::Math::LinAgl::Matrix::ZeroToeplitzBlock<Complex>(fl, block_size);

    for (size_t idx = 0; idx < fl; ++idx) {
        for (size_t jdx = 0; jdx < fl; ++jdx) {
            tp_mat.get_block(idx, jdx) = MatrixXc::Random(block_size, block_size);
        }
    }
    return tp_mat;
}

decltype(auto) make_double_toeplitz_block(size_t block_size, size_t fl, size_t sl) {
    auto tp_mat = EMW::Math::LinAgl::Matrix::ZeroDoubleToeplitzBlock<Complex>(block_size, fl, sl);

    for (size_t idx_s = 0; idx_s < sl; ++idx_s) {
        for (size_t jdx_s = 0; jdx_s < sl; ++jdx_s) {
            for (size_t idx_f = 0; idx_f < fl; ++idx_f) {
                for (size_t jdx_f = 0; jdx_f < fl; ++jdx_f) {
                    tp_mat.get_block(idx_s, jdx_s).get_block(idx_f, idx_f) =
                        MatrixXc::Random(block_size, block_size);
                }
            }
        }
    }
    return tp_mat;
}

using namespace EMW::Types;

class MatrixVectorBenchmark : public benchmark::Fixture {
public:
    void SetUp(const benchmark::State &state) override {
    openblas_set_num_threads(1);
        size_t inner_size = state.range(0);
        size_t fl = state.range(1);
        size_t sl = state.range(2);

        m_dense_2 = make_double_toeplitz_block(inner_size, fl, sl);
        m_dense_1 = make_toeplitz_block(inner_size, fl);

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
        m_dense_2.matvec(x_2, y_2);
        benchmark::DoNotOptimize(y_2.data());
        benchmark::ClobberMemory();
    }

    const auto n = static_cast<double>(state.range(0));
    state.SetComplexityN(n);
}

BENCHMARK_DEFINE_F(MatrixVectorBenchmark, RawDataMatvec)(benchmark::State &state) {
    openblas_set_num_threads(1);
    for (auto _ : state) {
        m_dense_2.matvec_wise(x_2.data(), x_2.size(), y_2.data(), y_2.size());
        benchmark::DoNotOptimize(y_2.data());
        benchmark::ClobberMemory();
    }
    const auto n = static_cast<double>(state.range(0));
    state.SetComplexityN(n);
};

BENCHMARK_REGISTER_F(MatrixVectorBenchmark, InplaceMatvec)
->Args({4 * 4 * 4 * 3, 10, 3})
->Args({4 * 4 * 4 * 3, 20, 3})
->Args({256, 40, 4})
->Args({512, 40, 2})
->Unit({benchmark::kMillisecond});

BENCHMARK_REGISTER_F(MatrixVectorBenchmark, RawDataMatvec)
->Args({4 * 4 * 4 * 3, 10, 3})
->Args({4 * 4 * 4 * 3, 20, 3})
->Args({256, 40, 4})
->Args({512, 40, 2})
->Unit({benchmark::kMillisecond});

BENCHMARK_MAIN();
