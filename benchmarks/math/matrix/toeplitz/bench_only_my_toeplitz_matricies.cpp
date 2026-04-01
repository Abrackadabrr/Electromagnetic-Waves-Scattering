//
// Created by evgen on 27.03.2026.
//

#include <benchmark/benchmark.h>

#include <Eigen/Dense>
#include <complex>
#include <cstdint>
#include <memory>
#include <random>

#include <cblas.h>

#include "math/matrix/Matrix.hpp"

using Complex = std::complex<double>;
using MatrixXc = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXc = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
using MyMatrix = EMW::Math::LinAgl::Matrix::ToeplitzBlock<Complex>;

namespace {

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

    for (size_t idx_f = 0; idx_f < fl; ++idx_f) {
        for (size_t jdx_f = 0; jdx_f < fl; ++jdx_f) {
            for (size_t idx_s = 0; idx_s < sl; ++idx_s) {
                for (size_t jdx_s = 0; jdx_s < sl; ++jdx_s) {
                    tp_mat.get_block(idx_s, jdx_s).get_block(idx_f, idx_f) =
                        MatrixXc::Random(block_size, block_size);
                }
            }
        }
    }
    return tp_mat;
}

VectorXc make_random_vector(size_t n) {
    std::mt19937_64 gen(1337);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    VectorXc x(n);
    for (Eigen::Index i = 0; i < x.size(); ++i) {
        x(i) = Complex(dist(gen), dist(gen));
    }
    return x;
}


class MatrixVectorBenchmark : public benchmark::Fixture {
public:
    void SetUp(const benchmark::State &state) override {
        openblas_set_num_threads(1);
        const auto n = static_cast<std::int64_t>(state.range(0));
        block_size_ = state.range(1);
        auto tp = make_toeplitz_block(block_size_, n);
        my_matrix_ = std::move(tp);
        x_ = make_random_vector(my_matrix_.rows());
        y_ = VectorXc::Zero(my_matrix_.rows());

    }

protected:
    VectorXc x_;
    VectorXc y_; // multiplication result;
    MyMatrix my_matrix_;
    size_t block_size_{};
};

BENCHMARK_DEFINE_F(MatrixVectorBenchmark, MyMatrixMatVecNew)(benchmark::State &state) {
    for (auto _ : state) {
        my_matrix_.matvec(x_, y_);
        benchmark::DoNotOptimize(y_.data());
        benchmark::ClobberMemory();
    }

    const auto n = static_cast<double>(state.range(0));
    state.SetComplexityN(n);
}

BENCHMARK_DEFINE_F(MatrixVectorBenchmark, MyMatrixMatVecWithoutCopy)(benchmark::State &state) {
    for (auto _ : state) {
        my_matrix_.matvec(x_.data(), x_.size(), y_);
        benchmark::DoNotOptimize(y_.data());
        benchmark::ClobberMemory();
    }

    const auto n = static_cast<double>(state.range(0));
    state.SetComplexityN(n);
}


BENCHMARK_DEFINE_F(MatrixVectorBenchmark, MyMatrixMatVecOld)(benchmark::State &state) {
    for (auto _ : state) {
        y_.noalias() = my_matrix_.matvec(x_);
        benchmark::DoNotOptimize(y_.data());
        benchmark::ClobberMemory();
    }

    const auto n = static_cast<double>(state.range(0));
    state.SetComplexityN(n);
}


BENCHMARK_REGISTER_F(MatrixVectorBenchmark, MyMatrixMatVecNew)
->ArgsProduct({
                                                                  benchmark::CreateDenseRange(100, 200, 10),
                                                                  benchmark::CreateRange(3, 3, 2)
                                                              })
    ->Unit(benchmark::kMillisecond);


BENCHMARK_REGISTER_F(MatrixVectorBenchmark, MyMatrixMatVecOld)
->ArgsProduct({
                                                                  benchmark::CreateDenseRange(100, 200, 10),
                                                                  benchmark::CreateRange(3, 3, 2)})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_REGISTER_F(MatrixVectorBenchmark, MyMatrixMatVecWithoutCopy)
->ArgsProduct({
                                                                          benchmark::CreateDenseRange(100, 200, 10),
                                                                          benchmark::CreateRange(3, 3, 2)})
    ->Unit(benchmark::kMillisecond);

} // namespace

BENCHMARK_MAIN();
