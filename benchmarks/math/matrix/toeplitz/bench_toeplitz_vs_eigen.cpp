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

MatrixXc make_random_toeplitz_matrix(std::int64_t n) {
    std::mt19937_64 gen(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    MatrixXc A(n, n);
    for (Eigen::Index i = 0; i < n; ++i) {
        for (Eigen::Index j = 0; j < n; ++j) {
            if (i >= j) {
                A(i, j) = i - j;
            } else {
                A(i, j) = i - j;
            }
        }
    }
    return A;
}

std::pair<MatrixXc, MyMatrix> make_pair_of_matricies(size_t n_blocks, size_t block_size) {
    auto tp_mat = EMW::Math::LinAgl::Matrix::ZeroToeplitzBlock<Complex>(n_blocks, block_size);
    auto mat = make_random_toeplitz_matrix(n_blocks * block_size);

    for (size_t idx = 0; idx < n_blocks; ++idx) {
        for (size_t jdx = 0; jdx < n_blocks; ++jdx) {
            const size_t rows_start = idx * block_size;
            const size_t cols_start = jdx * block_size;
            tp_mat.get_block(idx, jdx) = mat.block(rows_start, cols_start, block_size, block_size);
        }
    }
    // const auto dense_tp_m = tp_mat.to_dense();
    // std::cout << "Difference in norms = " << (mat - dense_tp_m).norm() / mat.norm() << std::endl;
    return {mat, tp_mat};
}

VectorXc make_random_vector(std::int64_t n) {
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
        auto [m, tp] = make_pair_of_matricies(n, block_size_);
        A_ = std::move(m);
        my_matrix_ = std::move(tp);
        x_ = make_random_vector(n * block_size_);
        y_ = VectorXc::Zero(n * block_size_);
    }

protected:
    MatrixXc A_;
    VectorXc x_;
    VectorXc y_; // multiplication result;
    MyMatrix my_matrix_;
    size_t block_size_{};
};

BENCHMARK_DEFINE_F(MatrixVectorBenchmark, EigenMatVec)(benchmark::State &state) {
    for (auto _ : state) {
        y_.noalias() = A_ * x_;
        benchmark::DoNotOptimize(y_.data());
        benchmark::ClobberMemory();
    }

    const auto n = static_cast<double>(state.range(0));
    state.SetComplexityN(n);
}

BENCHMARK_DEFINE_F(MatrixVectorBenchmark, MyMatrixMatVec)(benchmark::State &state) {
    for (auto _ : state) {
        y_.noalias() = my_matrix_.matvec(x_);
        benchmark::DoNotOptimize(y_.data());
        benchmark::ClobberMemory();
    }

    const auto n = static_cast<double>(state.range(0));
    state.SetComplexityN(n);
}

BENCHMARK_DEFINE_F(MatrixVectorBenchmark, BlasZgemvFromEigen)(benchmark::State& state) {
    const Complex alpha{1.0, 0.0};
    const Complex beta{0.0, 0.0};
    const auto n = static_cast<int>(A_.rows());

    for (auto _ : state) {
        cblas_zgemv(
            CblasColMajor,
            CblasNoTrans,
            n,
            n,
            &alpha,
            A_.data(),
            n,
            x_.data(),
            1,
            &beta,
            y_.data(),
            1
        );

        benchmark::DoNotOptimize(y_.data());
        benchmark::ClobberMemory();
    }

    state.SetComplexityN(state.range(0));
}

BENCHMARK_REGISTER_F(MatrixVectorBenchmark, EigenMatVec)
->ArgsProduct({
    benchmark::CreateDenseRange(1, 5, 1),
    benchmark::CreateRange(32, 2048, 2)
})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_REGISTER_F(MatrixVectorBenchmark, MyMatrixMatVec)
->ArgsProduct({
    benchmark::CreateDenseRange(1, 5, 1),
    benchmark::CreateRange(32, 2048, 2)
})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_REGISTER_F(MatrixVectorBenchmark, BlasZgemvFromEigen)
->ArgsProduct({
    benchmark::CreateDenseRange(1, 5, 1),
    benchmark::CreateRange(32, 2048, 2)
})
    ->Unit(benchmark::kMillisecond);

} // namespace

BENCHMARK_MAIN();

// Комментарии: оптимальное ускорение получается в том случае, когда несколько
// тёплицевых блоков одновременно помещаются в кеш, то есть при размерах блоков 32 - 512