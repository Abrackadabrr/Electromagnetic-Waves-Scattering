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
using MatrixXc = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using VectorXc = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
using MyMatrix = EMW::Math::LinAgl::Matrix::ToeplitzToeplitzBlock<Complex>;

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

std::pair<MatrixXc, MyMatrix> make_pair_of_matricies(size_t inner_size, size_t fl, size_t sl) {
    auto tp_mat = EMW::Math::LinAgl::Matrix::ZeroDoubleToeplitzBlock<Complex>(fl, sl, inner_size);
    auto mat = make_random_toeplitz_matrix(inner_size * fl * sl);

    for (size_t idx_f = 0; idx_f < fl; ++idx_f) {
        for (size_t jdx_f = 0; jdx_f < fl; ++jdx_f) {
            for (size_t idx_s = 0; idx_s < sl; ++idx_s) {
                for (size_t jdx_s = 0; jdx_s < sl; ++jdx_s) {
                    tp_mat.get_block(idx_s, jdx_s).get_block(idx_f, idx_f) = MatrixXc::Random(inner_size, inner_size);
                }
            }
        }
    }
    // в целом для бенчмарка нам наплевать какие именно данные лежат в матрице, главное -- чтобы был одинаковый размер
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
        auto [m, tp] = make_pair_of_matricies(block_size_, n, n);
        A_ = std::move(m);
        my_matrix_ = std::move(tp);
        x_ = make_random_vector(n * n * block_size_);
        y_ = VectorXc::Zero(n * n * block_size_);
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

BENCHMARK_DEFINE_F(MatrixVectorBenchmark, BlasZgemvFromEigen)(benchmark::State &state) {
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
    benchmark::CreateRange(32, 512, 2)
})
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oNSquared);

BENCHMARK_REGISTER_F(MatrixVectorBenchmark, MyMatrixMatVec)
->ArgsProduct({
    benchmark::CreateDenseRange(1, 5, 1),
    benchmark::CreateRange(32, 512, 2)
})
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oNSquared);

#ifdef BENCH_BLAS
BENCHMARK_REGISTER_F(MatrixVectorBenchmark, BlasZgemvFromEigen)
->ArgsProduct({
    benchmark::CreateDenseRange(1, 5, 1),
    benchmark::CreateRange(32, 512, 2)
})
    ->Unit(benchmark::kMillisecond)
    ->Complexity(benchmark::oNSquared);
#endif
} // namespace

BENCHMARK_MAIN();
