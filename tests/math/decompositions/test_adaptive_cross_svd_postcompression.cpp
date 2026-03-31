#include <gtest/gtest.h>

#include "math/matrix/decompositions/adaptive_cross.hpp"

using namespace EMW;

namespace {

using matrix_t = Types::MatrixXc;
using vector_t = Types::VectorXc;
using value_t = Types::complex_d;
using aca_t = Math::LinAgl::Decompositions::ACA<matrix_t, vector_t, value_t>;
using dynamic_factored_t = Math::LinAgl::Matrix::DynamicFactoredMatrix<matrix_t>;

matrix_t orthonormal_columns(Types::index rows, Types::index cols) {
    matrix_t random_matrix = matrix_t::Random(rows, cols);
    Eigen::HouseholderQR<matrix_t> qr(random_matrix);
    return qr.householderQ() * matrix_t::Identity(rows, cols);
}

dynamic_factored_t make_uv_transposed_factorization(Types::index rows, Types::index cols,
                                                    const Types::VectorXd& singular_values) {
    const Types::index rank = singular_values.size();
    const matrix_t Q_left = orthonormal_columns(rows, rank);
    const matrix_t Q_right = orthonormal_columns(cols, rank);

    const matrix_t U = Q_left;
    const matrix_t V_transposed = singular_values.cast<value_t>().asDiagonal() * Q_right.transpose();

    Containers::vector<matrix_t> factors;
    factors.reserve(2);
    factors.emplace_back(U);
    factors.emplace_back(V_transposed);
    return {std::move(factors), Containers::vector<bool>{false, false}};
}

Types::scalar relative_error(const matrix_t& reference, const matrix_t& approx) {
    return (reference - approx).norm() / reference.norm();
}

} // namespace

TEST(ADAPTIVE_CROSS_SVD_POSTCOMPRESSION, REDUCES_RANK_WITH_EXPECTED_TOLERANCE) {
    const Types::VectorXd singular_values =
        (Types::VectorXd(6) << 10.0, 1.0, 1e-2, 1e-5, 1e-7, 1e-9).finished();
    auto factorized = make_uv_transposed_factorization(40, 35, singular_values);
    const matrix_t dense_before = factorized.to_dense();

    const auto compressed = aca_t::svd_postcompression(std::move(factorized), 1e-4);
    const matrix_t dense_after = compressed.to_dense();

    EXPECT_EQ(compressed.template get<0>().cols(), 3);
    EXPECT_EQ(compressed.template get<1>().rows(), 3);
    EXPECT_LT(relative_error(dense_before, dense_after), 1e-3);
}

TEST(ADAPTIVE_CROSS_SVD_POSTCOMPRESSION, PRESERVES_MATRIX_FOR_TINY_TOLERANCE) {
    const Types::VectorXd singular_values =
        (Types::VectorXd(4) << 5.0, 2.0, 0.8, 0.3).finished();
    auto factorized = make_uv_transposed_factorization(25, 30, singular_values);
    const matrix_t dense_before = factorized.to_dense();

    const auto compressed = aca_t::svd_postcompression(std::move(factorized), 1e-14);
    const matrix_t dense_after = compressed.to_dense();

    EXPECT_EQ(compressed.template get<0>().cols(), 4);
    EXPECT_EQ(compressed.template get<1>().rows(), 4);
    EXPECT_LT(relative_error(dense_before, dense_after), 1e-12);
}

TEST(ADAPTIVE_CROSS_SVD_POSTCOMPRESSION, DECREASES_FACTOR_MEMORY_AFTER_COMPRESSION) {
    const Types::VectorXd singular_values =
        (Types::VectorXd(5) << 20.0, 0.6, 1e-3, 1e-6, 1e-9).finished();
    auto factorized = make_uv_transposed_factorization(60, 55, singular_values);
    const Types::scalar memory_before = factorized.memory_usage();
    const matrix_t dense_before = factorized.to_dense();

    const auto compressed = aca_t::svd_postcompression(std::move(factorized), 1e-3);
    const Types::scalar memory_after = compressed.memory_usage();
    const matrix_t dense_after = compressed.to_dense();

    EXPECT_LT(memory_after, memory_before);
    EXPECT_LT(relative_error(dense_before, dense_after), 2e-3);
}
