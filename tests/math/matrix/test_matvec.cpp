//
// Created by codex on 31.03.2026.
//

#include "gtest/gtest.h"
#include "math/matrix/Matrix.hpp"

#include <algorithm>

namespace {
using namespace EMW;

using complex_t = Types::complex_d;
using ToeplitzDense = Math::LinAgl::Matrix::ToeplitzBlock<complex_t>;
using Factor = Math::LinAgl::Matrix::DynamicFactoredMatrix<Types::MatrixX<complex_t>>;
using ToeplitzFactored = Math::LinAgl::Matrix::ToeplitzDynFactoredBlock<complex_t>;

Types::MatrixX<complex_t> make_dense_block(Types::index i, Types::index j, Types::index block_size)
{
    Types::MatrixX<complex_t> block(block_size, block_size);
    for (Types::index r = 0; r < block_size; ++r)
    {
        for (Types::index c = 0; c < block_size; ++c)
        {
            const double re = 0.01 * static_cast<double>((i + 1) * (r + 1) + (j + 2) * (c + 1));
            const double im = 0.001 * static_cast<double>((r == c ? 1 : -1) * (i + j + 1));
            block(r, c) = complex_t{re, im};
        }
    }
    return block;
}

Factor make_factored_block(Types::index i, Types::index j, Types::index block_size)
{
    auto left = make_dense_block(i, j, block_size);
    Types::MatrixX<complex_t> right = Types::MatrixX<complex_t>::Identity(block_size, block_size);
    for (Types::index k = 0; k < block_size; ++k)
    {
        right(k, k) = complex_t{1.0 + 0.05 * static_cast<double>(i + j + k + 1), 0.0};
    }
    return Factor{{std::move(left), std::move(right)}};
}

template<typename matrix_t>
void check_wise_equals_classic(matrix_t&& matrix, const Types::VectorX<complex_t>& x, double tol)
{
    const Types::VectorX<complex_t> expected = matrix.matvec(x);

    Types::VectorX<complex_t> x_mutable = x;
    Types::VectorX<complex_t> got = Types::VectorX<complex_t>::Zero(matrix.rows());
    matrix.matvec_wise(x_mutable.data(), static_cast<size_t>(x_mutable.size()),
                             got.data(), static_cast<size_t>(got.size()));

    const double rel_err = (expected - got).norm() / std::max(1.0, expected.norm());
    EXPECT_NEAR(rel_err, 0.0, tol);
}
} // namespace

TEST(MatvecWiseBlock, DenseToeplitzMatchesClassicMatvec)
{
    constexpr Types::index toeplitz_size = 6;
    constexpr Types::index block_size = 5;

    ToeplitzDense matrix(
        toeplitz_size, toeplitz_size,
        [](Types::index i, Types::index j) { return make_dense_block(i, j, block_size); });

    Types::VectorX<complex_t> x = Types::VectorX<complex_t>::Random(matrix.cols());
    check_wise_equals_classic(matrix, x, 1e-12);
}

TEST(MatvecWiseBlock, FactoredToeplitzMatchesClassicMatvec)
{
    constexpr Types::index toeplitz_size = 5;
    constexpr Types::index block_size = 4;

    ToeplitzFactored matrix(
        toeplitz_size, toeplitz_size,
        [](Types::index i, Types::index j) { return make_factored_block(i, j, block_size); });

    Types::VectorX<complex_t> x = Types::VectorX<complex_t>::Random(matrix.cols());
    check_wise_equals_classic(matrix, x, 1e-11);
}
