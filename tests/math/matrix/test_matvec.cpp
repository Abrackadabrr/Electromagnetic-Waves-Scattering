//
// Created by codex on 31.03.2026.
//

#include "gtest/gtest.h"
#include "math/matrix/Matrix.hpp"

#include <algorithm>

using namespace EMW;

using complex_t = Types::complex_d;
using ToeplitzDense = Math::LinAgl::Matrix::ToeplitzBlock<complex_t>;
using Factor = Math::LinAgl::Matrix::DynamicFactoredMatrix<Types::MatrixX<complex_t>>;
using ToeplitzFactored = Math::LinAgl::Matrix::ToeplitzDynFactoredBlock<complex_t>;


decltype(auto) make_toeplitz_block(size_t block_size, size_t fl) {
    auto tp_mat = EMW::Math::LinAgl::Matrix::ZeroToeplitzBlock<Types::complex_d>(fl, block_size);

    for (size_t idx = 0; idx < fl; ++idx) {
        for (size_t jdx = 0; jdx < fl; ++jdx) {
            tp_mat.get_block(idx, jdx) = Types::MatrixXc::Random(block_size, block_size);
        }
    }
    return tp_mat;
}

decltype(auto) make_double_toeplitz_block(size_t block_size, size_t fl, size_t sl) {
    auto tp_mat = EMW::Math::LinAgl::Matrix::ZeroDoubleToeplitzBlock<Types::complex_d>(fl, sl, block_size);

    for (size_t idx_s = 0; idx_s < sl; ++idx_s) {
        for (size_t jdx_s = 0; jdx_s < sl; ++jdx_s) {
            for (size_t idx_f = 0; idx_f < fl; ++idx_f) {
                for (size_t jdx_f = 0; jdx_f < fl; ++jdx_f) {
                    tp_mat.get_block(idx_s, jdx_s).get_block(idx_f, idx_f) =
                        Types::MatrixXc::Random(block_size, block_size);
                }
            }
        }
    }
    return tp_mat;
}

template <typename matrix_t>
void check_wise_equals_classic(matrix_t &&matrix, const Types::VectorX<complex_t> &x, double tol) {
    const Types::VectorX<complex_t> expected = matrix.matvec(x);

    Types::VectorX<complex_t> x_mutable = x;
    Types::VectorX<complex_t> got = Types::VectorX<complex_t>::Zero(matrix.rows());
    matrix.matvec_wise(x_mutable.data(), static_cast<size_t>(x_mutable.size()),
                       got.data(), static_cast<size_t>(got.size()));

    const double rel_err = (expected - got).norm() / std::max(1.0, expected.norm());
    EXPECT_NEAR(rel_err, 0.0, tol);
}

TEST(MatvecWiseBlock, DenseToeplitzMatchesClassicMatvec) {
    constexpr Types::index toeplitz_size = 6;
    constexpr Types::index block_size = 5;

    ToeplitzDense matrix = make_toeplitz_block(block_size, toeplitz_size);
    Types::VectorX<complex_t> x = Types::VectorX<complex_t>::Random(matrix.cols());
    check_wise_equals_classic(matrix, x, 1e-12);
}

TEST(MatvecWiseBlock, DenseDoubleToeplitzMatchesClassicMatvec) {
    constexpr Types::index fl = 6;
    constexpr Types::index sl = 6;
    constexpr Types::index block_size = 5;

    auto matrix = make_double_toeplitz_block(block_size, fl, sl);
    Types::VectorX<complex_t> x = Types::VectorX<complex_t>::Random(matrix.cols());
    check_wise_equals_classic(matrix, x, 1e-12);
}
