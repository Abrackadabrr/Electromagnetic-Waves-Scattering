//
// Created by evgen on 11.02.2025.
//

#include "math/matrix/Toeplitz.hpp"
#include "math/matrix/TwiceToeplitzMatrix.hpp"

#include "toeplitz_matrix_tests.hpp"
#include "gtest/gtest.h"

auto simple_toeplitz_matrix(Types::index i, Types::index j)->Types::scalar { return static_cast<integer>(i) - static_cast<integer>(j); };

Types::MatrixX<scalar> get_block(Types::index i, Types::index j) {
    Types::MatrixX<scalar> result = Types::MatrixX<scalar>::Zero(1, 1);
    result(0, 0) = simple_toeplitz_matrix(i, j);
    return result;
}

Types::MatrixX<scalar> get_toeplitz_matrix(Types::index rows, Types::index cols) {
    Types::MatrixX<scalar> result = Types::MatrixX<scalar>::Zero(rows, cols);
    for (Types::index i = 0; i < rows; i++)
        for (Types::index j = 0; j < cols; j++)
            result(i, j) = simple_toeplitz_matrix(i, j);
    return result;
}

TEST_F(TOEPLITZ_MATRIX_TESTS, BLOCK_TOEPLITZ_MULTIPLICATION_1) {
    Types::index N = 1500;
    Types::index M = 1400;
    const auto toeplitz_matrix = get_toeplitz_matrix(N, M);
    const Math::LinAgl::Matrix::BlockToeplitz<scalar> matrix(N, M, get_block);

    Types::VectorX<scalar> vec = Types::VectorX<scalar>::Zero(M);
    for (Types::index i = 0; i < M; i++)
        vec(i) = static_cast<scalar>(i);
    ASSERT_NEAR((toeplitz_matrix * vec - matrix * vec).norm(), 0, 0);
}

TEST_F(TOEPLITZ_MATRIX_TESTS, BLOCK_TOEPLITZ_MULTIPLICATION_2) {
    Types::index N = 1400;
    Types::index M = 1500;
    const auto toeplitz_matrix = get_toeplitz_matrix(N, M);
    const Math::LinAgl::Matrix::BlockToeplitz<scalar> matrix(N, M, get_block);

    Types::VectorX<scalar> vec = Types::VectorX<scalar>::Zero(M);
    for (Types::index i = 0; i < M; i++)
        vec(i) = static_cast<scalar>(i);
    ASSERT_NEAR((toeplitz_matrix * vec - matrix * vec).norm(), 0, 0);
}
