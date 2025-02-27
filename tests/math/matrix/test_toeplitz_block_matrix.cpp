//
// Created by evgen on 11.02.2025.
//

#include "math/matrix/Matrix.hpp"

#include "toeplitz_matrix_tests.hpp"
#include "gtest/gtest.h"

Types::scalar simple_toeplitz_matrix(Types::index i, Types::index j) {
    return static_cast<integer>(i) - static_cast<integer>(j);
};

Types::MatrixX<scalar> get_block(Types::index i, Types::index j) {
    return Types::MatrixX<scalar>{{{1, 2}, {2, 3}}} * simple_toeplitz_matrix(i, j);
}

Types::MatrixX<scalar> get_toeplitz_matrix(Types::index rows, Types::index cols) {
    Types::MatrixX<scalar> result = Types::MatrixX<scalar>::Zero(2 * rows, 2 * cols);
    for (Types::index i = 0; i < rows; i++)
        for (Types::index j = 0; j < cols; j++)
            result.block(2 * i, 2 * j, 2, 2) = get_block(i, j);
    return result;
}

TEST_F(TOEPLITZ_MATRIX_TESTS, TOEPLITZ_BLOCK_MATVEC_1) {
    Types::index N = 1500;
    Types::index M = 1400;
    const Types::MatrixX<scalar> toeplitz_matrix = get_toeplitz_matrix(N, M);
    const Math::LinAgl::Matrix::ToeplitzBlock<scalar> matrix(N, M, get_block);

    Types::VectorX<scalar> vec = Types::VectorX<scalar>::Zero(2 * M);
    for (Types::index i = 0; i < 2 * M; i++)
        vec(i) = static_cast<scalar>(i);
    ASSERT_NEAR((toeplitz_matrix * vec - matrix * vec).norm(), 0, 0);
}

TEST_F(TOEPLITZ_MATRIX_TESTS, TOEPLITZ_BLOCK_MULL_BY_VALUE) {
    Types::index N = 1500;
    Types::index M = 1400;
    Types::scalar value = 10;
    Math::LinAgl::Matrix::ToeplitzBlock<scalar> matrix(N, M, get_block);
    Types::VectorX<scalar> vec = Types::VectorX<scalar>::Zero(M);
    for (Types::index i = 0; i < M; i++)
        vec(i) = static_cast<scalar>(i);
    const auto res_1 = matrix * vec;
    matrix *= value;
    const auto res_2 = matrix * vec;

    ASSERT_NEAR((value * res_1 - res_2).norm(), 0, 0);
}

using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzBlock<scalar>;
using block_t = Math::LinAgl::Matrix::ToeplitzBlock<scalar>::block_type;
using toeplitz_block_t = Math::LinAgl::Matrix::ToeplitzToeplitzBlock<block_t>::block_type;

TEST_F(TOEPLITZ_MATRIX_TESTS, TOEPLITZ_BLOCK_CONSTRUCT_EASILY) {
    // Описываем структуру матрицы и создаем первую
    Types::index N = 1500;
    Types::index M = 1400;
    Types::scalar value = 10;
    ToeplitzBlock matrix(N, M, get_block);

    // Создаем вторую матрицу
    Containers::vector<block_t> blocks;
    Types::index size = ToeplitzBlock::get_size_of_container(N, M); blocks.reserve(size);
    for (int i = 0; i < M; i++) {
        blocks.emplace_back(get_block(0, i));
    }
    for (int i = 1; i < N; i++) {
        blocks.emplace_back(get_block(i, 0));
    }
    ToeplitzBlock matrix_2(N, M, std::move(blocks));
    ASSERT_EQ(blocks.size(), 0);

    // Создаем вектор и умножаемся
    Types::VectorX<scalar> vec = Types::VectorX<scalar>::Zero(2 * M);
    for (Types::index i = 0; i < M; i++)
        vec(i) = static_cast<scalar>(i);

    // Собираем вектор и умножаем
    ASSERT_NEAR((matrix * vec - matrix_2 * vec).norm(), 0, 0);
}
