//
// Created by evgen on 11.02.2025.
//

#include "math/matrix/Matrix.hpp"

#include "toeplitz_matrix_tests.hpp"
#include "gtest/gtest.h"

TEST_F(TOEPLITZ_MATRIX_TESTS, TOEPLITZ_BLOCK_MATVEC_BENCHMARK) {
    Types::index tl_rows = 5000;
    Types::index tl_cols = 5000;
    Types::index block_size = 2;
    Types::index N = block_size * tl_rows;
    Types::index M = block_size * tl_cols;
    const Types::MatrixX<complex_d> toeplitz_matrix = get_toeplitz_matrix(tl_rows, tl_cols, block_size);
    const Math::LinAgl::Matrix::ToeplitzBlock<complex_d> matrix(
        tl_rows, tl_cols, [block_size](Types::index i, Types::index j) { return get_block(i, j, block_size); });

    VectorXc vec(N);
    vec.setRandom();

    final_check_for_vectors(toeplitz_matrix, matrix, vec);
}

TEST_F(TOEPLITZ_MATRIX_TESTS, TOEPLITZ_BLOCK_MULL_BY_VALUE) {
    Types::index tl_rows = 10;
    Types::index tl_cols = 10;
    Types::index block_size = 10;
    Types::index N = block_size * tl_rows;
    Types::index M = block_size * tl_cols;
    Math::LinAgl::Matrix::ToeplitzBlock<complex_d> matrix(
        tl_rows, tl_cols, [block_size, this](Types::index i, Types::index j) { return get_block(i, j, block_size); });

    VectorXc vec(M);
    vec.setRandom();

    const auto res_1 = matrix * vec;
    const complex_d value{10., 0};
    matrix *= value;
    const auto res_2 = matrix * vec;

    ASSERT_NEAR((value * res_1 - res_2).norm() / res_1.norm(), 0, 1e-14);
}

using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzBlock<complex_d>;
using block_t = ToeplitzBlock::block_type;

TEST_F(TOEPLITZ_MATRIX_TESTS, TOEPLITZ_BLOCK_CONSTRUCT_EASILY) {
    // Описываем структуру матрицы и создаем первую
    Types::index tl_rows = 10;
    Types::index tl_cols = 10;
    Types::index block_size = 10;
    Types::index N = block_size * tl_rows;
    Types::index M = block_size * tl_cols;
    const Math::LinAgl::Matrix::ToeplitzBlock<complex_d> matrix(
        tl_rows, tl_cols, [block_size, this](Types::index i, Types::index j) { return get_block(i, j, block_size); });

    // Создаем вторую матрицу
    Containers::vector<block_t> blocks;
    Types::index size = ToeplitzBlock::get_size_of_container(tl_rows, tl_cols);
    blocks.reserve(size);
    for (int i = 0; i < tl_rows; i++) {
        blocks.emplace_back(get_block(0, i, block_size));
    }
    for (int i = 1; i < tl_cols; i++) {
        blocks.emplace_back(get_block(i, 0, block_size));
    }
    const ToeplitzBlock matrix_2(tl_rows, tl_cols, std::move(blocks));
    ASSERT_EQ(blocks.size(), 0);

    // Создаем вектор и умножаемся
    VectorXc vec(2 * N);
    vec.setRandom();
    // Собираем вектор и умножаем
    ASSERT_NEAR((matrix * vec - matrix_2 * vec).norm(), 0, 0);
}
