//
// Created by evgen on 10.02.2025.
//

#ifndef TOEPLITZ_MATRIX_TESTS_HPP
#define TOEPLITZ_MATRIX_TESTS_HPP

#include "meshes/plate/PlateGrid.hpp"
#include "gtest/gtest.h"

using namespace EMW;
using namespace EMW::Types;

class TOEPLITZ_MATRIX_TESTS : public testing::Test {
  protected:
    Types::scalar static simple_toeplitz_matrix(Types::index i, Types::index j) {
        return static_cast<scalar>(i) - static_cast<scalar>(j);
    };

    Types::MatrixX<complex_d> static get_block(Types::index i, Types::index j, Types::index size) {
        MatrixX<scalar> block(size, size);
        for (Types::index k = 0; k < size; ++k)
            for (Types::index l = 0; l < size; ++l)
            block(k, l) = static_cast<scalar>(k + l)/100.0;

        return complex_d{1, 1} * block * (simple_toeplitz_matrix(i, j) + (i == j ? 100 : 0));
    }

    Types::MatrixX<complex_d> get_toeplitz_matrix(Types::index tl_row, Types::index tl_cols, Types::index block_size) {
        Types::MatrixX<complex_d> result = Types::MatrixX<scalar>::Zero(tl_row * block_size, tl_cols * block_size);
#pragma omp parallel for collapse(2) num_threads(14)
        for (Types::index i = 0; i < tl_row; i++)
            for (Types::index j = 0; j < tl_cols; j++)
                result.block(block_size * i, block_size * j, block_size, block_size) = get_block(i, j, block_size);
        return result;
    }

    template <typename matrix_t>
    void final_check_for_vectors(const Types::MatrixX<complex_d> &matrix, const matrix_t &toeplitz,
                                 const Types::VectorX<complex_d> &vec) {
        auto start = std::chrono::system_clock::now();

        const Types::VectorX<complex_d> result_straightforward = matrix * vec;

        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Прямое произведение: " << elapsed.count() << '\n';

        start = std::chrono::system_clock::now();

        const Types::VectorX<complex_d> result_special = toeplitz * vec;

        end = std::chrono::system_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Произведение с новой матрицей: " << elapsed.count() << '\n';

        // Сравниваем результаты умножения
        ASSERT_NEAR((result_straightforward - result_special).norm() / result_straightforward.norm(), 0, 1e-14);
    }
};

#endif // TOEPLITZ_MATRIX_TESTS_HPP
