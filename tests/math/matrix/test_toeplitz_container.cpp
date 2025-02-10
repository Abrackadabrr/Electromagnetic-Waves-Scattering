//
// Created by evgen on 10.02.2025.
//

#include "toeplitz_matrix_tests.hpp"
#include "math/matrix/ToeplitzContainer.hpp"

auto gilbert_matrix(Types::index i, Types::index j)->Types::scalar { return 1. / (i + j + 1); };

TEST_F(TOEPLITZ_MATRIX_TESTS, EIGEN_COMPARISON) {
    constexpr Types::index N = 50;
    constexpr Types::index M = 60;
    const Math::LinAgl::Matrix::ToeplitzContainer<Types::scalar> toeplitz(50, 60, gilbert_matrix);
    // Проверка на правильность заполнения чисел
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            ASSERT_EQ(toeplitz(i, j), gilbert_matrix(i, j));
        }
    }
    // Проверка на правильность хранения чисел (свойство теплицевости)
    for (int i = 0; i < N + M - 1; i++) {
        const Types::scalar reference_element = (i >= N ? toeplitz(0, i) : toeplitz(0, i - N + 1));
        for (int row = 0; row < M; row++) {
            for (int col = 0; col < N; col++) {
                if (row - col == i) {}
                ASSERT_EQ(toeplitz(row, col), reference_element);
            }
        }
    }
}
