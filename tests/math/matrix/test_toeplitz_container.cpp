//
// Created by evgen on 10.02.2025.
//

#include "toeplitz_matrix_tests.hpp"
#include "math/matrix/ToeplitzContainer.hpp"

auto simple_toeplitz_matrix(Types::index i, Types::index j)->Types::scalar { return static_cast<integer>(i) - j; };

TEST_F(TOEPLITZ_MATRIX_TESTS, TOEPLITZ_CONTAINER_1) {
    constexpr Types::index N = 50;
    constexpr Types::index M = 60;
    const Math::LinAgl::Matrix::ToeplitzContainer<Types::scalar> toeplitz(N, M, simple_toeplitz_matrix);
    // Проверка на правильность заполнения чисел
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            ASSERT_EQ(toeplitz(i, j), simple_toeplitz_matrix(i, j));
        }
    }
}

TEST_F(TOEPLITZ_MATRIX_TESTS, TOEPLITZ_CONTAINER_2) {
    constexpr Types::index N = 5000;
    constexpr Types::index M = 23454;
    const Math::LinAgl::Matrix::ToeplitzContainer<Types::scalar> toeplitz(N, M, simple_toeplitz_matrix);
    // Проверка на свойство тёплицевости
    for (Types::integer i = -static_cast<integer>(M) + 1; i < N; i++) {
        // i -- это номер диагонали
        const scalar ref_el= (i + M - 1) < N ? toeplitz(i + M - 1, M - 1) : toeplitz(i + M - 1 - N, 0);
        for (integer c = M - 1; c >= 0; c--) {
            // c -- номер столбца
            integer r = c + i;  // r -- это номер строки
            if (r >= 0 && r < N) {
                ASSERT_EQ(toeplitz(r, c), ref_el);
            }
         }
    }
}
