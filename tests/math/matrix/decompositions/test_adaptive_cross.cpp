//
// Created by evgen on 15.05.2025.
//

#include <gtest/gtest.h>

#include "math/matrix/decompositions/Decompositions.hpp"

#include "mat_decomp.hpp"

TEST_F(MATRIX_DECOMPOSITIONS_TESTS, ACA_AUX_TEST1) {
    const Types::index N = 1000;
    const Types::index M = 1200;
    const auto mat = get_sinus_mat(N, M);
    const auto function_for_element = [&mat](const Types::index i, const Types::index j) { return mat(i, j); };

    // Совпадение стоблцов
    for (int i = 0; i < M; i++) {
        const auto&& col = Math::LinAgl::Decompositions::RealACA::get_col(i, N, M, function_for_element);
        ASSERT_NEAR((col - mat.col(i)).norm(), 0, 1e-14);
    }

    // Совпадение строк
    for (int i = 0; i < N; i++) {
        const auto&& row = Math::LinAgl::Decompositions::RealACA::get_row(i, N, M, function_for_element);
        ASSERT_NEAR((row - mat.row(i).transpose()).norm(), 0, 1e-14);
    }

    // 1. Проверка расчета строчек и столбцов от произведения матриц
    const Types::MatrixXd UVT = mat * (mat.transpose());  // матрицы N x N

    // Совпадение стоблцов (всего в итоговой матрице N стоблцов)
    for (int i = N; i < N; i++) {
        const auto&& col = Math::LinAgl::Decompositions::RealACA::get_col_UV(mat, mat, i);
        ASSERT_NEAR((col - UVT.col(i)).norm() / col.norm(), 0, 1e-14);
    }

    // Совпадение строк
    for (int i = 0; i < N; i++) {
        const auto&& row =  Math::LinAgl::Decompositions::RealACA::get_row_UV(mat, mat, i);
        ASSERT_NEAR((row - UVT.row(i).transpose()).norm() / row.norm(), 0, 1e-14);
    }

    // 2. Вычисление нормы матрицы


}

TEST_F(MATRIX_DECOMPOSITIONS_TESTS, ACA_SINUS) {
    const Types::index N = 3;
    const auto mat = get_sinus_mat(N);

    const auto function_for_element = [&mat](const Types::index i, const Types::index j) { return mat(i, j); };

    const Types::scalar error = 1e-6;
    const auto factored_matrix = Math::LinAgl::Decompositions::RealACA::compute(function_for_element, N, N, error);

    const Types::scalar err = (factored_matrix.compute() - mat).norm() / mat.norm();
    std::cout << err << std::endl;
}
