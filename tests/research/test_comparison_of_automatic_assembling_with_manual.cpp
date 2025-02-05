//
// Created by evgen on 04.02.2025.
//

#include "gtest/gtest.h"

#include "research/lattice/Equations.hpp"
#include "research/lattice/GeneralizedEquations.hpp"

#include "lattice_tests.hpp"

/**
 * Тест на сравнение матриц, получаемых в результате разных по форме алгоритмов
 * расчета матрицы системы уравнений для двух рядом стоящих волноводов
 */
TEST_F(AssemblingTests, AutomaticAssembled_vs_ManuallyAssembled) {
    // Для начала создадим периодическую структуру,
    constexpr Types::index N1 = 2;
    constexpr Types::index N2 = 1;
    constexpr Types::index N1_x_N2 = N1 * N2;

    const Geometry::PeriodicStructure<N1, N2> geometry{0.1, 0.1, mesh_base};

    const auto& mesh_1 = geometry.get_mesh_matrix()[0][0];
    const auto& mesh_2 = geometry.get_mesh_matrix()[1][0];

    // Далее сделаем расчет двух матриц
    // 1) Ручной сбор
    const auto matrix1 = Lattice::getSLAE(mesh_1, mesh_2, a, k);
    // 2) Автоматический сбор
    const auto matrix2 = GeneralizedEquations::getMatrix(geometry, a, k);

    ASSERT_NEAR((matrix1 - matrix1).norm(), 0, 1e-14);
}
