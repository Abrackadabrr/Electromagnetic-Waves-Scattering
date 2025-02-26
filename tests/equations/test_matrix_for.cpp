//
// Created by evgen on 06.02.2025.
//
#include "gtest/gtest.h"

#include "equations/EquationsOverGeometry.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "geometry/PeriodicStructure.hpp"
#include "math/MathConstants.hpp"
#include "math/fields/SurfaceVectorField.hpp"
#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "research/lattice/GeneralizedEquations.hpp"
#include "types/Types.hpp"

#include "tests/research/lattice_tests.hpp"

using namespace EMW;

constexpr Types::index N1 = 2;
constexpr Types::index N2 = 2;
constexpr Types::index N1_x_N2 = N1 * N2;

using slae_matrix = Equations::MatrixFor<Geometry::PeriodicStructure<N1, N2>, WaveGuideWithActiveSection::diagonal_t,
                                         WaveGuideWithActiveSection::submatrix_t>;

TEST(AutomaticalAssembling, Generalized_MATRIX_TEST) {
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/tests/"
                                  "meshes_for_tests/lattice_redused/320_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/tests/"
                                  "meshes_for_tests/lattice_redused/80_cells.csv";
    const EMW::Types::index nNodes = 320;
    const EMW::Types::index nCells = 80;

    // собираем сетки
    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};

    const Geometry::PeriodicStructure<N1, N2> geometry{0.1, 0.1, mesh_base};

    // собираем матрицу для расчета волновода
    const auto matrix = slae_matrix::compute(geometry, WaveGuideWithActiveSection::diagonal,
                                             WaveGuideWithActiveSection::submatrix, {a, k});

    std::cout << "matrix" << std::endl;
    // далее проверяем матрицу по блокам
    const Types::index block_size = matrix.rows() / N1_x_N2;

    ASSERT_EQ(block_size, 320);

    // 1) Проверка равенства всех диагональных блоков
    const Types::MatrixXc &matrix_11 = matrix.block(0, 0, block_size, block_size);
    for (int i = 0; i < N1_x_N2; i++) {
        const Types::MatrixXc &matrix_ii = matrix.block(i * block_size, i * block_size, block_size, block_size);
        ASSERT_NEAR((matrix_11 - matrix_ii).norm(), 0, 1e-12);
    }

    // 2) Полная проверка на дважды теплицевость (руками)
    // A_big_11 = A_big_22
    const Types::MatrixXc &matrix_big_11 = matrix.block(0, 0, 2 * block_size, 2 * block_size);
    const Types::MatrixXc &matrix_big_22 = matrix.block(2 * block_size, 2 * block_size, 2 * block_size, 2 * block_size);
    ASSERT_NEAR((matrix_big_11 - matrix_big_22).norm(), 0, 1e-12);

    // 3) Проверка внутри недиагональных блоков
    const Types::MatrixXc &matrix_13 = matrix.block(0, 2 * block_size, block_size, block_size);
    const Types::MatrixXc &matrix_24 = matrix.block(block_size, 3 * block_size, block_size, block_size);
    ASSERT_NEAR((matrix_13 - matrix_24).norm(), 0, 1e-12);

    const Types::MatrixXc &matrix_31 = matrix.block(2 * block_size, 0, block_size, block_size);
    const Types::MatrixXc &matrix_42 = matrix.block(3 * block_size, block_size, block_size, block_size);
    ASSERT_NEAR((matrix_31 - matrix_42).norm(), 0, 1e-12);
}

TEST_F(AssemblingTests, GENERAL_REAL_MATRIX_TEST) {

    // Задаем геометрию
    const Geometry::PeriodicStructure<N1, N2> geometry{0.1, 0.1, mesh_base};

    // собираем матрицу для расчета волновода
    const auto matrix = slae_matrix::compute(geometry, WaveGuideWithActiveSection::diagonal,
                                             WaveGuideWithActiveSection::submatrix, {a, k});
    std::cout << "matrix" << std::endl;
    // далее проверяем матрицу по блокам
    const Types::index block_size = matrix.rows() / N1_x_N2;

    ASSERT_EQ(block_size, 4400);

    // 1) Проверка равенства всех диагональных блоков
    const Types::MatrixXc &matrix_11 = matrix.block(0, 0, block_size, block_size);
    for (int i = 0; i < N1_x_N2; i++) {
        const Types::MatrixXc &matrix_ii = matrix.block(i * block_size, i * block_size, block_size, block_size);
        ASSERT_NEAR((matrix_11 - matrix_ii).norm(), 0, 1e-12);
    }

    // 2) Проверка на теплицевость (руками)
    // A_big_11 = A_big_22
    const Types::MatrixXc &matrix_big_11 = matrix.block(0, 0, 2 * block_size, 2 * block_size);
    const Types::MatrixXc &matrix_big_22 = matrix.block(2 * block_size, 2 * block_size, 2 * block_size, 2 * block_size);
    ASSERT_NEAR((matrix_big_11 - matrix_big_22).norm(), 0, 1e-12);

    // 3) Проверка внутри недиагональных блоков
    const Types::MatrixXc &matrix_13 = matrix.block(0, 2 * block_size, block_size, block_size);
    const Types::MatrixXc &matrix_24 = matrix.block(block_size, 3 * block_size, block_size, block_size);
    ASSERT_NEAR((matrix_13 - matrix_24).norm(), 0, 1e-12);

    const Types::MatrixXc &matrix_31 = matrix.block(2 * block_size, 0, block_size, block_size);
    const Types::MatrixXc& matrix_42 = matrix.block(3 * block_size, block_size, block_size, block_size);
    ASSERT_NEAR((matrix_31 - matrix_42).norm(), 0, 1e-12);
}