//
// Created by evgen on 18.01.2025.
//

#include "gtest/gtest.h"

#include "math/MathConstants.hpp"
#include "math/fields/SurfaceVectorField.hpp"
#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "research/lattice/GeneralizedEquations.hpp"
#include "types/Types.hpp"
#include "experiment/PhysicalCondition.hpp"

using namespace EMW;

TEST(WAVEGIUDE, MATRIX_TEST) {
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/8000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/2000_cells.csv";
    const EMW::Types::index nNodes = 8000;
    const EMW::Types::index nCells = 2000;

    // собираем сетки
    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};

    constexpr Types::index N1 = 2;
    constexpr Types::index N2 = 2;
    constexpr Types::index N1_x_N2 = N1 * N2;

    const Geometry::PeriodicStructure<N1, N2> geometry{0.1, 0.1, mesh_base};

    // собираем матрицу для расчета волновода
    const Types::scalar a = 0.07;
    const Types::scalar freq = Math::Constants::c / 1e8;
    const Types::complex_d k{Physics::get_k_on_frquency(freq), 0};
    const Types::complex_d beta = std::sqrt(k * k - (EMW::Math::Constants::PI_square<Types::scalar>() / (a * a)));

    const auto matrix = GeneralizedEquations::getMatrix(geometry, a, k);
    std::cout << "matrix" << std::endl;
    // далее проверяем матрицу по блокам
    const Types::index block_size = matrix.size() / N1_x_N2;

    // 1) Проверка равенства всех диагональных блоков
    const Types::MatrixXc matrix_11 = matrix.block(0, 0, block_size, block_size);
    for (int i = 0; i < N1_x_N2; i++) {
        const Types::MatrixXc matrix_ii = matrix.block(i * block_size, i * block_size, block_size, block_size);
        ASSERT_NEAR((matrix_11 - matrix_ii).norm(), 0, 1e-12);
    }

    // 2) Проверка на теплицевость (руками)
    // A_big_11 = A_big_22
    const Types::MatrixXc matrix_big_11 = matrix.block(0, 0, 2 * block_size, 2 * block_size);
    const Types::MatrixXc matrix_big_22 = matrix.block(2 * block_size, 2 * block_size, 2 * block_size, 2 * block_size);
    ASSERT_NEAR((matrix_big_11 - matrix_big_22).norm(), 0, 1e-12);

    // 3) Провека внутри недиагональных блоков
    const Types::MatrixXc matrix_13 = matrix.block(0, 2 * block_size, block_size, block_size);
    const Types::MatrixXc matrix_24 = matrix.block(block_size, 2 * block_size, block_size, block_size);
    ASSERT_NEAR((matrix_big_11 - matrix_big_22).norm(), 0, 1e-12);

    const Types::MatrixXc matrix_31 = matrix.block(2 * block_size, 0, block_size, block_size);
    const Types::MatrixXc matrix_42 = matrix.block(3 * block_size, block_size, block_size, block_size);
    ASSERT_NEAR((matrix_big_11 - matrix_big_22).norm(), 0, 1e-12);
}
