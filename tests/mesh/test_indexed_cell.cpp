//
// Created by evgen on 14.02.24.
//
#include "gtest/gtest.h"

#include "operators/Operators.hpp"
#include "operators/Functions.hpp"
#include "types/Types.hpp"
#include "mesh/MeshTypes.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "integration/newton_cotess/Rectangular.hpp"

using namespace EMW;

class CELL_TESTS : public testing::Test {};

TEST_F(CELL_TESTS, PRAMETRIZATION) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::vector<Mesh::Point> points = {
            Mesh::Point{0, 0, 0},
            Mesh::Point{1, 0, 0},
            Mesh::Point{1, 1, 0},
            Mesh::Point{0, 1, 0}
    };
    Mesh::IndexedCell cell{{0, 1, 2, 3}, points};

    Types::scalar h = 0.01;
    for (int i = 0; i <= 100; i++)
        for (int j = 0; j <= 100; j++) {
            const Mesh::Point paramPoint = cell.parametrization(i * h, j * h);
            const Types::scalar error = (Types::Vector3d{i * h, j * h, 0} - paramPoint).norm();
            ASSERT_NEAR(error, 0, 1e-14) << i << ' ' << j;
            const Types::scalar mul = cell.multiplier(i * h, j * h);
            ASSERT_NEAR(mul, 1, 1e-14);
        }
}
