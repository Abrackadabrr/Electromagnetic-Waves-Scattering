//
// Created by evgen on 25.11.2024.
//

#include <vector>
#include <utility>
#include <tuple>
#include <ranges>

#include "gtest/gtest.h"
#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/MeshTypes.hpp"
#include "meshes/plate/PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"

using namespace EMW;
using namespace EMW::Types;

class R_MATRIX : public testing::Test {
protected:
    scalar tolerance = 1e-15;
    complex_d k{1, 0};
};

TEST_F(R_MATRIX, ZERO_VALUES_ON_A_PLANE) {
    int N = 40;
    scalar h = 1. / (N - 1);
    const auto mesh = Examples::Plate::generateRectangularMesh(N, N, h, h);

    const auto m = Matrix::getMatrixR(k, mesh);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // std::cout << m(i, j) << " ";
            ASSERT_NEAR(std::abs(m(i, j)), 0, tolerance);
        }
    }
}
