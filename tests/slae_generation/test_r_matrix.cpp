//
// Created by evgen on 25.11.2024.
//

#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "meshes/plate/PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "types/Types.hpp"
#include "gtest/gtest.h"

#include <mesh/Parser.hpp>

using namespace EMW;
using namespace EMW::Types;

class R_MATRIX : public testing::Test {
  protected:
    std::string nodes = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/examples/cylinder/"
                        "mesh/2002_nodes.csv";
    std::string cells = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/examples/cylinder/"
                        "mesh/2050_cells.csv";
    static constexpr EMW::Types::index nNodes = 2002;
    static constexpr EMW::Types::index nCells = 2050;
    scalar tolerance = 1e-11;
    complex_d k{1, 0};
    Mesh::SurfaceMesh cylinder_mesh;

    void SetUp() override {
        auto parser_output = EMW::Parser::parseMesh(nodes, cells, nNodes, nCells);
        cylinder_mesh = Mesh::SurfaceMesh{parser_output.first, parser_output.second};
    }
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

TEST_F(R_MATRIX, SIMPLE_COINSIDANCE_BETWEEN_2_IMPLEMENTATION) {
    int N = 20;
    scalar h = 1. / (N - 1);
    const auto mesh = Examples::Plate::generateRectangularMesh(N, N, h, h);

    const auto m = Matrix::getMatrixR(k, mesh);
    const auto new_m = Matrix::getMatrixR(k, mesh, mesh);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ASSERT_NEAR(std::abs(m(i, j) - new_m(i, j)), 0, tolerance);
        }
    }
}

TEST_F(R_MATRIX, COINSIDANCE_BETWEEN_2_IMPLEMENTATION_ON_CYLINDER) {
    const auto m = Matrix::getMatrixR(k, cylinder_mesh);
    const auto new_m = Matrix::getMatrixR(k, cylinder_mesh, cylinder_mesh);
    const MatrixXc diff = m - new_m;
    std::cout << diff.norm() << std::endl;
    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j <  m.rows(); j++) {
            ASSERT_NEAR(std::abs(m(i, j) - new_m(i, j)), 0, tolerance);
        }
    }
}
