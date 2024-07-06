//
// Created by evgen on 30.01.24.
//

#include "examples/plate/PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/MathConstants.hpp"
#include "examples/pathes.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <iostream>

using namespace EMW;
using namespace EMW::Types;

int main() {
    int N1 = 41;
    scalar h1 = 1. / (N1 - 1);
    int N2 = 31;
    scalar h2 = 1. / (N2 - 1);

    // сетка
    auto *surfaceMesh = new Mesh::SurfaceMesh{EMW::Examples::Plate::generateRectangularMesh(N1, N2, h1, h2)};

    surfaceMesh->setName("surface_mesh_test_" + std::to_string(N1) + "_x_" + std::to_string(N2));

    VTK::surface_snapshot(1, *surfaceMesh, Pathes::examples + "plane/rectangular/");
}
