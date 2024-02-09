//
// Created by evgen on 24.01.24.
//

#include "PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "Eigen/IterativeLinearSolvers"
#include <iostream>

using namespace EMW;
using namespace EMW::Types;

struct physicalConditions {
    // polarization
    Vector3d E0;
    // wave number
    complex_d k;
    // wave vector
    Vector3d k_vec;
};

int main() {
    int N = 81;
    scalar h = 0.0125;
    // сетка
    auto *mesh3 = new Mesh::SurfaceMesh{EMW::Examples::Plate::generatePlatePrimaryMesh(N, h)};

    mesh3->setName("mesh_" + std::to_string(N-1) + "_1");

    // поляризация
    physicalConditions physics{
            .E0 = Vector3d{0, 1, 0},
            .k = complex_d{1, 0},
            .k_vec = Vector3d{-1, 0, 0}
    };

    // на мелкой
    const MatrixXd *A3 = new MatrixXd{Matrix::getMatrix(physics.k, mesh3->getCells()).real()};
    const VectorXc *b3 = new VectorXc{Matrix::getRHS(physics.E0, physics.k_vec, mesh3->getCells())};
    const auto cg3 = new Eigen::ConjugateGradient<MatrixXd, Eigen::Lower | Eigen::Upper>{*A3};
    const auto *j = new VectorXd{cg3->solve(b3->real())};

    mesh3->fillJ(*j);

    VTK::test_snapshot(0, *mesh3,
                       "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plate/OXY/");

}
