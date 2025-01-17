//
// Created by evgen on 24.01.24.
//

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
    // сетка
//    Mesh::SurfaceMesh *mesh1 = new Mesh::SurfaceMesh{EMW::Examples::Plate::generatePlatePrimaryMesh(9, 0.27/2)};
//    Mesh::SurfaceMesh *mesh2 = new Mesh::SurfaceMesh{EMW::Examples::Plate::generatePlatePrimaryMesh(25, 0.09/2)};
    Mesh::SurfaceMesh *mesh3 = new Mesh::SurfaceMesh{EMW::Examples::Plate::generatePlatePrimaryMesh(73, 0.03 / 2)};

//    mesh1->setName("mesh_8");
//    mesh2->setName("mesh_24");
    mesh3->setName("mesh_121");

    // поляризация
    physicalConditions phisycs{
            .E0 = Vector3d{0, 1, 0},
            .k = complex_d{1, 0},
            .k_vec = Vector3d{-1, 0, 0}
    };

    // волновое число
    const complex_d k{1, 0};

//    // расчет на грубой сетке
//    const MatrixXd A1 = Matrix::getMatrix(k, mesh1->getCells()).real();
//    const VectorXc b1 = Matrix::getRHS(E0, k, mesh1->getCells());
//    const Eigen::ConjugateGradient<MatrixXd, Eigen::Lower | Eigen::Upper> cg1{A1};
//    j[0] = cg1.solve(b1.real());

//    // на средней
//    const MatrixXd A2 = Matrix::getMatrix(k, mesh2->getCells()).real();
//    const VectorXc b2 = Matrix::getRHS(E0, k, mesh2->getCells());
//    const Eigen::ConjugateGradient<MatrixXd, Eigen::Lower | Eigen::Upper> cg2{A2};
//    j[1] = cg2.solve(b2.real());
    // на мелкой
    const MatrixXd *A3 = new MatrixXd{Matrix::getMatrix(k, mesh3->getCells()).real()};
    const VectorXc *b3 = new VectorXc{Matrix::getRHS(phisycs.E0, phisycs.k_vec, mesh3->getCells())};
    const auto cg3 = new Eigen::ConjugateGradient<MatrixXd, Eigen::Lower | Eigen::Upper>{*A3};
    const auto *j = new VectorXd{cg3->solve(b3->real())};

//    mesh1->fillJ(j[0]);
//    mesh2->fillJ(j[1]);
    mesh3->fillJ(*j);
//
//    VTK::test_snapshot(0, *mesh1,
//                       "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plate/");
//
//    VTK::test_snapshot(0, *mesh2,
//                       "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plate/");

    VTK::surface_snapshot(0, *mesh3,
                          "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plate/");

}
