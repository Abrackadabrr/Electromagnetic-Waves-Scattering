//
// Created by evgen on 22.02.24.
//

#include "mesh/Parser.hpp"
#include "types/Types.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/MeshTypes.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "mesh/Parser.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/MathConstants.hpp"
#include "examples/pathes.hpp"
#include "experiment/ESA.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <iostream>
#include <fstream>

#include <omp.h>

using namespace EMW;
using namespace EMW::Types;

template<template<typename> typename Container, typename T>
void to_csv(const Container<T> &cont1, const Container<T> &cont2, std::ofstream &str) {
    str << "sigma,alpha\n";
    for (int i = 0; i < cont1.size(); i++) {
        str << cont1[i] << ',' << cont2[i] << '\n';
    }
}

int main() {

    Eigen::initParallel();

    const std::string nodesFile = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/8102_nodes.csv";
    const std::string cellsFile = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/8200_cells.csv";
    const EMW::Types::index nNodes = 8102;
    const EMW::Types::index nCells = 8200;
    Mesh::SurfaceMesh surfaceMesh = Mesh::SurfaceMesh{EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells)};
    surfaceMesh.setName("cylinder_" + std::to_string(nNodes));

    // физика
    const scalar frequency = 2;  // GHz
    // частота разделяется на 10, поскольку геометрия у нас в дециметрах
    const scalar frec_div_c = (frequency / 0.299792458) / 10;
    const scalar omega_div_c = frec_div_c * (2 * Math::Constants::PI<scalar>());
    const scalar k = omega_div_c;
    std::cout << "Wavelenght = " << (2 * Math::Constants::PI<scalar>()) / k << std::endl;

    // сбор матрицы
    const MatrixXc A = Matrix::getMatrix(k, surfaceMesh.getCells());
    // preconditioning
    const VectorXc JacobiAux = A.diagonal().cwiseInverse();
    const MatrixXc A_prec = JacobiAux.asDiagonal() * A;

    std::cout << "Matrix has been built" << std::endl;

    // цикл по векторам правой части
    int samples = 120;
    Containers::vector_d esas;
    esas.resize(samples);
    Containers::vector_d angles;
    angles.resize(samples);

    int tasks_done = 0;

#pragma omp parallel for schedule(dynamic) num_threads(7) firstprivate(surfaceMesh)
    for (int i = 0; i < samples; i++) {
        // угол в радианах
        Types::scalar angle = i * 2 * Math::Constants::PI<scalar>() / (samples);
        // вектор k
        Types::Vector3d k_vec = {-std::cos(angle), -std::sin(angle), 0};
        // вектор E
        const Vector3d E0_V{0, 0, 1};
        const Vector3d E0_H = Vector3d{std::sin(angle), -std::cos(angle), 0}.normalized();
        // правая часть
        const VectorXc b = VectorXc{Matrix::getRHS(E0_V, k * k_vec, surfaceMesh.getCells())};
        // предобуславливание
        const VectorXc b_prec = JacobiAux.cwiseProduct(b);

        int rank = omp_get_thread_num();

        auto method = Eigen::GMRES<MatrixXc>{};
        method.setMaxIterations(10000);
        method.setTolerance(1e-3);
        method.set_restart(1000);
        method.compute(A_prec);
        const auto j = VectorXc{method.solveWithGuess(b_prec, b_prec)};
        std::cout << "#" + std::to_string(rank) + " total error: " << method.error() << std::endl;
        std::cout << "#" + std::to_string(rank) + " total iterations: " << method.iterations() << std::endl;
        std::cout << "#" + std::to_string(rank) + " Info: " << static_cast<int>(method.info()) << std::endl;
        surfaceMesh.fillJ(j);
        esas[i] = (ESA::calculateESA(-k_vec, k, surfaceMesh));
        angles[i] = (angle);
        tasks_done++;
        std::cout << "Tasks done: " << tasks_done << std::endl << std::endl;
    } 

    std::ofstream sigmaFile(
            "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/cylinder/sigma_back_" +
            std::to_string(nCells) + "_2GG_V_big_mesh_1_em3.csv");

    to_csv(esas, angles, sigmaFile);
}
