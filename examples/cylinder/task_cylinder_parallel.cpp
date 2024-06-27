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
#include "experiment/ESA.hpp"
#include "math/MathConstants.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <iostream>
#include <fstream>
#include <omp.h>

using namespace EMW;
using namespace EMW::Types;

struct physicalConditions {
    // polarization
    Vector3d E0;
    // wavelenght
    scalar lambda;
    // wave number
    complex_d k;
    // wave vector
    Vector3d k_vec;
};

template<template<typename> typename Container, typename T>
void to_csv(const Container<T> &cont1, const Container<T> &cont2, std::ofstream &str) {
    str << "sigma,alpha\n";
    for (int i = 0; i < cont1.size(); i++) {
        str << cont1[i] << ',' << cont2[i] << '\n';
    }
}

int main() {
    std::atomic<int> tasks_done = 0;
    const std::string nodes = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/2002_nodes.csv";
    const std::string cells = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/2050_cells.csv";
    const EMW::Types::index nNodes = 2002;
    const EMW::Types::index nCells = 2050;
    auto * surfaceMesh = new Mesh::SurfaceMesh{EMW::Parser::parseMesh(nodes, cells, nNodes, nCells)};

    surfaceMesh->setName("surface_mesh_" + std::to_string(2002));

    // физика
    const scalar frec_div_c = 2 / 0.299792458;
    const scalar omega_div_c = frec_div_c * (2 * Math::Constants::PI<scalar>());
    const scalar k = omega_div_c;
    std::cout << "Wavelenght = " << (2 * Math::Constants::PI<scalar>()) / k << std::endl;


    // матрица
    const MatrixXc A3 = Matrix::getMatrix(k, surfaceMesh->getCells());
    // preconditioning
    const VectorXc JacobiAux = A3.diagonal().cwiseInverse();
    const MatrixXc newMatrix = JacobiAux.asDiagonal() * A3;

    // цикл по векторам правой части
    int samples = 10;
    Containers::vector_d esas;
    esas.resize(samples);
    Containers::vector_d angles;
    angles.resize(samples);
#pragma omp parallel for num_threads(7)
    for (int i = 0; i < samples; i++) {
        // surfaceMesh copy
        Mesh::SurfaceMesh meshCopy = *surfaceMesh;
        // угол в радианах
        Types::scalar angle = i * Math::Constants::PI<scalar>() / samples;
        // вектор k
        Types::Vector3d k_vec = {-std::cos(angle), -std::sin(angle), 0};
        // вектор E
        const Vector3d E0 = Vector3d{std::sin(angle), -std::cos(angle), 0}.normalized();
        // расчет обратного эпр
        const VectorXc b3 = VectorXc{Matrix::getRHS(E0, k * k_vec, meshCopy.getCells())};
        const VectorXc newRHS = JacobiAux.cwiseProduct(b3);

        int rank = omp_get_thread_num();
        std::cout << "Matrix " << i << " has been built by thread #" << rank << std::endl;

        auto method = Eigen::GMRES<MatrixXc>{};
        method.setMaxIterations(10000);
        // std::cout << method.maxIterations() << std::endl;
        method.setTolerance(1e-7);
        method.set_restart(5000);
        method.compute(A3);
        const auto j = VectorXc{method.solveWithGuess(b3, newRHS)};
        std::cout << "total error: " << method.error() << std::endl;
        std::cout << "Info: " << static_cast<int>(method.info()) << std::endl << std::endl;
        meshCopy.fillJ(j);
        esas[i] = (ESA::calculateESA(-k_vec, k, meshCopy));
        angles[i] = (angle);
        tasks_done++;
        std::cout << "Tasks done: " << tasks_done << std::endl;
    }

    std::ofstream sigmaFile(
            "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/cylinder/sigma_back_" +
            std::to_string(2002) + "_2GG_small_mesh_29_may_18_28.csv");

    to_csv(esas, angles, sigmaFile);
}
