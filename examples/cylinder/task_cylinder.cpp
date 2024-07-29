//
// Created by evgen on 20.02.24.
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

using namespace EMW;
using namespace EMW::Types;

struct physicalConditions {
    // polarization
    Vector3d E0;
    // wavelenght
    scalar lambda;
    // wave number
    scalar k;
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
    const std::string nodes = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/8102_nodes.csv";
    const std::string cells = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/8200_cells.csv";
    const EMW::Types::index nNodes = 8102;
    const EMW::Types::index nCells = 8200;
    auto * surfaceMesh = new Mesh::SurfaceMesh{EMW::Parser::parseMesh(nodes, cells, nNodes, nCells)};

    surfaceMesh->setName("surface_mesh_" + std::to_string(8102));

    // физика
    const scalar frec_div_c = 8 / 0.299792458;
    const scalar omega_div_c = frec_div_c * (2 * Math::Constants::PI<scalar>());
    const scalar k = omega_div_c;
    const Vector3d E0 = Vector3d{0, 0, 1}.normalized();
    std::cout << "Wavelenght = " << (2 * Math::Constants::PI<scalar>()) / k << std::endl;


    // матрица
    const MatrixXc A3 = Matrix::getMatrix(k, surfaceMesh->getCells());
    // preconditioning
    const VectorXc JacobiAux = A3.diagonal().cwiseInverse();
    const MatrixXc newMatrix = JacobiAux.asDiagonal() * A3;

    // цикл по векторам правой части
    int samples = 180;
    Containers::vector_d esas;
    esas.reserve(samples);
    Containers::vector_d angles;
    angles.reserve(samples);

    for (int i = 0; i < samples; i++) {
    	// surfaceMesh copy
        // угол в радианах
        Types::scalar angle = i * Math::Constants::PI<scalar>() * 2 / samples;
        // вектор k
        Types::Vector3d k_vec = {-std::cos(angle), -std::sin(angle), 0};
        // расчет обратного эпр
        const VectorXc b3 = VectorXc{Matrix::getRHS(E0, k * k_vec, surfaceMesh->getCells())};
        const VectorXc newRHS = JacobiAux.cwiseProduct(b3);

        std::cout << "Matrix " << i << " has been built" << std::endl;

        auto method = Eigen::GMRES<MatrixXc>{};
        method.setMaxIterations(20000);
        std::cout << method.maxIterations() << std::endl;
        method.setTolerance(5e-5);
        method.set_restart(150);
        method.compute(A3);
        const auto j = VectorXc{method.solveWithGuess(b3, newRHS)};
        std::cout << "total iterations: " << method.iterations() << std::endl;
        std::cout << "total error: " << method.error() << std::endl;
        std::cout << "Info: " << static_cast<int>(method.info()) << std::endl << std::endl;
        surfaceMesh->fillJ(j);
        esas.push_back(ESA::calculateESA(-k_vec, k, *surfaceMesh));
        angles.push_back(angle);
    }

    std::ofstream sigmaFile(
            "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/cylinder/sigma_back_" +
            std::to_string(2002) + ".csv");

    to_csv(esas, angles, sigmaFile);
}
