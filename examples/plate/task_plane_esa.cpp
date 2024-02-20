//
// Created by evgen on 15.02.24.
//

#include "PlateGrid.hpp"
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

template<typename Range1, typename Range2, typename OutputIterator>
void cartesian_productXY(Range1 const &r1, Range2 const &r2, OutputIterator out, Types::index N,
                         Types::scalar h) {
    using std::begin;
    using std::end;

    for (auto i = begin(r1); i != end(r1); ++i) {
        for (auto j = begin(r2); j != end(r2); ++j) {

            *out++ = Types::Vector3d{static_cast<Types::scalar>(*i) - static_cast<Types::scalar>(N / 2.),
                                     static_cast<Types::scalar>(*j) - static_cast<Types::scalar>(N / 2.),
                                     0} * h;

        }
    }
}

int main() {
    int N = 41;
    scalar h = 1. / (N-1);
    // сетка
    auto *surfaceMesh = new Mesh::SurfaceMesh{EMW::Examples::Plate::generatePlatePrimaryMesh(N, h)};

    surfaceMesh->setName("surface_mesh_" + std::to_string(N));

    // физика
    physicalConditions physics{
            .E0 = Vector3d{0, 0, 1}.normalized(),
            .k = complex_d{8 * Math::Constants::PI<scalar>(), 0},
            .k_vec = Vector3d{-1, 1, 0}.normalized()
    };

    // на мелкой
    const VectorXc b3 = VectorXc{Matrix::getRHS(physics.E0, physics.k.real() * physics.k_vec, surfaceMesh->getCells())};
    const MatrixXc A3 = Matrix::getMatrix(physics.k, surfaceMesh->getCells());
    
    auto method = Eigen::GMRES<MatrixXc>{};
    method.setMaxIterations(20000);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(1e-10);
    method.set_restart(1000);
    method.compute(A3);
    const auto j = VectorXc{method.solve(b3)};
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    surfaceMesh->fillJ(j);

    int samples = 3600;
    Containers::vector_d esas;
    esas.reserve(samples);
    Containers::vector_d angles;
    angles.reserve(samples);

    for (int i = 0; i < samples; i++) {
        Types::scalar angle = i * Math::Constants::PI<scalar>() * 2 / samples;
        Types::Vector3d tau = {std::cos(angle), std::sin(angle), 0};
        esas.push_back(ESA::calculateESA(tau, physics.k, *surfaceMesh));
        angles.push_back(angle);
    }
    std::ofstream sigma(
            "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plane/sigmas/sigma_" +
            std::to_string(N) + ".csv");

    to_csv(esas, angles, sigma);
}
