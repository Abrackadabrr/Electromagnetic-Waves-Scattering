//
// Created by evgen on 06.07.24.
//
#include "mesh/Parser.hpp"
#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/MeshTypes.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "mesh/Parser.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/MathConstants.hpp"
#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "examples/Utils.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <iostream>
#include <fstream>
#include <ranges>

using namespace EMW;
using namespace EMW::Types;

int main() {
    const std::string nodesFile = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/8102_nodes.csv";
    const std::string cellsFile = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/8200_cells.csv";
    const EMW::Types::index nNodes = 8102;
    const EMW::Types::index nCells = 8200;
    auto * surfaceMesh = new Mesh::SurfaceMesh{EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells)};

    surfaceMesh->setName("surface_mesh__" + std::to_string(8102));

    // физика
    EMW::Physics::planeWaveCase physics(Vector3d{-1, 1, 0}.normalized(),     // polarization
                                        (13 + 1./3) * Math::Constants::PI<scalar>(),  // wave figure
                                        Vector3d{1, 1, 0}.normalized());   // wave unit vector
    std::cout << "lambda = " << 2 * Math::Constants::PI<scalar>() / physics.k << std::endl;
    // окружающая сетка
    int N_volume = 161;
    scalar h_volume = 0.04;
    Containers::vector<Mesh::point_t> nodes;
    nodes.reserve(N_volume * N_volume);
    Utils::cartesian_product(std::ranges::views::iota(0, N_volume), std::ranges::views::iota(0, N_volume),
                      std::ranges::views::iota(0, 1), N_volume, N_volume, 1, h_volume, h_volume, 1,
                      std::back_inserter(nodes));
    const auto cellView = nodes | std::views::transform([](const Mesh::point_t &p) { return Mesh::Node{p}; });
    Mesh::VolumeMesh volumeMesh{*surfaceMesh, {cellView.begin(), cellView.end()}};
    volumeMesh.setName("volume_mesh__" + std::to_string(N_volume));

    // расчет
    const VectorXc b = VectorXc{Matrix::getRHS(physics.E0, physics.k_vec, surfaceMesh->getCells())};
    const MatrixXc A = Matrix::getMatrix(physics.k, surfaceMesh->getCells());

    // preconditioning
    const VectorXc JacobiAux = A.diagonal().cwiseInverse();
    const VectorXc newRHS = A.diagonal().cwiseInverse().cwiseProduct(b);
    const MatrixXc newMatrix = JacobiAux.asDiagonal() * A;

    auto method = Eigen::GMRES<MatrixXc>{};
    method.setMaxIterations(20000);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(1e-3);
    method.set_restart(2000);
    method.compute(newMatrix);
    const auto j = VectorXc{method.solveWithGuess(newRHS, newRHS)};
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    surfaceMesh->fillJ(j);

    // считаем поле вокруг
    volumeMesh.calculateAll(physics.E0, physics.k_vec, physics.k);

    VTK::surface_snapshot(1, *surfaceMesh, Pathes::examples + "cylinder/");
    VTK::volume_snapshot(1, volumeMesh, Pathes::examples + "cylinder/");
}
