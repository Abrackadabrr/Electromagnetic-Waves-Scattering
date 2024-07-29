//
// Created by evgen on 08.02.24.
//

#include "PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/MathConstants.hpp"
#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <iostream>

using namespace EMW;
using namespace EMW::Types;

struct physicalConditions {
    // polarization
    Vector3d E0;
    // wave figure in media without conductivity (vacuum of air)
    scalar k;
    // wave vector
    Vector3d k_vec;
};

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
    int N_volume = 41;
    scalar h_volume = 0.075;

    int N = 41;
    scalar h = 1. / (N - 1);
    // сетка
    auto *surfaceMesh = new Mesh::SurfaceMesh{EMW::Examples::Plate::generatePlatePrimaryMesh(N, h)};

    surfaceMesh->setName("quadratic_random_basises" + std::to_string(N));

    // физика
    EMW::Physics::planeWaveCase physics{Vector3d{0, 1, 0}.normalized(),
                                        4 * Math::Constants::PI<scalar>(),
                                        Vector3d{0, 0, 1}.normalized()};

    // на мелкой
    const VectorXc b = VectorXc{Matrix::getRHS(physics.E0, physics.k_vec, surfaceMesh->getCells())};
    const MatrixXc A = Matrix::getMatrix(physics.k, surfaceMesh->getCells());
    // preconditioning
    const VectorXc JacobiAux = A.diagonal().cwiseInverse();
    const VectorXc newRHS = JacobiAux.cwiseProduct(b);
    const MatrixXc newMatrix = JacobiAux.asDiagonal() * A;

    auto method = Eigen::GMRES<MatrixXc>{};
    method.setMaxIterations(20000);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(1e-5);
    method.set_restart(2000);
    method.compute(newMatrix);
    const auto j = VectorXc{method.solveWithGuess(newRHS, newRHS)};
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    surfaceMesh->fillJ(j);

    // создаем окружающую сетку

    Containers::vector<Mesh::Point> nodes;
    nodes.reserve(N_volume * N_volume);
    cartesian_productXY(std::ranges::views::iota(0, N_volume),
                        std::ranges::views::iota(0, N_volume),
                        std::back_inserter(nodes), N_volume, h_volume);

    const auto cellView = nodes | std::views::transform([](const Mesh::Point &p) { return Mesh::Node{p}; });

    Mesh::VolumeMesh volumeMesh{*surfaceMesh, {cellView.begin(), cellView.end()}};
    volumeMesh.setName("volume_mesh_" + std::to_string(N_volume));
    volumeMesh.calculateAll(physics.E0, physics.k_vec, physics.k);

    VTK::surface_snapshot(1, *surfaceMesh, Pathes::examples + "plane_basis_experiment/");

    VTK::volume_snapshot(1, volumeMesh, Pathes::examples + "plane_basis_experiment/");

    delete surfaceMesh;

    return 0;
}
