//
// Created by evgen on 08.02.24.
//

#include "PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "Eigen/IterativeLinearSolvers"
#include "mesh/VolumeMesh.hpp"

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

template<typename Range1, typename Range2, typename Range3, typename OutputIterator>
void cartesian_product3D(Range1 const &r1, Range2 const &r2, Range3 const &r3, OutputIterator out, Types::index N,
                         Types::scalar h) {
    using std::begin;
    using std::end;

    for (auto i = begin(r1); i != end(r1); ++i) {
        for (auto j = begin(r2); j != end(r2); ++j) {
            for (auto k = begin(r3); k != end(r3); ++k) {
                *out++ = Types::Vector3d{static_cast<Types::scalar>(*i) - static_cast<Types::scalar>(N / 2.),
                                         static_cast<Types::scalar>(*j) - static_cast<Types::scalar>(N / 2.),
                                         static_cast<Types::scalar>(*k) - static_cast<Types::scalar>(N / 2.)} * h;
            }
        }
    }
}

int main() {
    int N_volume = 21;
    scalar h_volume = 0.07;

    int N = 51;
    scalar h = 0.02;
    // сетка
    auto *surfaceMesh = new Mesh::SurfaceMesh{EMW::Examples::Plate::generatePlatePrimaryMesh(N, h)};

    surfaceMesh->setName("mesh_" + std::to_string(N) + "_difraction");

    // поляризация
    physicalConditions physics{
            .E0 = Vector3d{0, 0, 1},
            .k = complex_d{1, 0},
            .k_vec = Vector3d{-1, 0, 0}
    };

    const MatrixXd *A3 = new MatrixXd{Matrix::getMatrix(physics.k, surfaceMesh->getCells()).real()};
    const VectorXc *b3 = new VectorXc{Matrix::getRHS(physics.E0, physics.k_vec, surfaceMesh->getCells())};
    const auto cg3 = new Eigen::ConjugateGradient<MatrixXd, Eigen::Lower | Eigen::Upper>{*A3};
    const auto *j = new VectorXd{cg3->solve(b3->real())};

    surfaceMesh->fillJ(*j);

    // создаем окружающую сетку

    Containers::vector<Mesh::Point> nodes;
    nodes.reserve(N_volume * N_volume * N_volume);
    cartesian_product3D(std::ranges::views::iota(0, N_volume), std::ranges::views::iota(0, N_volume), std::ranges::views::iota(0, N_volume),
                        std::back_inserter(nodes), N_volume, h_volume);

    const auto cellView = nodes | std::views::transform([](const Mesh::Point &p) { return Mesh::Node{p}; });

    Mesh::VolumeMesh volumeMesh{*surfaceMesh, {cellView.begin(), cellView.end()}};
    volumeMesh.calculateAll(physics.E0, physics.k_vec, physics.k);

    VTK::volume_snapshot(1, volumeMesh,
                         "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plate/volume/");
}
