//
// Created by evgen on 08.02.24.
//

#include "PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/MathConstants.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

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
    int N_volume = 81;
    scalar h_volume = 0.075;

    int N = 31;
    scalar h = 1. / (N-1);
    // сетка
    auto *surfaceMesh = new Mesh::SurfaceMesh{EMW::Examples::Plate::generatePlatePrimaryMesh(N, h)};

    surfaceMesh->setName("surface_mesh_" + std::to_string(N));

    // физика
    physicalConditions physics{
            .E0 = Vector3d{0, 0, 1}.normalized(),
            .k = complex_d{4 * Math::Constants::PI<scalar>(), 0},
            .k_vec = Vector3d{-1, -1, 0}.normalized()
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

    // создаем окружающую сетку

    Containers::vector<Mesh::Point> nodes;
    nodes.reserve(N_volume * N_volume * N_volume);
    cartesian_productXY(std::ranges::views::iota(0, N_volume),
                        std::ranges::views::iota(0, N_volume),
                        std::back_inserter(nodes), N_volume, h_volume);

    const auto cellView = nodes | std::views::transform([](const Mesh::Point &p) { return Mesh::Node{p}; });

    Mesh::VolumeMesh volumeMesh{*surfaceMesh, {cellView.begin(), cellView.end()}};
    volumeMesh.setName("volume_mesh_" + std::to_string(N_volume));
    volumeMesh.calculateAll(physics.E0, physics.k_vec, physics.k);

    VTK::test_snapshot(1, *surfaceMesh,
                       "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plane/volume/");

    VTK::volume_snapshot(1, volumeMesh,
                         "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plane/volume/");
}
