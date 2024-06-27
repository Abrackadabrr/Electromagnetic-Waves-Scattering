//
// Created by evgen on 27.06.24.
//

#include "examples/plate/PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/MathConstants.hpp"
#include "mesh/Parser.hpp"

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

template<typename Range1, typename Range2, typename OutputIterator>
void cartesian_productXY(Range1 const &r1, Range2 const &r2, OutputIterator out, Types::index N,
                         Types::scalar h) {
    using std::begin;
    using std::end;

    for (auto i = begin(r1); i != end(r1); ++i) {
        for (auto j = begin(r2); j != end(r2); ++j) {

            *out++ = Types::Vector3d{0,
                                     static_cast<Types::scalar>(*i) - static_cast<Types::scalar>(N / 2.),
                                     static_cast<Types::scalar>(*j) - static_cast<Types::scalar>(N / 2.)} * h;
        }
    }
}


int main() {
    int N_volume = 81;
    scalar h_volume = 0.075;

    // сетка на пластинке
    const std::string nodesFile = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/plate/triangulated/144_nodes.csv";
    const std::string cellsFile = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/plate/triangulated/246_cells.csv";
    const EMW::Types::index nNodes = 144;
    const EMW::Types::index nCells = 246;
    auto *surfaceMesh = new Mesh::SurfaceMesh{EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells)};

    surfaceMesh->setName("surface_mesh_triangular" + std::to_string(nNodes));

    // физика
    physicalConditions physics{
            .E0 = Vector3d{0, 1, 0}.normalized(),
            .k = complex_d{4 * Math::Constants::PI<scalar>(), 0},
            .k_vec = Vector3d{0, 0, 1}.normalized()
    };

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
    volumeMesh.setName("volume_mesh_triangular" + std::to_string(N_volume));
    volumeMesh.calculateAll(physics.E0, physics.k_vec, physics.k);

    VTK::test_snapshot(1, *surfaceMesh,
                       "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plane/volume/");

    VTK::volume_snapshot(1, volumeMesh,
                         "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plane/volume/");
}
