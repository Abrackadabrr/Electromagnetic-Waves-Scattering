//
// Created by evgen on 27.06.24.
//

#include "meshes/plate/PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/MathConstants.hpp"
#include "mesh/Parser.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "examples/pathes.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <iostream>

using namespace EMW;
using namespace EMW::Types;

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

std::array<Types::Vector3d, 3> equalLocalBasis_1(const Mesh::IndexedCell &cell) {
    return {
        Vector3d{1, 0, 0},
        Vector3d{0, 1, 0},
        cell.normal
    };
}

std::array<Types::Vector3d, 3> equalLocalBasis_2(const Mesh::IndexedCell &cell) {
    return {
        Vector3d{1, 1, 0}.normalized(),
        Vector3d{-1, 1, 0}.normalized(),
        cell.normal
    };
}

Math::SurfaceField solve(const Math::SurfaceField &incidentField, const Physics::planeWaveCase& physics) {
    const Mesh::SurfaceMesh &surfaceMesh = incidentField.getManifold();
    const VectorXc b3 = incidentField.asSLAERHS();
    const MatrixXc A3 = Matrix::getMatrix(physics.k, surfaceMesh);

    auto method = Eigen::GMRES<MatrixXc>{};
    method.setMaxIterations(30000);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(1e-5);
    method.set_restart(3000);
    method.compute(A3);
    const VectorXc j_vec = VectorXc{method.solve(b3)};
    auto j = Math::SurfaceField::TangentField(surfaceMesh, j_vec);
    j.setName("j");
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    return j;
}

int main() {
    // сетка на пластинке
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/plate/triangulated/1191_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/plate/triangulated/2256_cells.csv";
    const EMW::Types::index nNodes = 1191;
    const EMW::Types::index nCells = 2256;
    auto surfaceMesh = Mesh::SurfaceMesh{EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells)};

    surfaceMesh.setName("surface_mesh_triangular_basis_1" + std::to_string(nNodes));
    surfaceMesh.customLocalBasis(equalLocalBasis_1);

    // физика
    EMW::Physics::planeWaveCase physics{
            Vector3d{0, 1, 0}.normalized(),
            4 * Math::Constants::PI<scalar>(),
            Vector3d{0, 0, 1}.normalized()
    };
    const auto initial_field_function = [physics](const Mesh::point_t & point) -> Vector3c {
        return physics.value(point);
    };

    // след падающего поля на расчетной поверхности
    const Math::SurfaceField incidentField(surfaceMesh, initial_field_function);
    const auto j = solve(incidentField, physics);

    VTK::united_snapshot(surfaceMesh, {j}, Pathes::examples + "plane/new_discretization/");

    // изменение параметров сетки и повторное решение методом
    surfaceMesh.customLocalBasis(equalLocalBasis_2);
    surfaceMesh.setName("surface_mesh_triangular_basis_2" + std::to_string(nNodes));
    const auto j2 = solve(incidentField, physics);

    auto diff = j - j2;
    diff.setName("diff");

    VTK::united_snapshot(surfaceMesh, {j2, diff}, Pathes::examples + "plane/new_discretization/");
}
