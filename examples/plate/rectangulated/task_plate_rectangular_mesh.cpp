//
// Created by evgen on 08.02.24.
//

#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "math/MathConstants.hpp"
#include "math/fields/SurfaceVectorField.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/VolumeMesh.hpp"
#include "meshes/plate/PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <iostream>
#include <ranges>

using namespace EMW;
using namespace EMW::Types;

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

std::array<Types::Vector3d, 3> equalLocalBasis_1(const Mesh::IndexedCell &cell) {
    return {Vector3d{1, 0, 0}, Vector3d{0, 1, 0}, cell.normal};
}

std::array<Types::Vector3d, 3> equalLocalBasis_2(const Mesh::IndexedCell &cell) {
    return {Vector3d{1, 1, 0}.normalized(), Vector3d{-1, 1, 0}.normalized(), cell.normal};
}

Math::SurfaceVectorField solve(const Math::SurfaceVectorField &incidentField, const Physics::planeWaveCase &physics) {
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
    auto j = Math::SurfaceVectorField::TangentField(surfaceMesh, j_vec);
    j.setName("j");
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    return j;
}

int main() {
    int N1 = 41;
    scalar h1 = 1. / (N1 - 1);
    int N2 = 31;
    scalar h2 = 1. / (N2 - 1);

    // физика
    EMW::Physics::planeWaveCase physics{Vector3d{0, 1, 0}.normalized(), 4 * Math::Constants::PI<scalar>(),
                                        Vector3d{0, 0, 1}.normalized()};

    const auto initial_field = [physics](const Mesh::point_t & point) {
        return physics.value(point);
    };

    auto surfaceMesh = Examples::Plate::generateRectangularMesh(N1, N2, h1, h2);

    surfaceMesh.setName("surface_mesh_rect_basis_1");
    surfaceMesh.customLocalBasis(equalLocalBasis_1);

    // след падающего поля на расчетной поверхности
    const Math::SurfaceVectorField incidentField(surfaceMesh, initial_field);
    const Math::SurfaceVectorField j = solve(incidentField, physics);

    VTK::united_snapshot({j}, {}, surfaceMesh, Pathes::examples + "plane/new_discretization/");

    /*// изменение параметров сетки и повторное решение методом
    surfaceMesh.customLocalBasis(equalLocalBasis_2);
    surfaceMesh.setName("surface_mesh_rect_basis_2");
    const auto j2 = solve(incidentField, physics);

    auto diff = j - j2;
    auto diff_norm = diff.fieldNorm();
    diff.setName("diff");
    diff_norm.setName("diff_norm");

    VTK::united_snapshot({j2, diff}, {diff_norm}, surfaceMesh, Pathes::examples + "plane/new_discretization/");
    */

}
