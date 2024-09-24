//
// Created by evgen on 08.02.24.
//

#include "slae_generation/MatrixGeneration.hpp"
#include "meshes/plate/PlateGrid.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "math/fields/SurfaceField.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/MathConstants.hpp"
#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"

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

int main() {
    int N_volume = 41;
    scalar h_volume = 0.075;

    int N1 = 41;
    scalar h1 = 1. / (N1 - 1);
    int N2 = 31;
    scalar h2 = 1. / (N2 - 1);

    // физика
    EMW::Physics::planeWaveCase physics(Vector3d{0, 1, 0}.normalized(),
                                        4 * Math::Constants::PI<scalar>(),
                                        Vector3d{1, 0, 0}.normalized());
    const auto initial_field = [physics](const Mesh::point_t & point) {
        return physics.value(point);
    };

    // сетка
    auto surfaceMesh = Mesh::SurfaceMesh{EMW::Examples::Plate::generateRectangularMesh(N1, N2, h1, h2)};
    surfaceMesh.setName("surface_mesh_" + std::to_string(N1) + "_x_" + std::to_string(N2));

    // след падающего поля на расчетной поверхности
    const Math::SurfaceField incidentField(surfaceMesh, initial_field);

    const VectorXc b3 = incidentField.asSLAERHS();
    const MatrixXc A3 = Matrix::getMatrix(physics.k, surfaceMesh);

    auto method = Eigen::GMRES<MatrixXc>{};
    method.setMaxIterations(20000);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(1e-5);
    method.set_restart(1000);
    method.compute(A3);
    const VectorXc j_vec = VectorXc{method.solve(b3)};
    auto j = Math::SurfaceField::TangentField(surfaceMesh, j_vec);
    j.setName("j");
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    surfaceMesh.fillJ(j_vec);

    VTK::united_snapshot(surfaceMesh, {j}, Pathes::examples + "plane/new_discretization/");
}
