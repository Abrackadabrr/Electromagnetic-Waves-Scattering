//
// Created by evgen on 08.02.24.
//

#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "math/MathConstants.hpp"
#include "meshes/plate/PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <iostream>
#include <math/integration/gauss_quadrature/GaussLegenderPoints.hpp>
#include <operators/OperatorK.hpp>

using namespace EMW;
using namespace EMW::Types;

template <typename Range1, typename Range2, typename OutputIterator>
void cartesian_productXY(Range1 const &r1, Range2 const &r2, OutputIterator out, Types::index N, Types::scalar h) {
    using std::begin;
    using std::end;

    for (auto i = begin(r1); i != end(r1); ++i) {
        for (auto j = begin(r2); j != end(r2); ++j) {

            *out++ = Types::Vector3d{static_cast<Types::scalar>(*i) - static_cast<Types::scalar>(N / 2.),
                                     static_cast<Types::scalar>(*j) - static_cast<Types::scalar>(N / 2.), 0} *
                     h;
        }
    }
}

inline Types::Vector3c operatorK_in_point(const Math::SurfaceVectorField &field, const Types::complex_d k,
                                          const Mesh::point_t &point) {
    return EMW::OperatorK::K1<DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>>(point, k, field) +
           EMW::OperatorK::K0<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(point, k, field);
}

inline Types::Vector3c getE_in_point(const Math::SurfaceVectorField &j_e, const Types::complex_d k,
                                     const Mesh::point_t &point) {
    const Types::complex_d mul = Math::Constants::i / k;
    const auto value = operatorK_in_point(j_e, k, point);
    // std::cout << value << std::endl;
    return value;
}

int main() {

    // физика
    EMW::Physics::planeWaveCase physics(Vector3d{0, 1, 0}.normalized(),
                                        Types::complex_d{4 * Math::Constants::PI<scalar>(), 0},
                                        Vector3d{1, 0, 0}.normalized());
    const auto initial_field_function = [physics](const Mesh::point_t &point) -> Vector3c {
        return physics.value(point);
    };

    int N_volume = 41;
    scalar h_volume = 0.075;

    int N = 41;
    scalar h = 1. / (N - 1);
    // сетка
    auto surfaceMesh = Mesh::SurfaceMesh{EMW::Examples::Plate::generatePlatePrimaryMesh(N, h)};

    surfaceMesh.setName("mesh" + std::to_string(N));

    // след падающего поля на расчетной поверхности
    const Math::SurfaceVectorField incidentField(surfaceMesh, initial_field_function);

    // на мелкой
    const VectorXc b = -incidentField.asVector();
    const MatrixXc A = Matrix::getMatrixK(physics.k, surfaceMesh);
    //    std::cout << A3.operatorNorm() << std::endl << std::endl;
    //    std::cout << newMatrix << std::endl << std::endl;
    //    std::cout << newMatrix.inverse().norm() << std::endl;
    std::cout << "Matrix has been built" << std::endl;
    //    std::cout << "Matrix conditional value: " << (A3.operatorNorm()) * (A3.inverse().operatorNorm()) << std::endl;
    //    std::cout << "Matrix result value: " << newMatrix.operatorNorm() * (newMatrix.inverse().operatorNorm()) <<
    //    std::endl;

    auto method = Eigen::GMRES<MatrixXc>{};
    method.setMaxIterations(1000);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(1e-2);
    method.set_restart(1000);
    method.compute(A);
    const auto j = VectorXc{method.solve(b)};
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;

    auto j_field = Math::SurfaceVectorField::TangentField(surfaceMesh, j);
    j_field.setName("j");

    // создаем окружающую сетку

    Containers::vector<Mesh::point_t> nodes;
    nodes.reserve(N_volume * N_volume);
    cartesian_productXY(std::ranges::views::iota(0, N_volume), std::ranges::views::iota(0, N_volume),
                        std::back_inserter(nodes), N_volume, h_volume);

    const auto calculated_field_view =
        nodes | std::views::transform([&](const Types::Vector3d& p)->Types::Vector3c { return getE_in_point(j_field, physics.k, p) + initial_field_function(p); });
    const Containers::vector<Types::Vector3c> calculated_field{calculated_field_view.begin(),
                                                               calculated_field_view.end()};

    const std::string path_to_dir = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/"
                                    "vtk_files/examples/plane/verification_of_operator_realisation/";

    VTK::united_snapshot({j_field}, {}, surfaceMesh, path_to_dir);
    VTK::field_in_points_snapshot({calculated_field}, {"E"}, nodes, "verif_E_mesh", path_to_dir);

    return 0;
}
