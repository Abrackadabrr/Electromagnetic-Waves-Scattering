//
// Created by evgen on 29.07.24.
//

#include "meshes/plate/PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "math/SurfaceField.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/MathConstants.hpp"
#include "mesh/Parser.hpp"
#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "operators/Operators.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <iostream>

using namespace EMW;
using namespace EMW::Types;

int main() {
    // сетка на пластинке
    int N1 = 82;
    scalar h1 = 1. / (N1 - 1);
    int N2 = 82;
    scalar h2 = 1. / (N2 - 1);
    auto surfaceMesh = Mesh::SurfaceMesh{EMW::Examples::Plate::generateRectangularMesh(N1, N2, h1, h2)};
    surfaceMesh.setName("surface_mesh_" + std::to_string(N1) + "_x_" + std::to_string(N2));

    // начальное поле на пластинке
    // задаем дальта-функцию в точке на поверхности: это касательный вектор, который по модулю равен 1/площадь ячейки
    // тогда интеграл от такой функции будет как раз равен единице
    Containers::vector<Types::Vector3c> delta_function(surfaceMesh.getCells().size(), Types::Vector3c::Zero());
    Types::index delta_support = surfaceMesh.getCells().size() / 2;
    delta_function[delta_support] =
            (1 / surfaceMesh.getCells()[delta_support].area_) *
            surfaceMesh.getCells()[delta_support].tau[0];
//    std::cout << delta_function[delta_support] << std::endl;
    Math::SurfaceField rhs(surfaceMesh, delta_function);
    rhs.setName("delta_rhs");

    const scalar k = 4 * Math::Constants::PI<scalar>();
    // собираем СЛАУ
    const VectorXc b = rhs.asSLAERHS();
    const MatrixXc A = Matrix::getMatrix(k, surfaceMesh);
    // предобуславливание
    const VectorXc JacobiAux = A.diagonal().cwiseInverse();
    const MatrixXc A_prec = JacobiAux.asDiagonal() * A;
    const VectorXc b_prec = JacobiAux.cwiseProduct(b);

    // решаем
    auto method = Eigen::GMRES<MatrixXc>{};
    method.setMaxIterations(20000);
    method.setTolerance(5e-4);
    method.set_restart(10000);
    method.compute(A_prec);
    const VectorXc j_vec = VectorXc{method.solveWithGuess(b_prec, b_prec)};
    auto j = Math::SurfaceField::TangentField(surfaceMesh, j_vec);
    j.setName("solution");
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;

    VTK::surface_snapshot(surfaceMesh, Pathes::studies + "plane/analytical_solution/");

    VTK::field_snapshot(j, Pathes::studies + "plane/analytical_solution/");

    // строим аналитическое решение
    const auto analytical = [rhs, k] (const Types::Vector3d &point) {
        return EMW::Operators::K1<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(point, k, rhs) +
        EMW::Operators::K0<DefiniteIntegrals::GaussLegendre::Quadrature<8>>(point, k, rhs);
    };

    auto analytical_solution = Math::SurfaceField(surfaceMesh, analytical);
    analytical_solution.setName("analytical");
    VTK::field_snapshot(analytical_solution, Pathes::studies + "plane/analytical_solution/");
}
