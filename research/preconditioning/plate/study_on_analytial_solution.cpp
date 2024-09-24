//
// Created by evgen on 29.07.24.
//

#include "meshes/plate/PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "math/fields/SurfaceField.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/MathConstants.hpp"
#include "mesh/Parser.hpp"
#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "operators/Operators.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "mesh/Algorithms.hpp"

#include "Functions.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <iostream>

using namespace EMW;
using namespace EMW::Types;

Vector3c discrete_delta(const Vector3d &e, const Mesh::IndexedCell &cell) {
    const Mesh::point_t support{0, 0, 0};
    return e * Mesh::Algorithm::PointInTriangle(support, cell.getVertex()) / cell.area_;
}

Math::SurfaceField solve(const Math::SurfaceField &rhs, scalar k) {
    const auto &surfaceMesh = rhs.getManifold();
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

    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;

    auto j = Math::SurfaceField::TangentField(surfaceMesh, j_vec);
    return j;
}

Math::SurfaceField operatorK(const Math::SurfaceField &field, scalar k) {
    const auto analytical = [field, k](const Types::Vector3d &point) {
        return EMW::Operators::K1_singularityExtraction<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8 >>(point, k,
                                                                                                             field) +
               EMW::Operators::K0<
                       DefiniteIntegrals::GaussLegendre::Quadrature<
                               8 >>(point, k, field);
    };

    return Math::SurfaceField(field.getManifold(), analytical);
}

int main() {
    // треугольная сетка не пластинке
    const std::string nodesFile = "/media/evgen/SecondLinuxDisk/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/plate/triangulated/1191_nodes.csv";
    const std::string cellsFile = "/media/evgen/SecondLinuxDisk/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/plate/triangulated/2256_cells.csv";
    const EMW::Types::index nNodes = 1191;
    const EMW::Types::index nCells = 2256;
    auto triangularMesh = Mesh::SurfaceMesh{EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells)};
    triangularMesh.setName("triangular_mesh_" + std::to_string(nCells));

    const Vector3d e{0, 1, 0};
    // поле на треугольной сетке
    auto rhs_triangular = Math::SurfaceField(triangularMesh,
                                             [e](const Mesh::IndexedCell &cell) { return discrete_delta(e, cell); });
    rhs_triangular.setName("delta_rhs");

    // расчет
    const scalar k = 4 * Math::Constants::PI<scalar>();
//    auto j_tri = solve(rhs_triangular, k);
//    j_tri.setName("solution");

    // строим аналитическое решение
    auto analytical_solution_tri = operatorK(rhs_triangular, k);
    analytical_solution_tri.setName("an_num");
    auto analytical_tri = Math::SurfaceField(triangularMesh,
                                             [e, k](const Mesh::IndexedCell &cell) {
                                                 return analyticalInfinitePlate(e, cell, k);
                                             }) + rhs_triangular;
    analytical_tri.setName("an_math");

    auto analytical_difference = analytical_solution_tri - analytical_tri;
    analytical_difference.setName("difference_an");

    // сохраняем результаты
    VTK::united_snapshot(triangularMesh,
                         {rhs_triangular, analytical_solution_tri, analytical_tri, analytical_difference},
                         Pathes::studies + "plane/analytical_solution/");
}
