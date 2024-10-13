//
// Created by evgen on 29.07.24.
//

#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "math/MathConstants.hpp"
#include "math/fields/Utils.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "mesh/Algorithms.hpp"
#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "meshes/plate/PlateGrid.hpp"
#include "operators/Operators.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"

#include "Functions.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <unsupported/Eigen/IterativeSolvers>

using namespace EMW;
using namespace EMW::Types;

/*
 * Сеточный аналог векторной дельта-функции, заданной на ячейках сетки
 */
Vector3c discrete_delta(const Vector3d &e, const Mesh::IndexedCell &cell) {
    const Mesh::point_t support{0, 0, 0};
    return e * Mesh::Algorithm::PointInTriangle(support, cell.getVertex()) / cell.area_;
}

Math::SurfaceVectorField operatorK(const Math::SurfaceVectorField &field, scalar k) {
    const auto analytical = [field, k](const Types::Vector3d &point) -> Vector3c {
        return EMW::Operators::K1_singularityExtraction<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(point, k,
                                                                                                            field) +
               EMW::Operators::K0<DefiniteIntegrals::GaussLegendre::Quadrature<8>>(point, k, field);
    };

    return Math::SurfaceVectorField(field.getManifold(), analytical);
}

std::array<Types::Vector3d, 3> equalLocalBasis_1(const Mesh::IndexedCell &cell) {
    return {Vector3d{1, 0, 0}, Vector3d{0, 1, 0}, cell.normal};
}

Types::scalar multiplier(const Types::Vector3d &x) {
    if (x.norm() > 0.02)
        return 1;
    return 0;
}
/*
 * Эта программа сравнивает аналитическую формулу, которая получается, если в правой части стоит дельта-функция,
 * с непосредственным расчетом оператора K от дельта-функции.
 * 1) Для сетки из кватратов:
 * Решения получаются схожими в мнимой части значений,
 * при этом в некоторой окрестности носителя дельты имеется относительное расхождение порядка 0.2
 * Связано ли это в дискретизацией оператора пока что неясно
 *
 * 2) Для сетки из треугольников:
 *
 */
int main() {
    // треугольная сетка не пластинке
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "plate/triangulated/3021_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "plate/triangulated/5840_cells.csv";
    const EMW::Types::index nNodes = 3021;
    const EMW::Types::index nCells = 5840;

    auto surfaceMesh = Mesh::SurfaceMesh{EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells)};
    surfaceMesh.setName("triangular_mesh_" + std::to_string(nCells));

    const Vector3d e{0, 1, 0};
    // падающее поле на треугольной сетке
    auto discrete_delta_field =
        Math::SurfaceVectorField(surfaceMesh, [e](const Mesh::IndexedCell &cell) { return discrete_delta(e, cell); });
    discrete_delta_field.setName("delta_discrete_delta");

    // расчет
    const scalar k = 4 * Math::Constants::PI<scalar>();

    // дискретизованный оператор K, примененный полю на сетке
    auto analytical_solution_numerical = operatorK(discrete_delta_field, k);
    analytical_solution_numerical.setName("K_numerical");
    analytical_solution_numerical.multiply(multiplier);

    // аналитическое применение оператора К к дельта функции на сетке
    auto analytical_solution_analytical =
        Math::SurfaceVectorField(surfaceMesh,
                                 [e, k](const Mesh::IndexedCell &cell) { return kernelFromPaper(e, cell, k); }) +
        discrete_delta_field;
    analytical_solution_analytical.setName("K_analytical");
    analytical_solution_analytical.multiply(multiplier);

    // расчитываем разницу между полями
    auto analytical_difference = analytical_solution_numerical - analytical_solution_analytical;
    analytical_difference.setName("difference");
    // и относительную ошибку по каждой из компонент
    auto relativeError = Math::FieldUtils::relativeError(analytical_solution_numerical, analytical_solution_analytical);
    relativeError.setName("relativeError");

    // сохраняем результаты
    VTK::united_snapshot({discrete_delta_field, analytical_solution_numerical,
        analytical_solution_analytical, analytical_difference},
        {relativeError}, surfaceMesh,
        Pathes::studies + "plane/analytical_solution/triangular/");
}
