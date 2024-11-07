//
// Created by evgen on 25.10.2024.
//

#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "math/MathConstants.hpp"
#include "math/fields/Utils.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "mesh/Algorithms.hpp"
#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/Utils.hpp"
#include "meshes/plate/PlateGrid.hpp"
#include "operators/Operators.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"

#include "Functions.hpp"
#include "slae_generation/Preconditioning.hpp"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SVD>
#include <unsupported/Eigen/IterativeSolvers>

#include <iostream>

using namespace EMW;
using namespace EMW::Types;

Math::SurfaceVectorField operatorK(const Math::SurfaceVectorField &field, scalar k) {
    const auto analytical = [field, k](const Types::Vector3d &point) -> Vector3c {
        return EMW::Operators::K1_singularityExtraction<DefiniteIntegrals::GaussLegendre::Quadrature<2, 2>>(point, k,
                                                                                                            field) +
               EMW::Operators::K0<DefiniteIntegrals::GaussLegendre::Quadrature<2>>(point, k, field);
    };

    return {field.getManifold(), analytical};
}

/**
 * Локальные базисы на пластине
 * @param cell -- объект ячейки
 * @return массив векторов, составляющих правую тройку
 */
std::array<Types::Vector3d, 3> basis(const Mesh::IndexedCell &cell) {
    return {Vector3d{1, 0, 0}, Vector3d{0, 1, 0}, cell.normal};
}

/**
 * Точки коллокации для треугольных элементов на ячейке cell
 * @param cell -- объект ячейки
 * @return точку коллокации
 */
Types::Vector3d collocationPoint(const Mesh::IndexedCell &cell) {
    const auto vertex = cell.getVertexAsArray();
    const Vector3d center = (1. / 3.) * (vertex[0] + vertex[1] + vertex[2]);
    return center;
}

Types::scalar multiplier(const Types::Vector3d &x) {
    if (x.norm() > 0.02)
        return 1;
    return 0;
}

/**
 * В этой функции я хочу посмотреть на разность ((предобуславливатель * матрица) - единичная)
 * @return форбениусова норма ((предобуславливатель * матрица) - единичная)
 */
Types::scalar getDifferenceForbenius(const Mesh::SurfaceMesh &surfaceMesh, const scalar k) {
    std::cout << "Сборка А" << std::endl;
    const MatrixXc A = Matrix::getMatrix(k, surfaceMesh);
    std::cout << "Сборка В" << std::endl;
    const MatrixXc B = Matrix::Preconditioning::getPreconditiotner(surfaceMesh, 0, k);
    std::cout << "Перемножение" << std::endl;
    Eigen::setNbThreads(14);
    const MatrixXc diff = B * A;
    std::cout << "Считаем норму" << std::endl;
    std::cout << "Норма А: " << A.norm() << std::endl;
    return (diff - MatrixXc::Identity(B.rows(), B.cols())).eval().norm();
}

std::pair<Types::scalar, Types::scalar> getDifferenceSpectral(const Mesh::SurfaceMesh &surfaceMesh, const scalar k) {
    std::cout << "Сборка А" << std::endl;
    const MatrixXc A = Matrix::getMatrix(k, surfaceMesh);
    std::cout << "Сборка В" << std::endl;
    const MatrixXc B = Matrix::Preconditioning::getPreconditiotner(surfaceMesh, 0, k);
    std::cout << "Перемножение" << std::endl;
    Eigen::setNbThreads(14);
    const MatrixXc diff = B * A;
    std::cout << "Считаем норму" << std::endl;
    Eigen::JacobiSVD<MatrixXc> svd(A);
    std::cout << svd.singularValues() << std::endl;
    double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
    std::cout << cond << std::endl;

    // Eigen::JacobiSVD<MatrixXc> svd_d{diff};
    // std::cout << "svd_d computed: " << std::endl;
    // Eigen::JacobiSVD<MatrixXc> svd_a{A};
    // std::cout << svd_d.singularValues() << std::endl;
    // std::cout << "A: Max singular_value: " << svd_a.singularValues().cwiseAbs().maxCoeff()
    // << " ; Min singular_value: " << svd_a.singularValues().cwiseAbs().minCoeff() << std::endl;

    return {0, 0};
}

std::pair<VectorXc, int> solve(const MatrixXc &A, const VectorXc &b, scalar tolerance) {
    auto method = Eigen::GMRES<MatrixXc>{};

    method.setMaxIterations(30000);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(tolerance);
    method.set_restart(30000);

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    method.compute(A);
    VectorXc j_vec = VectorXc{method.solve(b)};

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::ratio<1, 1>> elapsed_seconds = end - start;
    std::cout << "Время решения GMRES: " << elapsed_seconds.count() << std::endl;
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "tolerance: " << method.tolerance() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    return {j_vec, method.iterations()};
}

/*
 * Сеточный аналог векторной дельта-функции, заданной на ячейках сетки
 */
Vector3c discrete_delta(const Mesh::IndexedCell &cell) {
    const Mesh::point_t support{0.1, 0, 0};
    return Vector3d{0, 1, 0} * Mesh::Algorithm::PointInTriangle(support, cell.getVertex()) / cell.area_;
}

int main() {
    Eigen::setNbThreads(14);

    // треугольная сетка не пластинке
    auto surfaceMesh = Mesh::Utils::loadTriangularMesh(1440, 2742, "1_1");
    surfaceMesh.customLocalBasis(basis);
    surfaceMesh.customCollocationPpoints(collocationPoint);

    // int N1 = 30;
    // scalar h1 = 1. / (N1 - 1);
    // int N2 = 30;
    // scalar h2 = 1. / (N2 - 1);
    // auto surfaceMesh = Examples::Plate::generateRectangularMesh(N1, N2, h1, h2);
    //
    // surfaceMesh.setName("rect" + std::to_string(N1) + "_x_" + std::to_string(N2));
    // surfaceMesh.customLocalBasis(equalLocalBasis_1);

    // падающее поле на треугольной сетке
    constexpr scalar k = 4 * Math::Constants::PI<scalar>();
    const Vector3d e = Vector3d{0, 1, 0}.normalized();

    // EMW::Physics::planeWaveCase physics{e, k, Vector3d{0, 0, -1}.normalized()};
    // const auto initial_field_function = [physics](const Mesh::point_t &point) -> Vector3c {
    //     return physics.value(point);
    // };

    auto rhsField =
        Math::SurfaceVectorField(surfaceMesh, discrete_delta);

    // Ставится задача K[u] = m x n, которая на бесконечной плоскости
    // имеет решение вида u = (4/k^2) n x K[m].
    const VectorXc b = rhsField.asSLAERHS();
    const MatrixXc A = Matrix::getMatrix(k, surfaceMesh);
    const scalar tolerance = 1e-2;
    std::cout << "Solving without preconditioner" << std::endl;
    const auto sol = solve(A, b, tolerance);

    // файл для записи
    std::ofstream fileout;
    fileout.open (Pathes::studies + "plane/preconditioning/triangular/1_x_1_2742_E_is_delta_a_little_shifted_log.csv");
    fileout << "r,i" << std::endl;

    for (scalar r = 0.3; r < 1.41; r += 0.1) {
        std::cout << "1x1, radius = " << r << std::endl;

        const MatrixXc P = Matrix::Preconditioning::getPreconditiotner(surfaceMesh, r, k);
        const MatrixXc AP = A*P; // potentially too long operation
        // такое умножение на константу нужно для того, чтобы не ум

        // 2) Решим систему с правой частью и посмотрим на количество итераций
        const auto sol_p = solve(AP, b, tolerance);
        const VectorXc solution = P * sol_p.first;
        // 3) Посмотрим на относительную ошибку в решениях
        const VectorXc diff = sol.first - solution;
        const scalar rel_err = diff.norm() / solution.norm();
        std::cout << "Rel error: " << rel_err << '\n' << std::endl;

        fileout << r << "," << sol_p.second << std::endl;

        // 4) Посмотрим на поля на повержности через Paraview
        // Math::SurfaceVectorField solField = Math::SurfaceVectorField::TangentField(surfaceMesh, sol);
        // solField.setName("solution_straightforward");
        // Math::SurfaceVectorField sol_pField = Math::SurfaceVectorField::TangentField(surfaceMesh, sol_p);
        // sol_pField.setName("solution_preconditioned");
        // VTK::united_snapshot({solField, sol_pField}, {}, surfaceMesh,
        //     Pathes::studies + "plane/preconditioning/triangular/");
    }
    fileout.close();
}
