//
// Created by evgen on 25.10.2024.
//

#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "math/MathConstants.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "mesh/Algorithms.hpp"
#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/Utils.hpp"
#include "meshes/plate/PlateGrid.hpp"
#include "operators/Operators.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"

#include "iterative_solvers/SimpleIterationMethod.hpp"
#include "slae_generation/Preconditioning.hpp"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SVD>
#include <unsupported/Eigen/IterativeSolvers>

#include <chrono>
#include <filesystem>
#include <iostream>

using namespace EMW;
using namespace EMW::Types;

Math::SurfaceVectorField operatorK(const Math::SurfaceVectorField &field, complex_d k) {
    const auto analytical = [field, k](const Types::Vector3d &point) -> Vector3c {
        return EMW::Operators::K1_singularityExtraction<DefiniteIntegrals::GaussLegendre::Quadrature<2, 2>>(point, k,
                                                                                                            field) +
               EMW::Operators::K0<DefiniteIntegrals::GaussLegendre::Quadrature<2>>(point, k, field);
    };

    return {field.getManifold(), analytical};
}

Types::scalar multiplier(const Types::Vector3d &x) {
    if (x.norm() > 0.02)
        return 1;
    return 0;
}

Types::Vector3c getConstantRHS(const Types::Vector3d &x) {
    if (x.norm() < 0.1)
        return {complex_d{1, 0}, {0, 0}, {0, 0}};
    return Vector3c::Zero();
}

std::pair<VectorXc, int> solve(const MatrixXc &A, const VectorXc &b, scalar tolerance) {
    auto method = Eigen::GMRES<MatrixXc, Eigen::IdentityPreconditioner>{};

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

std::pair<VectorXc, int> solveSimpleIterations(const MatrixXc &A, const MatrixXc &P, const VectorXc &b,
                                               scalar tolerance) {
    auto method = EMW::IterativeSolvers::SimpleIterationMethod<MatrixXc, VectorXc>{};
    method.tolerance = tolerance;

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    VectorXc j_vec = VectorXc{method.solve(P, A, b, VectorXc::Zero(A.rows()))};

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::ratio<1, 1>> elapsed_seconds = end - start;
    std::cout << "Время решения простой итерации: " << elapsed_seconds.count() << std::endl;
    std::cout << "total iterations: " << method.iterations << std::endl;
    std::cout << "tolerance: " << method.tolerance << std::endl;
    std::cout << "total error: " << method.error << std::endl;
    std::cout << "Info: " << method.info() << std::endl;
    return {j_vec, method.iterations};
}

int main() {
    Eigen::setNbThreads(14);

    int N1 = 81;
    scalar h1 = 1. / (N1 - 1);
    int N2 = 81;
    scalar h2 = 1. / (N2 - 1);
    auto surfaceMesh = Examples::Plate::generateRectangularMesh(N1, N2, h1, h2);


    surfaceMesh.setName("rect" + std::to_string(N1) + "_x_" + std::to_string(N2));

    // падающее поле на треугольной сетке
    constexpr complex_d k = 4 * Math::Constants::PI<scalar>() * complex_d{1, 0.05};
    const Vector3d e = Vector3d{1, 0, 0}.normalized();
    EMW::Physics::planeWaveCase physics{e, k, Vector3d{0, 0, -1}.normalized()};
    const auto initial_field_function = [physics](const Mesh::point_t &point) -> Vector3c {
        // std::cout << physics.value(point) << std::endl;
        return physics.value(point);
    };

    auto rhsField = Math::SurfaceVectorField(surfaceMesh, getConstantRHS);

    // Ставится задача K[u] = m x n, которая на бесконечной плоскости
    // имеет решение вида u = (4/k^2) n x K[m].
    const VectorXc b = rhsField.asSLAERHS();
    const MatrixXc A = Matrix::getMatrix(k, surfaceMesh);
    constexpr scalar tolerance = 1e-2;
    std::cout << "Solving without preconditioner" << std::endl;
    const auto sol = solve(A, b, tolerance);

    // файл для записи
    // std::ofstream fileout;
    const std::string path = Pathes::studies + "plane/preconditioning/rectangular/" + surfaceMesh.getName();
    // std::filesystem::create_directories(path);
    // fileout.open(path + "data.csv");
    // fileout << "r,i" << std::endl;

    for (scalar r = 1.409; r < 1.41; r += 0.1) {
        std::cout << "1x1, radius = " << r << std::endl;

        const MatrixXc P = (1./2) * (Matrix::Preconditioning::getPreconditiotner(surfaceMesh, r, k) + A);
        const MatrixXc AP =  A * P; // potentially too long operation

        // 2) Решим систему с правой частью и посмотрим на количество итераций
        const auto sol_p = solve(AP, b, tolerance);
        const VectorXc solution = P * sol_p.first;
        // 3) Посмотрим на относительную ошибку в решениях
        const VectorXc diff = sol.first - solution;
        const scalar rel_err = diff.norm() / solution.norm();
        std::cout << "Rel error: " << rel_err << '\n' << std::endl;

        // fileout << r << "," << sol_p.second << std::endl;

        // 4) Посмотрим на поля на повержности через Paraview
        Math::SurfaceVectorField solField = Math::SurfaceVectorField::TangentField(surfaceMesh, sol.first);
        solField.setName("solution_straightforward");
        Math::SurfaceVectorField sol_pField = Math::SurfaceVectorField::TangentField(surfaceMesh, solution);
        sol_pField.setName("solution_preconditioned");
        VTK::united_snapshot({solField, sol_pField, rhsField}, {}, surfaceMesh, path + "/");
    }
    // fileout.close();
}
