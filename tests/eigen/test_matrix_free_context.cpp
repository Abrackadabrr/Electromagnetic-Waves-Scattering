//
// Created by evgen on 27.02.2025.
//

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <unsupported/Eigen/IterativeSolvers>

#include "math/matrix/eigen_matrix_replacement/MatrixReplacement.hpp"
#include "math/matrix/eigen_matrix_replacement/MatrixTraits.hpp"

#include <chrono>

#include <gtest/gtest.h>

using namespace EMW;
using namespace EMW::Types;

using Eigen::SparseMatrix;

using MatrixReplacement = EMW::Wrappers::MatrixReplacement<MatrixXc>;

Types::scalar simple_toeplitz_matrix(Types::index i, Types::index j) {
    return static_cast<integer>(i) - static_cast<integer>(j);
};

Types::MatrixX<scalar> get_block(Types::index i, Types::index j) {
    return Types::MatrixX<scalar>{{{1, 2}, {2, 3}}} * simple_toeplitz_matrix(i, j);
}

Types::MatrixX<complex_d> get_toeplitz_matrix(Types::index rows, Types::index cols) {
    Types::MatrixX<complex_d> result = Types::MatrixX<scalar>::Zero(2 * rows, 2 * cols);
#pragma omp parallel for collapse(2) num_threads(14)
    for (Types::index i = 0; i < rows; i++)
        for (Types::index j = 0; j < cols; j++)
            result.block(2 * i, 2 * j, 2, 2) = complex_d{1, 1} * get_block(i, j);
    return result;
}

TEST(MATRIX_FREE_COMTEXT, TEST_TO_CHECK_IT_WORKS) {
    int n = 1000;
    MatrixXc S = get_toeplitz_matrix(n, n) + 0.002 * n * MatrixXc::Identity(2 * n, 2 * n);
    MatrixXc Jacobi = S.diagonal().cwiseInverse().asDiagonal();
    S = S * Jacobi;
    std::cout << "Считаем число обусловленности" << std::endl;
    std::cout << S.norm() * S.inverse().norm() << std::endl;

    MatrixReplacement A;
    A.attachMyMatrix(S);

    EMW::Types::VectorXc b(2 * n), x;
    b.setRandom();

    // Solve Ax = b using various iterative solver with matrix-free version:
    {
        Eigen::BiCGSTAB<MatrixReplacement, Eigen::IdentityPreconditioner> bicg;
        bicg.compute(A);
        x = bicg.solve(b);
        std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error()
                  << std::endl;
    }

    {
        Eigen::GMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;
        gmres.compute(A);
        x = gmres.solve(b);
        std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error()
                  << std::endl;
    }
}


TEST(MATRIX_FREE_COMTEXT, BENCHMARK_1) {
    int n = 8000;
    MatrixXc S = get_toeplitz_matrix(n, n) + 0.002 * n * MatrixXc::Identity(2 * n, 2 * n);
    // MatrixXc Jacobi = S.diagonal().cwiseInverse().asDiagonal();
    // S = S * Jacobi;
    // std::cout << "Считаем число обусловленности" << std::endl;
    // std::cout << S.norm() * S.inverse().norm() << std::endl;

    MatrixReplacement A;
    A.attachMyMatrix(S);

    EMW::Types::VectorXc b(2 * n), x;
    b.setRandom();
    // Test case 1
    {
        Eigen::GMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;

        auto start = std::chrono::system_clock::now();

        gmres.compute(A);
        x = gmres.solve(b);

        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "GMRES with wrapped matrix: " << elapsed.count() << '\n';

        std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error()
                  << std::endl;
        std::cout << (A.my_matrix() * x - b).norm() << std::endl;
    }
    // Test case 2
    {
        Eigen::GMRES<MatrixXc, Eigen::IdentityPreconditioner> gmres;

        auto start = std::chrono::system_clock::now();

        gmres.compute(S);
        x = gmres.solve(b);

        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "GMRES with general matrix: " << elapsed.count() << '\n';

        std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error()
                  << std::endl;
        std::cout << (S * x - b).norm() << std::endl;
    }
}