//
// Created by evgen on 26.02.2025.
//

#include "math/matrix/Matrix.hpp"
#include "math/matrix/eigen_matrix_replacement/MatrixReplacement.hpp"
#include "math/matrix/eigen_matrix_replacement/MatrixTraits.hpp"

#include "toeplitz_matrix_tests.hpp"

#include <unsupported/Eigen/IterativeSolvers>

#include "gtest/gtest.h"


using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzBlock<complex_d>;
using block_t = ToeplitzBlock::block_type;

template <typename MatrixType>
Types::VectorXc solve(const MatrixType &A, const Types::VectorXc &b, Types::scalar tolerance) {

    Eigen::GMRES<MatrixType, Eigen::IdentityPreconditioner> method;

    method.setTolerance(tolerance);
    VectorXc x{};

    auto start = std::chrono::high_resolution_clock::now();

    method.compute(A);
    x = method.solve(b);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Время решения GMRES: " << elapsed_seconds.count() << std::endl;
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "tolerance: " << method.tolerance() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    return x;
}

TEST_F(TOEPLITZ_MATRIX_TESTS, TOEPLITZ_BLOCK_IN_EIGEN_SOLVER) {
    // Описываем структуру матрицы и создаем первую
    Types::index tl_rows = 16;
    Types::index tl_cols = 16;
    Types::index block_size =  1000;
    Types::index N = block_size * tl_rows;
    Types::index M = block_size * tl_cols;
    const Math::LinAgl::Matrix::ToeplitzBlock<complex_d> toeplitz_matrix(
        tl_rows, tl_cols, [block_size, this](Types::index i, Types::index j) { return get_block(i, j, block_size); });

    const auto &full_matrix = get_toeplitz_matrix(tl_rows, tl_cols, block_size);
    ASSERT_EQ(toeplitz_matrix.rows(), N);
    ASSERT_EQ(toeplitz_matrix.cols(), M);
    // Создаем вектор и умножаемся
    EMW::Types::VectorXc vec(M);
    vec.setRandom();
    // Собираем вектор и умножаем
    ASSERT_NEAR((toeplitz_matrix * vec - full_matrix * vec).norm() / (full_matrix * vec).norm(), 0, 1e-14);

    constexpr Types::scalar tolerance = 1e-2;
    std::cout << "With normal matrix" << std::endl;
    const auto solution_with_full_matrix = solve(full_matrix, full_matrix * vec, tolerance);
    std::cout << "With toeplitz matrix" << std::endl;
    const auto solution_with_test_matrix = solve(Wrappers::MatrixReplacement{toeplitz_matrix}, toeplitz_matrix * vec, tolerance);

    ASSERT_NEAR((solution_with_full_matrix - solution_with_test_matrix).norm() / solution_with_test_matrix.norm(), 0, 1e-14);
}
