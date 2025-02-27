//
// Created by evgen on 26.02.2025.
//

#include "math/matrix/Matrix.hpp"

#include "toeplitz_matrix_tests.hpp"
#include "gtest/gtest.h"
#include <unsupported/Eigen/IterativeSolvers>

Types::scalar simple_toeplitz_matrix(Types::index i, Types::index j) {
    return static_cast<integer>(i) - static_cast<integer>(j);
};

Types::MatrixX<scalar> get_block(Types::index i, Types::index j) {
    return Types::MatrixX<scalar>{{{1, 2}, {2, 3}}} * simple_toeplitz_matrix(i, j);
}

Types::MatrixX<scalar> get_toeplitz_matrix(Types::index rows, Types::index cols) {
    Types::MatrixX<scalar> result = Types::MatrixX<scalar>::Zero(2 * rows, 2 * cols);
    for (Types::index i = 0; i < rows; i++)
        for (Types::index j = 0; j < cols; j++)
            result.block(2 * i, 2 * j, 2, 2) = get_block(i, j);
    return result;
}

using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzBlock<scalar>;
using block_t = Math::LinAgl::Matrix::ToeplitzBlock<scalar>::block_type;

template <typename MatrixType>
Types::VectorXc solve(const MatrixType &A, const Types::VectorXc &b, Types::scalar tolerance) {

    auto method = Eigen::GMRES<MatrixType>{};
    Types::index max_iterations = 1000;
    // use matrix-free context of gmras method
    method.setMaxIterations(max_iterations);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(tolerance);
    method.set_restart(max_iterations);

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    method.compute(A);
    auto j_vec = Types::VectorXc{method.solve(b)};

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::ratio<1, 1>> elapsed_seconds = end - start;
    std::cout << "Время решения GMRES: " << elapsed_seconds.count() << std::endl;
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "tolerance: " << method.tolerance() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    return j_vec;
}

TEST_F(TOEPLITZ_MATRIX_TESTS, TOEPLITZ_BLOCK_IN_EIGEN_SOLVER) {
    // Описываем структуру матрицы и создаем первую
    Types::index N = 1500;
    Types::index M = 1400;
    Types::scalar value = 10;
    // Делаем обычную тёплицеву матрицу
    const auto &full_matrix = get_toeplitz_matrix(N, M);
    // Создаем матрицу с тёплицевой структурой
    Containers::vector<block_t> blocks;
    Types::index size = ToeplitzBlock::get_size_of_container(N, M);
    blocks.reserve(size);
    for (int i = 0; i < M; i++) {
        blocks.emplace_back(get_block(0, i));
    }
    for (int i = 1; i < N; i++) {
        blocks.emplace_back(get_block(i, 0));
    }
    ToeplitzBlock test_matrix(N, M, std::move(blocks));
    ASSERT_EQ(blocks.size(), 0);
    ASSERT_EQ(test_matrix.rows(), N * 2);
    ASSERT_EQ(test_matrix.cols(), M * 2);
    // Создаем вектор и умножаемся
    Types::VectorX<scalar> vec = Types::VectorX<scalar>::Zero(2 * M);
    for (Types::index i = 0; i < M; i++)
        vec(i) = static_cast<scalar>(i);

    // Собираем вектор и умножаем
    ASSERT_NEAR((test_matrix * vec - full_matrix * vec).norm(), 0, 1e-14);

    constexpr Types::scalar tolerance = 1e-2;
    // const auto solution_with_full_matrix = solve(full_matrix, full_matrix * vec, tolerance);
    // const auto solution_with_test_matrix = solve(test_matrix, vec, tolerance);

    //ASSERT_NEAR((solution_with_full_matrix - solution_with_test_matrix).norm() / solution_with_test_matrix.norm(), 0, 1e-14);
}
