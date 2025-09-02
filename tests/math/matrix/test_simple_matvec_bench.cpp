//
// Created by evgen on 01.09.2025.
//

#include <gtest/gtest.h>
#include <eigen3/Eigen/Dense>

#include <chrono>

TEST(MATVEC, RESHAPE_TEST) {
    const int N = 4;
    const int K = 3;
    const Eigen::MatrixXd A = Eigen::MatrixXd::Random(N, N);
    const Eigen::VectorXd B = Eigen::VectorXd::Random(K * N);

    std::cout << "Vec: " << B.transpose() << std::endl;
    std::cout << B.reshaped(N, K).reshaped(N * K, 1).transpose() << std::endl;
}

TEST(MATVEC, SIMPLE_TEST) {
    const int N = 3000;
    const int K = 40;
    const Eigen::MatrixXd A = Eigen::MatrixXd::Random(N, N);
    const Eigen::VectorXd B = Eigen::VectorXd::Random(K * N);

    // 1) Прямое перемножение матрицы на блоки вектора

    auto start = std::chrono::high_resolution_clock::now();

    Eigen::VectorXd res_1 = Eigen::VectorXd::Zero(K * N);

    for (int i = 0; i < K; i++) {
        res_1.block(i * N, 0, N, 1) = A * B.block(i * N, 0, N, 1);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_1 = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // 2) Решейп в матрицу и перемножение
    start = std::chrono::high_resolution_clock::now();

    const Eigen::VectorXd res_2 = (A * B.reshaped(N, K)).reshaped(N * K, 1);

    end = std::chrono::high_resolution_clock::now();
    auto elapsed_2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // 3) Смотрим на одинаковость результатов
    const auto error = (res_1 - res_2).norm();
    std::cout << "Elapsed time 1 : " << elapsed_1.count() << " us" << std::endl;
    std::cout << "Elapsed time 2 : " << elapsed_2.count() << " us" << std::endl;
    std::cout << "Error : " << error << std::endl;
}