//
// Created by evgen on 18.03.2025.
//

#include <gtest/gtest.h>

#include "math/matrix/DynamicFactoredMatrix.hpp"
#include "types/Types.hpp"

using namespace EMW;

class DYNAMIC_FACTORED_MATRIX_TEST : public ::testing::Test {
  protected:
    using f_mat = Types::MatrixXd;
    using c_mat = Types::MatrixXc;
};

TEST_F(DYNAMIC_FACTORED_MATRIX_TEST, CREATION) {
    int n = 10;
    const f_mat identity = f_mat::Identity(n, n);
    const f_mat ones = f_mat::Ones(n, n);
    const Types::VectorXd vector = Types::VectorXd::Random(n);
    Math::LinAgl::Matrix::DynamicFactoredMatrix<f_mat> matrix1{{f_mat::Identity(n, n), f_mat::Identity(n, n), f_mat::Identity(n, n)}};
    Math::LinAgl::Matrix::DynamicFactoredMatrix<f_mat> matrix2{{f_mat::Identity(n, n), 2 * f_mat::Identity(n, n), 4 * f_mat::Identity(n, n)}};
    Math::LinAgl::Matrix::DynamicFactoredMatrix<f_mat> matrix3{{f_mat::Ones(n, n), f_mat::Identity(n, n), 7 * f_mat::Identity(n, n)}};

    ASSERT_EQ(identity, matrix1.get<0>());
    ASSERT_EQ(identity, matrix1.get<1>());
    ASSERT_EQ(identity, matrix1.get<2>());

    ASSERT_EQ(identity, matrix2.get<0>());
    ASSERT_EQ(2 * identity, matrix2.get<1>());
    ASSERT_EQ(4 * identity, matrix2.get<2>());

    ASSERT_EQ(ones, matrix3.get<0>());
    ASSERT_EQ(identity, matrix3.get<1>());
    ASSERT_EQ(7 * identity, matrix3.get<2>());
}

TEST_F(DYNAMIC_FACTORED_MATRIX_TEST, MANIPULATION) {
    int n = 10;
    const f_mat identity = f_mat::Identity(n, n);
    const f_mat ones = f_mat::Ones(n, n);
    const Types::VectorXd vector = Types::VectorXd::Random(n);
    const Math::LinAgl::Matrix::DynamicFactoredMatrix<f_mat> matrix1{{identity, identity}};
    const Math::LinAgl::Matrix::DynamicFactoredMatrix<f_mat> matrix2{{identity, 2 * identity, 4 * identity}};
    const Math::LinAgl::Matrix::DynamicFactoredMatrix<f_mat> matrix3{{ones, identity, 7 * identity}};

    // Проверка селекторов
    ASSERT_EQ(matrix1.factor_number(), 2);
    ASSERT_EQ(matrix2.factor_number(), 3);
    ASSERT_EQ(matrix3.factor_number(), 3);

    ASSERT_EQ(matrix1.rows(), n);
    ASSERT_EQ(matrix2.rows(), n);
    ASSERT_EQ(matrix3.rows(), n);

    ASSERT_EQ(matrix1.cols(), n);
    ASSERT_EQ(matrix2.cols(), n);
    ASSERT_EQ(matrix3.cols(), n);

    // Проверка матвека

    const auto res1 = matrix1 * vector;
    const auto res2 = matrix2 * vector;
    const auto res3 = matrix3 * vector;

    const Types::VectorXd& right_res1 = vector;
    const Types::VectorXd right_res2 = 2 * (4 * vector);
    const Types::VectorXd right_res3 = ones * ( identity * ((7 * identity) * vector));

    ASSERT_EQ(right_res1, res1);
    EXPECT_EQ(right_res2, res2);
    EXPECT_EQ(right_res3, res3);
}
