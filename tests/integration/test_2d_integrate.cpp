//
// Created by evgen on 25.01.24.
//

#include "gtest/gtest.h"

#include "types/Types.hpp"
#include "math/integration/decart/Integration.hpp"


using EMW::Types::scalar;
using namespace EMW;
namespace GL = DecartIntegration::GaussLegendre;
namespace NC = DecartIntegration::NewtonCotess;

scalar f4(scalar x, scalar y) {
    return x * x * y * y * y * y;;
}

scalar f5(scalar x, scalar y) {
    return x * x * y * y * y * y * y;
}

scalar f(scalar x) {
    return (x + 1) * x * x * x * x;
}

class GAUSSIAN_QUADRATURE : public ::testing::Test {
protected:
    scalar prec = 4e-15;
};

TEST_F(GAUSSIAN_QUADRATURE, TEST_FOURTH_ORDER) {
    const auto result = DecartIntegration::integrate<GL::Quadrature<3, 3>>(f4, {0, 0}, {2, 2});
    ASSERT_EQ(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2. / 15, result);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_FIFTH_ORDER) {
    const auto result = DecartIntegration::integrate<GL::Quadrature<3, 3>>(f5, {0, 0}, {2, 2});
    ASSERT_NEAR(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2 * 2. / 18 - result, 0, prec);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_FIFTH_ORDER_1D) {
    const auto result = DecartIntegration::integrate<GL::Quadrature<3>>(f, {0}, {2});
    ASSERT_NEAR(2 * 2 * 2 * 2 * 2 * 2. / 6 + 2 * 2 * 2 * 2 * 2. / 5 - result, 0, prec);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_LINEAR_1) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DecartIntegration::integrate<GL::Quadrature<1, 1>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_LINEAR_4) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DecartIntegration::integrate<GL::Quadrature<1, 1>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_LINEAR_2) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x * y; };

    const auto result = DecartIntegration::integrate<GL::Quadrature<2, 2>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_LINEAR_3) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x + x * y; };

    const auto result = DecartIntegration::integrate<GL::Quadrature<1, 1>>(lambda, {0, 0}, {2, 2});
    ASSERT_EQ(8, result);
}


TEST_F(GAUSSIAN_QUADRATURE, TEST_SECOND_ORDER) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x * x * y * y; };

    const auto result = DecartIntegration::integrate<GL::Quadrature<2, 2>>(lambda, {1, 1}, {2, 2});
    ASSERT_EQ(26 * 26. / 9, result);
}


class NC_QUADRATURE : public ::testing::Test {
protected:
    scalar prec = 4e-15;
};

TEST_F(NC_QUADRATURE, TEST_LINEAR_3) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x + x * y; };

    const auto result = DecartIntegration::integrate<NC::Quadrature<1, 1>>(lambda, {0, 0}, {2, 2});
    ASSERT_EQ(8, result);
}


TEST_F(NC_QUADRATURE, TEST_LINEAR_1) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DecartIntegration::integrate<NC::Quadrature<1, 1>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(NC_QUADRATURE, TEST_LINEAR_4) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DecartIntegration::integrate<NC::Quadrature<1, 1>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(NC_QUADRATURE, TEST_LINEAR_2) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x * y; };

    const auto result = DecartIntegration::integrate<NC::Quadrature<4, 4>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

class ADAPTIVE_QUADRATURE : public ::testing::Test {

};

TEST_F(ADAPTIVE_QUADRATURE, LINEAR_TEST) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DecartIntegration::integrate_with_decomposition<NC::Quadrature<4, 4>>(
        lambda, std::make_tuple(-1., -1.), std::make_tuple(2., 2.), 3);
    ASSERT_NEAR(0, result, 1e-15);
    const auto result_1 = DecartIntegration::integrate_with_decomposition<GL::Quadrature<4, 4>>(
        lambda, std::make_tuple(-1., -1.), std::make_tuple(2., 2.), 3);
    ASSERT_NEAR(0, result_1, 1e-15);
}

TEST_F(ADAPTIVE_QUADRATURE, TEST_LINEAR_3) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x + x * y; };

    const auto result = DecartIntegration::integrate_with_decomposition<NC::Quadrature<1, 1>>(
        lambda, std::make_tuple(0, 0), std::make_tuple(2, 2), 1);
    ASSERT_EQ(8, result);
}


TEST_F(ADAPTIVE_QUADRATURE, LINEAR_ADAPTIVE_TEST) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x + x * y; };

    const auto result = DecartIntegration::adaptive_integrate<GL::Quadrature<4, 4>>(
        lambda, std::make_tuple(0., 0.), std::make_tuple(2., 2.),
        [](double old_res, double new_res) { return std::abs(old_res - new_res) < 1e-10 * new_res + 1e-15; });
    ASSERT_NEAR(8, result.first, 1e-15);
    std::cout << result.second << std::endl;
}
