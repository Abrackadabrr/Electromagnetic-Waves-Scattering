//
// Created by evgen on 25.01.24.
//

#include "gtest/gtest.h"

#include "types/Types.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/integration/newton_cotess/Rectangular.hpp"
#include "math/integration/Quadrature.hpp"

using EMW::Types::scalar;
using namespace EMW;
namespace GL = DefiniteIntegrals::GaussLegendre;
namespace NC = DefiniteIntegrals::NewtonCotess;

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
    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<3, 3>>(f4, {0, 0}, {2, 2});
    ASSERT_EQ(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2. / 15, result);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_FIFTH_ORDER) {
    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<3, 3>>(f5, {0, 0}, {2, 2});
    ASSERT_NEAR(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2 * 2. / 18 - result, 0, prec);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_FIFTH_ORDER_1D) {
    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<3>>(f, {0}, {2});
    ASSERT_NEAR(2 * 2 * 2 * 2 * 2 * 2. / 6 + 2 * 2 * 2 * 2 * 2. / 5 - result, 0, prec);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_LINEAR_1) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<1, 1>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_LINEAR_4) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<1, 1>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_LINEAR_2) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x*y; };

    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<2, 2>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_LINEAR_3) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x + x*y; };

    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<1, 1>>(lambda, {0, 0}, {2, 2});
    ASSERT_EQ(8, result);
}


TEST_F(GAUSSIAN_QUADRATURE, TEST_SECOND_ORDER) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x*x * y*y; };

    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<2, 2>>(lambda, {1, 1}, {2, 2});
    ASSERT_EQ(26 * 26. / 9, result);
}


class NC_QUADRATURE : public ::testing::Test {
protected:
    scalar prec = 4e-15;
};

TEST_F(NC_QUADRATURE, TEST_LINEAR_3) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x + x*y; };

    const auto result = DefiniteIntegrals::integrate<NC::Quadrature<1, 1>>(lambda, {0, 0}, {2, 2});
    ASSERT_EQ(8, result);
}


TEST_F(NC_QUADRATURE, TEST_LINEAR_1) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DefiniteIntegrals::integrate<NC::Quadrature<1, 1>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(NC_QUADRATURE, TEST_LINEAR_4) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DefiniteIntegrals::integrate<NC::Quadrature<1, 1>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(NC_QUADRATURE, TEST_LINEAR_2) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x*y; };

    const auto result = DefiniteIntegrals::integrate<NC::Quadrature<4, 4>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}
