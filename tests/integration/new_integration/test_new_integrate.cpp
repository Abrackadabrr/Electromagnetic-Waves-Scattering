//
// Created by evgen on 25.01.24.
//

#include "gtest/gtest.h"

#include "types/Types.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/integration/newton_cotess/Rectangular.hpp"
#include "math/integration/Quadrature.hpp"

#include "EMW/math/integration/new_numerical_quadratures/base_routines/GeneralQuadrature.hpp"
#include "EMW/math/integration/new_numerical_quadratures/base_routines/gauss_quadrature/SegmentQuadraturePoints.hpp"
#include "EMW/math/integration/new_numerical_quadratures/base_routines/gauss_quadrature/TriangleQuadraturePoints.hpp"

#include <mesh/MeshTypes.hpp>

using EMW::Types::scalar;
using namespace EMW;
namespace GL = DefiniteIntegrals::GaussLegendre;
namespace NC = DefiniteIntegrals::NewtonCotess;
namespace NewGL = Math::Integration::Numerical::GaussianQuadratures;

scalar f4(scalar x, scalar y) {
    return x * x * y * y * y * y;;
}

scalar f5(scalar x, scalar y) {
    return x * x * y * y * y * y * y;
}

constexpr scalar f(scalar x) {
    return (x + 1) * x * x * x * x;
}

class NEW_GAUSSIAN_QUADRATURE : public ::testing::Test {
protected:
    scalar prec = 4e-15;
};

struct Segment {
    scalar a, b;

    constexpr scalar barycentric(scalar w) const {
        return w * a + (1 - w) * b;
    }

    constexpr scalar mes() const {
        return b - a;
    }

    using point_t = scalar;
};

struct Triangle {
    using point_t = Types::point_t;
    point_t a, b, c;

    constexpr scalar mes() const {
        return (b - a).cross(b - c).norm() / 2;
    }

    constexpr point_t barycentric(const Containers::array<scalar, 2> &coord) const {
        return coord[0] * a + coord[1] * b + (1 - coord[0] - coord[1]) * c;
    }

    constexpr point_t barycentric(scalar u, scalar v) const {
        return u * a + v * b + (1 - v - u) * c;
    }
};

TEST_F(NEW_GAUSSIAN_QUADRATURE, TEST_SEGMENT) {
    constexpr double result = Math::Integration::Numerical::integrate<NewGL::GaussianPoints<1, 11>>(f, Segment{0, 1});
    constexpr double result2 = Math::Integration::Numerical::quadrature_sum<NewGL::GaussianPoints<
        1, 11>>(f, 0., 1.);
    ASSERT_NEAR(11./30, result, prec);
    ASSERT_NEAR(11./30, result2, prec);
}

TEST_F(NEW_GAUSSIAN_QUADRATURE, TEST_TRIANGLE) {
    const auto integrand = [](const Types::point_t &x) { return 0.1; };
    const double result = Math::Integration::Numerical::integrate<NewGL::GaussianPoints<2, 5>>(
        integrand, Triangle{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}});
    const double result2 = Math::Integration::Numerical::quadrature_sum<NewGL::GaussianPoints<2, 5>>(
        integrand, Types::point_t{0, 0, 0}, Types::point_t{1, 0, 0}, Types::point_t{0, 1, 0});

    ASSERT_NEAR(0.05, result, prec);
    ASSERT_NEAR(0.1, result2, prec);
    std::cout.precision(20);
    std::cout << result2 << std::endl;
    std::cout << 0.1 << std::endl;
    std::cout << std::nextafter(0.1, 1)<< std::endl;
}

TEST_F(NEW_GAUSSIAN_QUADRATURE, TEST_FIFTH_ORDER) {
    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<3, 3>>(f5, {0, 0}, {2, 2});
    ASSERT_NEAR(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2 * 2. / 18 - result, 0, prec);
}

TEST_F(NEW_GAUSSIAN_QUADRATURE, TEST_FIFTH_ORDER_1D) {
    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<3>>(f, {0}, {2});
    ASSERT_NEAR(2 * 2 * 2 * 2 * 2 * 2. / 6 + 2 * 2 * 2 * 2 * 2. / 5 - result, 0, prec);
}

TEST_F(NEW_GAUSSIAN_QUADRATURE, TEST_LINEAR_1) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<1, 1>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(NEW_GAUSSIAN_QUADRATURE, TEST_LINEAR_4) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x; };

    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<1, 1>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(NEW_GAUSSIAN_QUADRATURE, TEST_LINEAR_2) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x * y; };

    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<2, 2>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}

TEST_F(NEW_GAUSSIAN_QUADRATURE, TEST_LINEAR_3) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x + x * y; };

    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<1, 1>>(lambda, {0, 0}, {2, 2});
    ASSERT_EQ(8, result);
}

TEST_F(NEW_GAUSSIAN_QUADRATURE, TEST_SECOND_ORDER) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x * x * y * y; };

    const auto result = DefiniteIntegrals::integrate<GL::Quadrature<2, 2>>(lambda, {1, 1}, {2, 2});
    ASSERT_EQ(26 * 26. / 9, result);
}

class NC_QUADRATURE : public ::testing::Test {
protected:
    scalar prec = 4e-15;
};

TEST_F(NC_QUADRATURE, TEST_LINEAR_3) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x + x * y; };

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
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x * y; };

    const auto result = DefiniteIntegrals::integrate<NC::Quadrature<4, 4>>(lambda, {-1, -1}, {2, 2});
    ASSERT_EQ(0, result);
}
