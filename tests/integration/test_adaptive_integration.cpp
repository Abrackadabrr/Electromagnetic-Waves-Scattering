//
// Created by evgen on 06.05.2026.
//

#include <gtest/gtest.h>

#include "types/Types.hpp"
#include "math/integration/decart/Integration.hpp"


using EMW::Types::scalar;
using namespace EMW;
namespace GL = DecartIntegration::GaussLegendre;
namespace NC = DecartIntegration::NewtonCotess;

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

    for (int level = 1; level <= 3; ++level) {
        const auto result = DecartIntegration::integrate_with_decomposition<NC::Quadrature<1, 1>>(
            lambda, std::make_tuple(0, 0), std::make_tuple(2, 2), level);
        ASSERT_EQ(8, result);
    }
}


TEST_F(ADAPTIVE_QUADRATURE, LINEAR_ADAPTIVE_TEST) {
    constexpr auto lambda = [&](const scalar x, const scalar y) { return x + x * y; };

    constexpr scalar rTol = 1e-10;
    constexpr scalar aTol = 1e-10;

    const auto result = DecartIntegration::adaptive_integrate<GL::Quadrature<4, 4>>(
        lambda, {0., 0.}, {2., 2.},
        [rTol, aTol](double old_res, double new_res) {
            return std::abs(old_res - new_res) < rTol * new_res + aTol;
        }, 3);
    ASSERT_NEAR(8, result.first, rTol);
}

using namespace std::complex_literals;

TEST_F(ADAPTIVE_QUADRATURE, EXPONENT_ADAPTIVE_TEST) {
    constexpr scalar k = 1;
    const auto lambda = [&k](const scalar x, const scalar y, const scalar z) { return std::exp(1i * k * (x + y + z)); };
    const auto single_res = (std::exp(k * 1i) - 1.) / (1i * k);
    constexpr scalar rTol = 1e-7;
    constexpr scalar aTol = 1e-15;

    const auto result = DecartIntegration::adaptive_integrate<GL::Quadrature<3, 3, 3>>(
        lambda, {0., 0., 0}, {1., 1., 1},
        [rTol, aTol](Types::complex_d old_res, Types::complex_d new_res) {
            return std::abs(old_res - new_res) < rTol * std::abs(new_res) + aTol;
        }, 7);
    std::cout << result.second << std::endl;
    ASSERT_NEAR(std::abs(single_res * single_res * single_res - result.first), 0, rTol);
}

TEST_F(ADAPTIVE_QUADRATURE, EXPONENT_ADAPTIVE_TEST_TWODIM) {
    constexpr scalar k = 1;
    const auto lambda = [&k](const scalar x, const scalar y) { return std::exp(1i * k * (x + y)); };
    const auto single_res = (std::exp(k * 1i) - 1.) / (1i * k);
    constexpr scalar rTol = 1e-7;
    constexpr scalar aTol = 1e-15;

    const auto result = DecartIntegration::adaptive_integrate<GL::Quadrature<3, 3>>(
        lambda, {0., 0.}, {1., 1.},
        [rTol, aTol](Types::complex_d old_res, Types::complex_d new_res) {
            return std::abs(old_res - new_res) < rTol * std::abs(new_res) + aTol;
        }, 7);
    std::cout << result.second << std::endl;
    ASSERT_NEAR(std::abs(single_res * single_res - result.first), 0, rTol);
}


TEST_F(ADAPTIVE_QUADRATURE, EXPONENT_ADAPTIVE_TEST_ONEDIM) {
    constexpr scalar k = 1;
    const auto lambda = [&k](const scalar x) { return std::exp(1i * k * x); };
    const auto single_res = (std::exp(k * 1i) - 1.) / (1i * k);
    constexpr scalar rTol = 1e-7;
    constexpr scalar aTol = 1e-15;

    const auto result = DecartIntegration::adaptive_integrate<GL::Quadrature<3>>(
        lambda, {0.}, {1.},
        [rTol, aTol](Types::complex_d old_res, Types::complex_d new_res) {
            return std::abs(old_res - new_res) < rTol * std::abs(new_res) + aTol;
        }, 7);
    std::cout << result.second << std::endl;
    ASSERT_NEAR(std::abs(single_res - result.first), 0, rTol);
}
