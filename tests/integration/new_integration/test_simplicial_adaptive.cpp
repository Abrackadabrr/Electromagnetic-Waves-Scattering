//
// Created by evgen on 18.05.2026.
//

#include <gtest/gtest.h>

#include "math/integration/simplicial/Integration.hpp"
#include "types/Types.hpp"

using EMW::Types::scalar;
using namespace EMW;

class SIMPLICIAL_ADAPTIVE_QUADRATURE : public ::testing::Test {};

TEST_F(SIMPLICIAL_ADAPTIVE_QUADRATURE, LINEAR_TEST) {
    constexpr auto lambda = [&](const Types::point_t &r) { return r.x(); };

    const Containers::array verts = {Types::point_t{0, 0, 0}, Types::point_t{1, 0, 0}, Types::point_t{0, 1, 0}};

    const auto result =
        0.5 *
        SimplicialIntegration::triangle_quadrature_sum_with_decomposition<SimplicialIntegration::GaussianPoints<2, 5>>(
            lambda, 2, verts[0], verts[1], verts[2]);
    ASSERT_NEAR(1. / 6., result, 1e-15);
}

TEST_F(SIMPLICIAL_ADAPTIVE_QUADRATURE, TEST_LINEAR_3) {
    constexpr auto lambda = [&](const Types::point_t &r) { return r.x() + r.y() * r.x(); };

    const Containers::array verts = {Types::point_t{0, 0, 0}, Types::point_t{1, 0, 0}, Types::point_t{0, 1, 0}};

    for (int level = 0; level <= 3; ++level) {
        const auto result =
            0.5 * SimplicialIntegration::triangle_quadrature_sum_with_decomposition<
                      SimplicialIntegration::GaussianPoints<2, 5>>(lambda, 1 << level, verts[0], verts[1], verts[2]);
        ASSERT_NEAR(5. / 24., result, 1e-15);
    }
}

TEST_F(SIMPLICIAL_ADAPTIVE_QUADRATURE, LINEAR_ADAPTIVE_TEST) {
    constexpr auto lambda = [&](const Types::point_t &r) { return r.x() + r.y() * r.x(); };
    constexpr scalar rTol = 1e-10;
    constexpr scalar aTol = 1e-10;
    // трекгольник
    const Containers::array verts = {Types::point_t{0, 0, 0}, Types::point_t{1, 0, 0}, Types::point_t{0, 1, 0}};

    const auto result = SimplicialIntegration::adaptive_quadrature_sum<SimplicialIntegration::GaussianPoints<2, 5>>(
        lambda,
        [rTol, aTol](double old_res, double new_res) { return std::abs(old_res - new_res) < rTol * new_res + aTol; },
        verts[0], verts[1], verts[2], 4);
    std::cout << result.level << std::endl;
    ASSERT_NEAR(5. / 24., 0.5 * result.value, rTol);
}

using namespace std::complex_literals;

TEST_F(SIMPLICIAL_ADAPTIVE_QUADRATURE, EXPONENT_ADAPTIVE_TEST) {
    constexpr scalar k = 10;
    const auto lambda = [&k](const Types::point_t &p) { return std::exp(1i * k * (p.x() + p.y())); };
    const auto single_res = (std::exp(k * 1i) - 1.) / (1i * k);
    constexpr scalar rTol = 1e-7;
    constexpr scalar aTol = 1e-15;

    const Containers::array verts = {Types::point_t{0, 0, 0}, Types::point_t{1, 0, 0}, Types::point_t{1, 1, 0},
                                     Types::point_t{0, 1, 0}};
    auto result1 = SimplicialIntegration::adaptive_quadrature_sum<SimplicialIntegration::GaussianPoints<2, 5>>(
        lambda,
        [rTol, aTol](Types::complex_d old_res, Types::complex_d new_res) {
            return std::abs(old_res - new_res) < rTol * std::abs(new_res) + aTol;
        },
        verts[0], verts[1], verts[3], 7);

    auto result2 = SimplicialIntegration::adaptive_quadrature_sum<SimplicialIntegration::GaussianPoints<2, 5>>(
                        lambda,
                        [rTol, aTol](Types::complex_d old_res, Types::complex_d new_res) {
                            return std::abs(old_res - new_res) < rTol * std::abs(new_res) + aTol;
                        },
                        verts[1], verts[2], verts[3], 7);

    std::cout << result1.level << std::endl;
    std::cout << result2.level << std::endl;
    ASSERT_NEAR(std::abs(single_res * single_res - 0.5 * (result1.value + result2.value)), 0, rTol);
}
