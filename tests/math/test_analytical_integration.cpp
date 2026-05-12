//
// Created by evgen on 23.08.24.
//

#include "gtest/gtest.h"

#include "math/integration/analytical/SingularIntegration.hpp"

#include "math/integration/decart/Integration.hpp"

using namespace EMW;
using namespace EMW::Types;

class ANALYTICAL_INTEGRATION_TESTS : public testing::Test {};

TEST_F(ANALYTICAL_INTEGRATION_TESTS, ONE_DIV_R_OVER_CELL) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::vector points = {Mesh::point_t{-1, -1, 0}, Mesh::point_t{1, -1, 0}, Mesh::point_t{1, 1, 0},
                                       Mesh::point_t{-1, 1, 0}};
    Mesh::IndexedCell cell{{0, 1, 2, 3}, points};

    scalar res = Math::Integration::Analytical::integrate_1_div_r(Mesh::point_t{0, 0, 0}, cell);
    ASSERT_NEAR(std::abs((res - 8 * std::asinh(1.)) / res), 0, 1e-15);
}

TEST_F(ANALYTICAL_INTEGRATION_TESTS, ONE_DIV_R_OVER_CELL_BASE_ROUTINE) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::array verts = {Mesh::point_t{-1, -1, 0}, Mesh::point_t{1, -1, 0}, Mesh::point_t{1, 1, 0},
                                     Mesh::point_t{-1, 1, 0}};

    scalar res = Math::Integration::Analytical::integrate_1_div_r(point_t{0, 0, 0}, verts, Vector3d{0, 0, 1});
    ASSERT_NEAR(std::abs((res - 8 * std::asinh(1.)) / res), 0, 3e-16);
}

TEST_F(ANALYTICAL_INTEGRATION_TESTS, ONE_DIV_R_OVER_CELL_VIA_VERTEXES_ONLY) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::array verts = {Mesh::point_t{-1, -1, 0}, Mesh::point_t{1, -1, 0}, Mesh::point_t{1, 1, 0},
                                     Mesh::point_t{-1, 1, 0}};

    scalar res = Math::Integration::Analytical::integrate_1_div_r(point_t{0, 0, 0}, verts);
    ASSERT_NEAR(std::abs((res - 8 * std::asinh(1.)) / res), 0, 3e-16);
}

Types::scalar numerical_newton_potential(Types::point_t r, Types::scalar cl) {
    // интегрирование старое (квадратура поверх 1/r по квадрату)
    const auto cut_of_cube_integral = [cl](Types::scalar x, Types::scalar y, Types::scalar z, Types::scalar z_dash) {
        // Интегрирование по сечению куба на высоте z_dash
        Containers::vector<Types::point_t> points1 = {
            {-cl / 2., -cl / 2., z_dash},
            {cl / 2., -cl / 2., z_dash},
            {cl / 2., cl / 2., z_dash},
            {-cl / 2., cl / 2., z_dash},
        };

        return Math::Integration::Analytical::integrate_1_div_r(Types::point_t{x, y, z}, points1,
                                                                Types::Vector3d{0, 0, 1});
    };
    const auto integrand = [x = r.x(), y = r.y(), z = r.z(), cut_of_cube_integral](Types::scalar z_dash) {
        return cut_of_cube_integral(x, y, z, z_dash);
    };
    auto [result, level] = DecartIntegration::adaptive_integrate<DecartIntegration::GaussLegendre::Quadrature<7>>(
        integrand, {-cl / 2.}, {cl}, [](scalar x, scalar y) { return std::abs(x - y) < 1e-10 * std::abs(x) + 1e-16; },
        2048);
    // std::cout << level << std::endl;
    return result;
}

TEST_F(ANALYTICAL_INTEGRATION_TESTS, NEWTON_POTENTIAL_OF_CUBE) {
    constexpr Types::scalar cube_length = 1;
    constexpr size_t N = 300;
    constexpr double h = 0.05;
    std::cout.precision(std::numeric_limits<double>::digits10);
    for (size_t idx = 0; idx < N; ++idx) {
        Types::point_t r = {0, idx * h, idx * h};
        auto res1 = Math::Integration::Analytical::newtonian_potential_of_parallelepiped(
            r, cube_length / 2, cube_length / 2, cube_length / 2);
        // std::cout << "Distance = " << r.norm() << std::endl;
        auto res2 = numerical_newton_potential(r, cube_length);
        // std::cout << res1 << ' ' << res2 << std::endl;
        if (r.norm() < 0.77) {
            ASSERT_NEAR(std::abs(res1 - res2), 0, 1e-7 * std::abs(res1));
        } else {
            ASSERT_NEAR(std::abs(res1 - res2), 0, 1e-10 * std::abs(res1));
        }
    }
}

Types::scalar cube_newtonian_energy(Types::scalar l) {
    const auto integrand = [l](Types::scalar x, Types::scalar y, Types::scalar z) {
        return Math::Integration::Analytical::newtonian_potential_of_parallelepiped({x, y, z}, l/2, l/2, l/2);
    };
    auto [res, level] =  Math::Integration::Numerical::Decart::adaptive_integrate<
        DecartIntegration::GaussLegendre::Quadrature<3, 3, 3>>(
        integrand, {-l / 2., -l / 2., -l / 2.}, {l, l, l},
        [](scalar x, scalar y) { return std::abs(x - y) < 1e-9 * std::abs(x) + 1e-16; }, 128);
    // std::cout << level << std::endl;
    return res;
}

TEST_F(ANALYTICAL_INTEGRATION_TESTS, CUBE__NEWTONIAN_ENERGY) {
    constexpr size_t N = 300;
    constexpr double h = 0.05;
    std::cout.precision(std::numeric_limits<double>::digits10);
    for (size_t idx = 0; idx < N; ++idx) {
        auto cube_length = (idx + 1) * h;
        auto res1 = Math::Integration::Analytical::self_newtonian_energy_over_cube(cube_length);
        auto res2 = cube_newtonian_energy(cube_length);

        // std::cout << res1 << " " << res2 << std::endl;
        // std::cout << res2 / res1 << std::endl;

        ASSERT_NEAR(std::abs(res1 - res2), 0, 1e-8 * std::abs(res1));
    }
}