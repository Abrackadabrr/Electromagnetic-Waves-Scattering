//
// Created by evgen on 27.01.2026.
//

#include "gtest/gtest.h"

#include "math/integration/Quadrature.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/integration/newton_cotess/Rectangular.hpp"
#include "types/Types.hpp"

using EMW::Types::scalar;
using namespace EMW;
namespace GL = DefiniteIntegrals::GaussLegendre;
namespace NC = DefiniteIntegrals::NewtonCotess;

scalar linear_x(scalar x, scalar y, scalar z) { return x; }
scalar linear_y(scalar x, scalar y, scalar z) { return y; }
scalar linear_z(scalar x, scalar y, scalar z) { return z; }

scalar f_pow(scalar x, scalar y, scalar z, scalar p, scalar q, scalar m) {
    return std::pow(x, p) * std::pow(y, q) + std::pow(z, m);
}

class GAUSSIAN_QUADRATURE : public ::testing::Test {
  protected:
    scalar prec = 4e-15;
};

TEST_F(GAUSSIAN_QUADRATURE, TEST_LINEAR) {
    const auto res_x = DefiniteIntegrals::integrate<GL::Quadrature<1, 1, 1>>(linear_x, {0, 0, 0}, {1, 1, 1});
    const auto res_y = DefiniteIntegrals::integrate<GL::Quadrature<1, 1, 1>>(linear_y, {0, 0, 0}, {1, 1, 1});
    const auto res_z = DefiniteIntegrals::integrate<GL::Quadrature<1, 1, 1>>(linear_z, {0, 0, 0}, {1, 1, 1});

    ASSERT_NEAR(res_x, 1. / 2, prec);
    ASSERT_NEAR(res_z, 1. / 2, prec);
    ASSERT_NEAR(res_y, 1. / 2, prec);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_2ND_ORDER) {
    const auto fx = [](scalar x, scalar y, scalar z) { return x * y; };
    const auto fy = [](scalar x, scalar y, scalar z) { return x * z; };
    const auto fz = [](scalar x, scalar y, scalar z) { return z * z; };
    const auto res_x = DefiniteIntegrals::integrate<GL::Quadrature<2, 2, 2>>(fx, {0, 0, 0}, {1, 1, 1});
    const auto res_y = DefiniteIntegrals::integrate<GL::Quadrature<2, 2, 2>>(fy, {0, 0, 0}, {1, 1, 1});
    const auto res_z = DefiniteIntegrals::integrate<GL::Quadrature<2, 2, 2>>(fz, {0, 0, 0}, {1, 1, 1});

    ASSERT_NEAR(res_x, 1. / 4, prec);
    ASSERT_NEAR(res_z, 1. / 3, prec);
    ASSERT_NEAR(res_y, 1. / 4, prec);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_3RD_ORDER) {
    const auto fx = [](scalar x, scalar y, scalar z) { return f_pow(x, y, z, 1, 2, 0); };
    const auto fy = [](scalar x, scalar y, scalar z) { return f_pow(x, y, z, 1, 0, 2); };
    const auto fz = [](scalar x, scalar y, scalar z) { return f_pow(x, y, z, 0, 0, 3); };
    const auto res_x = DefiniteIntegrals::integrate<GL::Quadrature<3, 3, 3>>(fx, {0, 0, 0}, {1, 1, 1});
    const auto res_y = DefiniteIntegrals::integrate<GL::Quadrature<3, 3, 3>>(fy, {0, 0, 0}, {1, 1, 1});
    const auto res_z = DefiniteIntegrals::integrate<GL::Quadrature<3, 3, 3>>(fz, {0, 0, 0}, {1, 1, 1});

    ASSERT_NEAR(res_x, 1. / 6, prec);
    ASSERT_NEAR(res_y, 1. / 6, prec);
    ASSERT_NEAR(res_z, 1. / 4, prec);
}

TEST_F(GAUSSIAN_QUADRATURE, TEST_GRAVITATIONAL_ENERGY) {
    const Types::scalar l = 3;
    // аналитика
    const Types::scalar an_res = Math::AnalyticalIntegration::self_newtonian_energy_over_cube(l);

    // интегрирование новое
    const auto integrand = [l](scalar x, scalar y, scalar z) {
        const Types::point_t point{x, y, z};
        return Math::AnalyticalIntegration::newtonian_potential_of_parallelepiped(point, l / 2., l / 2., l / 2.);
    };
    const auto res =
        DefiniteIntegrals::integrate<GL::Quadrature<5, 5, 5>>(integrand, {-l / 2, -l / 2, -l / 2}, {l, l, l});

    // интегрирование старое (квадратура в квадратуре в квадратуре)
    const auto cut_of_cube_integral = [l](scalar x, scalar y, scalar z, scalar z_dash) {
        // Интегрирование по сечениям куба
        const Containers::vector<Types::point_t> points1 = {
            {-l / 2, -l / 2, z_dash},
            {l / 2, -l / 2, z_dash},
            {l / 2, l / 2, z_dash},
            {-l / 2, l / 2, z_dash},
        };
        const Mesh::IndexedCell cell1({0, 1, 2, 3}, points1);
        return Math::AnalyticalIntegration::integrate_1_div_r(Types::point_t{x, y, z}, cell1);
    };
    const auto newtonian_potential_of_cube = [l, cut_of_cube_integral](scalar x, scalar y, scalar z) {
        const auto integrand = [x, y, z, cut_of_cube_integral](Types::scalar z_dash) {
            return cut_of_cube_integral(x, y, z, z_dash);
        };
        return DefiniteIntegrals::integrate<GL::Quadrature<4>>(integrand, {-l/2}, {l/2}) +
            DefiniteIntegrals::integrate<GL::Quadrature<4>>(integrand, {0}, {l/2});
    };
    const auto res_old =
        DefiniteIntegrals::integrate<GL::Quadrature<7, 7, 7>>(newtonian_potential_of_cube, {-l / 2, -l / 2, -l / 2}, {l, l, l});


    std::cout.precision(16);
    std::cout << res << std::endl;
    std::cout << an_res << std::endl;
    std::cout << res_old << std::endl;
}
