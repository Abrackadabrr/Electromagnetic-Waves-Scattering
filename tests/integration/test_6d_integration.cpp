//
// Created by evgen on 27.01.2026.
//
#include "gtest/gtest.h"

#include "math/integration/Quadrature.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"

#include "types/Types.hpp"

#include "operators/Functions.hpp"
#include "operators/volume/OperatorK.hpp"

#include <numeric>
#include <ranges>

using EMW::Types::scalar;
using namespace EMW;
namespace GL = DefiniteIntegrals::GaussLegendre;

class MORE_THAN_3D_INTEGRATION : public ::testing::Test {
  protected:
    scalar prec = 4e-15;
};

/**
 * Тест на интегрирование по квадратуре Гаусса и интегрирование по методу средней точки (одной единственной)
 * объемной части оператора K.
 *
 * Интересно, что при стремлении расстояния между объемами к бесконечности, ошибка численного интегрирования стабилизируется
 * на каком-то одном уровне и перестает падать, в противопоставление с аналогичным тестом для k = 0.
 */
TEST_F(MORE_THAN_3D_INTEGRATION, TEST_GRAVITATIONAL_ENERGY) {
    const Types::scalar l = 0.1;
    const Types::scalar l2 = l * l;
    Containers::vector<Types::scalar> distances; distances.resize(30);
    std::iota(distances.begin(), distances.end(), 0);
    for (auto &&distance : distances) {
        // Создаём два куба с заданным размером на заданном расстоянии друг от друга
        const Containers::vector<Types::point_t> points1 = {
            {-l / 2, -l / 2, -l / 2}, {l / 2, -l / 2, -l / 2}, {-l / 2, l / 2, -l / 2}, {l / 2, l / 2, -l / 2},
            {-l / 2, -l / 2, l / 2},  {l / 2, -l / 2, l / 2},  {-l / 2, l / 2, l / 2},  {l / 2, l / 2, l / 2},
        };
        Containers::vector<Types::point_t> points2;
        points2.resize(points1.size());
        std::transform(points1.begin(), points1.end(), points2.begin(), [distance](const Types::point_t &point) {
            return Types::point_t{point.x(), point.y(), point.z() + distance};
        });

        const Mesh::VolumeCells::IndexedCube cube1(points1, {0, 1, 2, 3, 4, 5, 6, 7});
        const Mesh::VolumeCells::IndexedCube cube2(points2, {0, 1, 2, 3, 4, 5, 6, 7});

        // Счёт интегралов по разным формулам
        // 1. Очень простой расчет по одной точке
        const Types::complex_d k = {2 * 3.1415926 / 10, 0};
        const auto far_zone_result = l2 * l2 * l2 * Helmholtz::F_bounded_part(k, cube1.center_, cube2.center_);

        // 2. Очень сложный расчет по шестимерной квадратуре
        const auto integrand = [k](Types::scalar x1, Types::scalar y1, Types::scalar z1, Types::scalar x2, Types::scalar y2,
                                   Types::scalar z2) {
            return Helmholtz::F_bounded_part(k, {x1, y1, z1}, {x2, y2, z2});
        };
        const auto complex_result = DefiniteIntegrals::integrate<GL::Quadrature<3, 3, 3, 3, 3, 3>>(
            integrand, {-l / 2, -l / 2, -l / 2, -l / 2, -l / 2, -l / 2 + distance},
            {l, l, l, l, l, l,});

        std::cout.precision(16);
        std::cout << "Dist = " << distance / l << " l-s; Rel error = " << std::abs(complex_result - far_zone_result) / std::abs(complex_result) << std::endl;
        std::cout << "Absolute value complex result: " << complex_result << std::endl;
        std::cout << "Absolute value farzone result: " << far_zone_result << std::endl;
    }
}

TEST_F(MORE_THAN_3D_INTEGRATION, TEST_GRAVITATIONAL_ENERGY_WITH_VS_WITHOUT_EXTRACTION) {
    const Types::scalar l = 1;
    Containers::vector<Types::scalar> distances; distances.resize(200);
    std::iota(distances.begin(), distances.end(), 0);
    std::transform(distances.begin(), distances.end(), distances.begin(), [l](Types::scalar value) { return (value * 0.1) * l; });
    for (auto &&distance : distances) {

        // Счёт интегралов по разным формулам
        const Types::complex_d k = {2 * 3.1415926 / 10, 0};
        // 1. Расчет с выделением особенности

        // 1.1 Собственно само выделение особенности через интегрирование Ньютонова потенциала куба
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
        // Расчет Ньютонова потенциала куба (численный через аналитику на сечениях)
        const auto newtonian_potential_of_cube = [l, cut_of_cube_integral](scalar x, scalar y, scalar z) {
            const auto integrand = [x, y, z, cut_of_cube_integral](Types::scalar z_dash) {
                return cut_of_cube_integral(x, y, z, z_dash);
            };
            return DefiniteIntegrals::integrate<GL::Quadrature<3>>(integrand, {-l/2}, {l/2}) +
                DefiniteIntegrals::integrate<GL::Quadrature<3>>(integrand, {0}, {l/2});;
        };

        const auto integrand_with = [k](Types::scalar x1, Types::scalar y1, Types::scalar z1, Types::scalar x2, Types::scalar y2,
                                   Types::scalar z2) {
            return Helmholtz::F_bounded_part(k, {x1, y1, z1}, {x2, y2, z2});
        };

        const auto res_with = DefiniteIntegrals::integrate<GL::Quadrature<2, 2, 2, 2, 2, 2>>(
            integrand_with, {-l / 2, -l / 2, -l / 2, -l / 2, -l / 2, -l / 2 + distance},
            {l, l, l, l, l, l,}) +
                Math::Constants::inverse_4PI<Types::scalar>() * DefiniteIntegrals::integrate<GL::Quadrature<2, 2, 2>>(newtonian_potential_of_cube, {-l / 2, -l / 2, -l / 2 + distance}, {l, l, l});

        // 2. Расчет без выделения особенности
        const auto integrand = [k](Types::scalar x1, Types::scalar y1, Types::scalar z1, Types::scalar x2, Types::scalar y2,
                                   Types::scalar z2) {
            return Helmholtz::F(k, {x1, y1, z1}, {x2, y2, z2});
        };
        const auto res_without = DefiniteIntegrals::integrate<GL::Quadrature<2, 2, 2, 2, 2, 2>>(
            integrand, {-l / 2, -l / 2, -l / 2, -l / 2, -l / 2, -l / 2 + distance},
            {l, l, l, l, l, l,});

        std::cout.precision(16);
        std::cout << std::abs(res_with - res_without) / std::abs(res_with) << "," << distance << std::endl;
    }
}
