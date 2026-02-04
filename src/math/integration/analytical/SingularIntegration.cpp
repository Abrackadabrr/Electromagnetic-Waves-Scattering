//
// Created by evgen on 23.08.24.
//

#include "SingularIntegration.hpp"
#include "math/MathConstants.hpp"

#include <cmath>
#include <iostream>

namespace EMW::Math::AnalyticalIntegration {

Types::scalar integrate_1_div_r(const Mesh::point_t &r, const Mesh::IndexedCell &cell) {
    // определяем основные геометрические характеритики
    const Types::Vector3d n = cell.normal;
    const auto vertex = cell.getVertexAsArray();
    const Types::scalar d = std::abs((r - cell.collPoint_).dot(n));
    Types::scalar result = 0;
    for (int i = 0; i != vertex.size(); i++) {
        const auto &rp = vertex[(i + 1) % vertex.size()];
        const auto &rm = vertex[i];
        // расчитываем геометрию, зависящую от сегмента, если этот сегмент не вырожден
        if ((rp - rm).norm() >= 1e-12) {
            const Types::Vector3d l = (rp - rm).normalized();
            // если особая точка лежит не на прямой, содержащей часть границы
            if (std::abs(std::abs((r - rm).dot(l)) - (r - rm).norm()) >= 1e-12) {
                const Types::Vector3d u = l.cross(n);
                const Types::scalar l_plus = (rp - r).dot(l);
                const Types::scalar l_minus = (rm - r).dot(l);
                const Types::scalar p0 = (rp - r).dot(u);
                const Types::Vector3d p0_vec = (rp - r - l_plus * l) / p0;
                const Types::scalar R_plus = (rp - r).norm();
                const Types::scalar R_minus = (rm - r).norm();
                const Types::scalar R0_sq = p0 * p0 + d * d;
                // расчет слагаемого от iй части границы
                result +=
                    p0_vec.dot(u) *
                    (p0 * std::log((R_plus + l_plus) / (R_minus + l_minus)) -
                     d * (std::atan2(p0 * l_plus, R0_sq + d * R_plus) - std::atan2(p0 * l_minus, R0_sq + d * R_minus)));
            }
        }
    }
    // std::cout << result << std::endl;
    return result;
}

Types::scalar self_newtonian_energy_over_cube(Types::scalar length) {
    const Types::scalar sqrt2 = std::sqrt(2.0);
    const Types::scalar sqrt3 = std::sqrt(3.0);
    const auto pi = Math::Constants::PI<Types::scalar>();

    const Types::scalar C = ((2 * sqrt3 - sqrt2 - 1) / 5) + pi / 3 + std::log((sqrt2 - 1) * (2 - sqrt3));

    const Types::scalar L5 = length * length * length * length * length;
    return -2 * L5 * C;
}

// Первообразная, после тройного интегрирования 1/|r| по кубу
static Types::scalar primitive_F(Types::scalar x, Types::scalar y, Types::scalar z) {
    constexpr Types::scalar eps = std::numeric_limits<Types::scalar>::epsilon();
    const Types::scalar R = std::sqrt(x*x + y*y + z*z);
    if (std::abs(R) < eps)
        return 0.0;

    Types::scalar res = 0.0;

    // yz*log(x+R) + zx*log(y+R) + xy*log(z+R)
    res += (y * z) * std::log(x + R);
    res += (z * x) * std::log(y + R);
    res += (x * y) * std::log(z + R);

    // -1/2 * x^2 * atan( (y z) / (x R) )
    res -= 0.5 * x * x * std::atan2(std::copysign(y * z, x * R), std::abs(x * R));
    res -= 0.5 * y * y * std::atan2(std::copysign(z * x, y * R), std::abs(y * R));
    res -= 0.5 * z * z * std::atan2(std::copysign(x * y, z * R), std::abs(z * R));

    return res;
}

// Ньютонов потенциал для параллелепипеда [-a, a] x [-b, b] x [-c, c]
// Эта штука не робит
Types::scalar newtonian_potential_of_parallelepiped(const Types::point_t &point, Types::scalar a, Types::scalar b,
                                                    Types::scalar c) {
    // Координаты точки расчета относительно параллелепипеда
    const auto X = point.x();
    const auto Y = point.y();
    const auto Z = point.z();
    // Вершины параллелепипеда для интегрирования
    const Types::scalar xs[2] = {-a - X, a - X};
    const Types::scalar ys[2] = {-b - Y, b - Y};
    const Types::scalar zs[2] = {-c - Z, c - Z};

    Types::scalar sum = 0.0;

    Types::scalar mul = 1;
    for (double x : xs) {
        mul *= -1;
        for (double y : ys) {
            mul *= -1;
            for (double z : zs) {
                mul *= -1;
                sum += mul * primitive_F(x, y, z);
            }
        }
    }
    return sum;
}

}
