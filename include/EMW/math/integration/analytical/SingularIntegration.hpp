//
// Created by evgen on 23.08.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_SINGULARINTEGRATION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_SINGULARINTEGRATION_HPP

#include "mesh/MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::Math::Integration::Analytical {

/**
 * Интегрирование 1/|r - r'| по кубу r' а затем по кубу r
 */
Types::scalar self_newtonian_energy_over_cube(Types::scalar length);

/**
 * Интегрирование 1/|r - r'| по прямоугольному параллелепипеду ([-a, a] x [-b, b] x [-c, c] в координатах r') в
 * зависимости от r, Для такого случая выведена аналитическая формула в статье "Gravitational potential and energy of
 * homogeneous rectangular parallelepiped" Zakir F. Seidov, P.I. Skvirsky, https://arxiv.org/abs/astro-ph/0002496v1
 *
 * Считается, что центр параллелепипеда лежит в (0, 0, 0), его стороны равны 2a, 2b, 2c,
 * и point -- это вектор относительно центра параллелепипеда
 */
Types::scalar newtonian_potential_of_parallelepiped(const Types::point_t &point, Types::scalar a, Types::scalar b,
                                                    Types::scalar c);

/**
 * Интегрирование функции 1/|r - r'| по поверхностной ячейке (треугольник или четырёхугольник)
 * @param point - r
 * @param cell - ячейка (r')
 * @return значение интеграла по этой ячейке
 */
Types::scalar integrate_1_div_r(const Mesh::point_t &point, const Mesh::IndexedCell &cell);

/**
 * Базовая c-style формула для расчета интеграла 1/|r-r'| по многоугольнику (аппроксимирующему гладкую поверхность)
 * Контейнер для вершин должен поддерживать operator[int]
 */
template <typename vertex_cont_t>
Types::scalar integrate_1_div_r(const Types::point_t &point, const vertex_cont_t &vertexes,
                                const Types::Vector3d &normal) {
    // определяем основные геометрические характеритики
    const Types::scalar d = std::abs((point - vertexes[0]).dot(normal));
    Types::scalar result = 0;
    auto n_vert = std::size(vertexes);
    for (int i = 0; i != n_vert; i++) {
        const auto &rp = vertexes[(i + 1) % n_vert];
        const auto &rm = vertexes[i];
        const Types::Vector3d rm_m_point = rm - point;
        const Types::Vector3d rp_m_point = rp - point;
        const Types::Vector3d l = (rp - rm).normalized();
        const Types::scalar R_minus = rm_m_point.norm();

        // если особая точка лежит не на прямой, содержащей часть границы,
        // то учитываем вклад от этой части границы
        if (std::abs(std::abs((rm_m_point).dot(l)) - R_minus) >= 1e-12) {
            const Types::Vector3d u = l.cross(normal);
            const Types::scalar l_plus = rp_m_point.dot(l);
            const Types::scalar l_minus = rm_m_point.dot(l);
            const Types::scalar p0 = rp_m_point.dot(u);
            const Types::scalar R_plus = rp_m_point.norm();
            const Types::scalar R0_sq = p0 * p0 + d * d;
            // расчет слагаемого от iй части границы
            result += p0 * std::log((R_plus + l_plus) / (R_minus + l_minus)) -
                      d * (std::atan2(p0 * l_plus, R0_sq + d * R_plus) - std::atan2(p0 * l_minus, R0_sq + d * R_minus));
        }
    }
    return result;
}

/**
 * Расчет интеграла 1/|r-r'| по плоскому многоугольнику
 * Контейнер для вершин должен поддерживать operator[int]
 * @warning При малых углах в многоугольнике может быть плохо рассчитана нормаль.  При плохих сетках считать нормаль
 * нужно через QR разложение.
 */
template <typename vertex_cont_t>
Types::scalar integrate_1_div_r(const Types::point_t &point, const vertex_cont_t &vertexes) {
    return integrate_1_div_r(
        point, vertexes, Types::Vector3d{(vertexes[1] - vertexes[0]).cross(vertexes[2] - vertexes[0]).normalized()});
}

}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_SINGULARINTEGRATION_HPP
