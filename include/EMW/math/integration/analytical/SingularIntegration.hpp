//
// Created by evgen on 23.08.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_SINGULARINTEGRATION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_SINGULARINTEGRATION_HPP

#include "mesh/MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::Math::AnalyticalIntegration {
/**
 * Интегрирование функции 1/|r - r'| по поверхностной ячейке (треугольник или четырёхугольник)
 * @param point - r
 * @param cell - ячейка (r')
 * @return значение интеграла по этой ячейке
 */
Types::scalar integrate_1_div_r(const Mesh::point_t &point, const Mesh::IndexedCell &cell);

/**
 * Интегрирование 1/|r - r'| по кубу r' а затем по кубу r
 */
Types::scalar self_newtonian_energy_over_cube(Types::scalar length);

/**
 * Интегрирование 1/|r - r'| по прямоугольному параллелепипеду ([-a, a] x [-b, b] x [-c, c] в координатах r') в зависимости от r,
 * Для такого случая выведена аналитическая формула в статье
 * "Gravitational potential and energy of homogeneous rectangular parallelepiped"
 * Zakir F. Seidov, P.I. Skvirsky, https://arxiv.org/abs/astro-ph/0002496v1
 *
 * Считается, что центр параллелепипеда лежит в (0, 0, 0), его стороны равны 2a, 2b, 2c,
 * и point -- это вектор относительно центра параллелепипеда
*/
Types::scalar newtonian_potential_of_parallelepiped(const Types::point_t &point, Types::scalar a, Types::scalar b, Types::scalar c);

}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_SINGULARINTEGRATION_HPP
