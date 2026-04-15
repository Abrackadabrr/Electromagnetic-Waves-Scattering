//
// Created by evgen on 11.02.24.
//

#ifndef NUMEDIE_MATH_INTEGRATORS_NUMERICAL_BASE_ROUTINES_QUADRATUREPOINTS_HPP
#define NUMEDIE_MATH_INTEGRATORS_NUMERICAL_BASE_ROUTINES_QUADRATUREPOINTS_HPP

#include "Types.hpp"

namespace numedie::math::integration::base_routines::quadrature_utils {

/*
 * Тип узла интегрирования в барицентрических координатах
 */
// TODO: разделить веса и точка в два отдельных массива для лучшей кэш-локальности
template <int dim> struct Node {
  using point_t = Containers::array_d<dim>;
  // барицентрические координаты (без последней координаты, которая вычисляется
  // по предыдущим)
  point_t point;
  // веса квадратуры в этих координатах
  Types::scalar weight;
};

template <Types::index dimention_, Types::index order_> struct BasePoints {
  static constexpr unsigned dim = dimention_;
  static constexpr unsigned order = order_;
};

}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_QUADRATUREPOINTS_HPP
