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
template <int dim>
using bar_point_t = Containers::array_d<dim>;


template <Types::index dimention_, Types::index order_> struct BasePoints {
  static constexpr unsigned dim = dimention_;
  static constexpr unsigned order = order_;
};

}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_QUADRATUREPOINTS_HPP
