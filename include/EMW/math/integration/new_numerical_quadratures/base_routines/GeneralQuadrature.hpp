//
// Created by evgen on 21.10.2025.
//

#ifndef GENERALQUADRATURE_HPP
#define GENERALQUADRATURE_HPP

#include "EMW/types/Types.hpp"

namespace EMW::Math::Integration::Numerical {

namespace detail {
/** Суммирование значения функции в точках */
template <typename quadrature, typename integrand_t, typename points_array_t,
          Types::index... Idx>
decltype(auto) sum_impl(integrand_t &&callable, points_array_t &&points,
                        const std::index_sequence<Idx...> & /**/) {
  return ((quadrature::weights[Idx] * callable(points[Idx])) + ...);
}

template <typename quadrature, typename cell_t, Types::index... Idx>
decltype(auto) get_points(cell_t &&cell,
                          const std::index_sequence<Idx...> & /**/) {
  return Containers::array{quadrature::nodes[Idx].to_domain(cell) ... };
}
} // namespace detail

/**
 * Интегрирование функции по заданной ячейке
 */
template <typename quadrature, typename integrand_t, typename cell_t>
decltype(auto) integrate(integrand_t &&callable, cell_t &&cell) {
  // // Базовая проверка на сравнение размерности квадратуры и размерности ячейки
  // static_assert(std::remove_cvref_t<cell_t>::dim ==
  //               quadrature::dim, "Dimension of quadrature must be the same as dimension of cell under integration");

  const auto &points = detail::get_points<quadrature>(
      std::forward<cell_t>(cell),
      std::make_index_sequence<quadrature::nodes.size()>{});
  return cell.mes() * detail::sum_impl<quadrature>(
                          std::forward<integrand_t>(callable), points,
                          std::make_index_sequence<quadrature::nodes.size()>{});
};

}
#endif // GENERALQUADRATURE_HPP
