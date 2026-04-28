//
// Created by evgen on 21.10.2025.
//

#ifndef GENERALQUADRATURE_HPP
#define GENERALQUADRATURE_HPP

#include "EMW/types/Types.hpp"

namespace EMW::Math::Integration::Numerical
{
    namespace detail
    {
        /** Суммирование значения функции в точках */
        template <typename quadrature, typename integrand_t, typename points_array_t,
                  Types::index... Idx>
        constexpr std::invoke_result_t<integrand_t, typename std::remove_cvref_t<points_array_t>::value_type> sum_impl(
            integrand_t&& callable, points_array_t&& points,
            const std::index_sequence<Idx...>& /**/)
        {
            return ((quadrature::weights[Idx] * callable(points[Idx])) + ...);
        }

        // TODO: добавить концепт сюда
        template <typename quadrature, typename cell_t, Types::index... Idx,
                  typename = std::void_t<typename std::remove_cvref_t<cell_t>::point_t>>
        constexpr decltype(auto) get_points(cell_t&& cell, const std::index_sequence<Idx...>& /**/)
        {
            return Containers::array{quadrature::nodes[Idx].to_domain(cell)...};
        }

        template <typename quadrature, typename... vertex_t, Types::index... Idx>
        constexpr decltype(auto) get_points(const std::index_sequence<Idx...>& /**/, vertex_t&&... vertexes)
        {
            return Containers::array{quadrature::nodes[Idx].to_domain(vertexes...)...};
        }
    } // namespace detail

    /**
     * Интегрирование функции по заданной ячейке
     * TODO: добавить сюда концепт
     */
    template <typename quadrature, typename integrand_t, typename cell_t, typename = std::void_t<typename
                  std::remove_cvref_t<cell_t>::point_t>>
    constexpr std::invoke_result_t<std::remove_cvref_t<integrand_t>, typename std::remove_cvref_t<cell_t>::point_t>
    integrate(integrand_t&& callable, cell_t&& cell)
    {
        const auto& points = detail::get_points<quadrature>(
            std::forward<cell_t>(cell),
            std::make_index_sequence<quadrature::nodes.size()>{});
        return cell.mes() * detail::sum_impl<quadrature>(
            std::forward<integrand_t>(callable), points,
            std::make_index_sequence<quadrature::nodes.size()>{});
    };

    template <typename quadrature, typename integrand_t, typename... vertex_t>
    constexpr decltype(auto) quadrature_sum(integrand_t&& callable, vertex_t&&... vertexes)
    {
        static_assert(quadrature::dim + 1 == sizeof...(vertexes),
                      "quadrature_sum could only be used for simplicial cells (segments, triangles, tetrahedrons,...)");

        decltype(auto) points = detail::get_points<quadrature>(
            std::make_index_sequence<quadrature::nodes.size()>{},
            std::forward<vertex_t>(vertexes)...);

        return detail::sum_impl<quadrature>(
            std::forward<integrand_t>(callable), points,
            std::make_index_sequence<quadrature::nodes.size()>{});
    };
}
#endif // GENERALQUADRATURE_HPP
