//
// Created by evgen on 21.10.2025.
//

#ifndef GENERALQUADRATURE_HPP
#define GENERALQUADRATURE_HPP

#include "EMW/types/Types.hpp"
#include "SimplicialRange.hpp"

namespace EMW::Math::Integration::Numerical::Simplicial {
namespace detail {
/** Суммирование значения функции в точках */
template <typename quadrature, typename integrand_t, typename points_array_t, Types::index... Idx>
constexpr std::invoke_result_t<integrand_t, typename std::remove_cvref_t<points_array_t>::value_type>
sum_impl(integrand_t &&callable, points_array_t &&points, const std::index_sequence<Idx...> & /**/) {
    return ((quadrature::weights[Idx] * callable(points[Idx])) + ...);
}

// TODO: добавить концепт сюда
template <typename quadrature, typename cell_t, Types::index... Idx,
          typename = std::void_t<typename std::remove_cvref_t<cell_t>::point_t>>
constexpr decltype(auto) get_points(cell_t &&cell, const std::index_sequence<Idx...> & /**/) {
    return Containers::array{quadrature::nodes[Idx].to_domain(cell)...};
}

template <typename quadrature, typename... vertex_t, Types::index... Idx>
constexpr decltype(auto) get_points(const std::index_sequence<Idx...> & /**/, vertex_t &&...vertexes) {
    return Containers::array{quadrature::nodes[Idx].to_domain(vertexes...)...};
}
} // namespace detail

/**
 * Интегрирование функции по заданной ячейке
 * TODO: добавить сюда концепт
 */
template <typename quadrature, typename integrand_t, typename cell_t,
          typename = std::void_t<typename std::remove_cvref_t<cell_t>::point_t, decltype(std::declval<cell_t>().mes())>>
constexpr std::invoke_result_t<std::remove_cvref_t<integrand_t>, typename std::remove_cvref_t<cell_t>::point_t>
integrate(integrand_t &&callable, cell_t &&cell) {
    const auto &points = detail::get_points<quadrature>(std::forward<cell_t>(cell),
                                                        std::make_index_sequence<quadrature::nodes.size()>{});
    return cell.mes() * detail::sum_impl<quadrature>(std::forward<integrand_t>(callable), points,
                                                     std::make_index_sequence<quadrature::nodes.size()>{});
};

/**
 * Квадратурная сумма для симплициальной области в пространстве произвольной размерности
 *
 * @param callable интегранд
 * @param vertexes вершины симплекса соотвествующей размерности
 */
template <typename quadrature, typename integrand_t, typename... vertex_t>
constexpr decltype(auto) quadrature_sum(integrand_t &&callable, vertex_t &&...vertexes) {
    static_assert((quadrature::dim + 1) == sizeof...(vertexes),
                  "quadrature_sum could only be used for simplicial cells (segments, triangles, tetrahedrons,...)");

    decltype(auto) points = detail::get_points<quadrature>(std::make_index_sequence<quadrature::nodes.size()>{},
                                                           std::forward<vertex_t>(vertexes)...);

    return detail::sum_impl<quadrature>(std::forward<integrand_t>(callable), points,
                                        std::make_index_sequence<quadrature::nodes.size()>{});
};

template <typename resType> struct adaptive_integration_result {
    resType value;
    size_t level;
};

/**
 * Квадратурная сумма для логически треугольной области
 *
 * @param callable интегранд
 * @param level уровень подразбиения ячейки
 * @param a, b, c вершины треугольника
 * @param vertexes вершины симплекса соотвествующей размерности
 */
template <typename quadrature, typename integrand_t, typename vertex_t>
constexpr decltype(auto) triangle_quadrature_sum_with_decomposition(integrand_t &&callable, size_t fineness,
                                                                    const vertex_t &a, const vertex_t &b,
                                                                    const vertex_t &c) {
    static_assert(quadrature::dim + 1 == 3, "the quadrature is not for triangle");

    std::invoke_result_t<integrand_t, std::remove_cvref_t<vertex_t>> result{};

    for (auto triangle : QuadratureUtils::TriangleRange(a, b, c, fineness)) {
        decltype(auto) points = detail::get_points<quadrature>(std::make_index_sequence<quadrature::nodes.size()>{},
                                                               triangle.a, triangle.b, triangle.c);

        result += detail::sum_impl<quadrature>(std::forward<integrand_t>(callable), points,
                                               std::make_index_sequence<quadrature::nodes.size()>{});
    }
    result = result / static_cast<Types::scalar>(fineness * fineness);

    return result;
}

// пока что работает только для треугольника
// Изменение стратегии подразбиения должно быть в этой функции захардкожено.
template <typename quadrature, typename integrand_t, typename stop_criterion_t, typename vertex_t>
constexpr decltype(auto) adaptive_quadrature_sum(integrand_t &&callable, stop_criterion_t &&stop_crit,
                                                 const vertex_t &a, const vertex_t &b, const vertex_t &c,
                                                 size_t max_level) {
    static_assert(quadrature::dim + 1 == 3, "the quadrature is not for triangle");
    auto result =
        triangle_quadrature_sum_with_decomposition<quadrature>(std::forward<integrand_t>(callable), 1, a, b, c);
    if (max_level == 0)
        return adaptive_integration_result{result, 1};
    decltype(result) curlev_result{};
    for (size_t cur_level = 1; cur_level <= max_level; ++cur_level) {
        // стратегия обновления мелкости подразбиения
        size_t fineness_of_decomposition = size_t{1} << cur_level;
        // итегрирование по подразбиению
        curlev_result = triangle_quadrature_sum_with_decomposition<quadrature>(
            std::forward<integrand_t>(callable), fineness_of_decomposition, a, b, c);
        // критерий остановки адаптивности
        if (std::invoke(std::forward<stop_criterion_t>(stop_crit), result, curlev_result)) {
            return adaptive_integration_result{curlev_result, cur_level};
        }
        // перезарядка на новую итерацию
        result = curlev_result;
    }
    return adaptive_integration_result{result, max_level};
}

}
#endif // GENERALQUADRATURE_HPP
