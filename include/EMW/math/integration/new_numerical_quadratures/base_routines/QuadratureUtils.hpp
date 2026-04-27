//
// Created by evgen on 11.02.24.
//

#ifndef NUMEDIE_MATH_INTEGRATORS_NUMERICAL_BASE_ROUTINES_QUADRATUREPOINTS_HPP
#define NUMEDIE_MATH_INTEGRATORS_NUMERICAL_BASE_ROUTINES_QUADRATUREPOINTS_HPP

#include "EMW/types/Types.hpp"

namespace EMW::Math::Integration::Numerical::QuadratureUtils
{
    /*
     * Базовая структура для квадратуры
     */
    template <Types::index dimention_, Types::index order_>
    struct BasePoints
    {
        static constexpr unsigned dim = dimention_;
        static constexpr unsigned order = order_;
        using vector_t = Types::Vector3d;
    };

    /*
     * Бариценртическая координата
     */
    template <typename scalar, size_t dim_>
    struct BarPoint
    {
        using vector_t = Types::Vector3d;
        using bar_point = Containers::array<scalar, dim_ + 1>;
        bar_point point;

        vector_t to_domain(const Containers::array<vector_t, dim_ + 1>& points) const
        {
            vector_t res = vector_t::Zero();
            for (unsigned i = 0; i < dim_ + 1; ++i)
            {
                res += point[i] * points[i];
            }
            return res;
        }

        // TODO: повесить сюда концепт
        template <typename cell_t, typename = std::void_t<std::invoke_result_t<
                      decltype(std::declval<cell_t>().barycentric), bar_point>>>
        vector_t to_domain(cell_t&& cell) const
        {
            return std::forward<cell_t>(cell).barycentric(point);
        }
    };

    template <typename scalar>
    struct BarPoint<scalar, 1>
    {
        scalar w;
        using vector_t = Types::point_t;

        vector_t to_domain(const Containers::array<vector_t, 2>& points) const
        {
            return w * points[0] + (1 - w) * points[1];
        }

        vector_t to_domain(const vector_t& v1, const vector_t& v2) const
        {
            return w * v1 + (1 - w) * v2;
        }

        template <typename cell_t, typename = std::void_t<
                      decltype(std::declval<std::remove_cvref_t<cell_t>>().barycentric(scalar{})),
                      typename std::remove_cvref_t<cell_t>::point_t>>
        typename std::remove_cvref_t<cell_t>::point_t to_domain(cell_t&& cell) const
        {
            return std::forward<cell_t>(cell).barycentric(w);
        }

        operator scalar()
        {
            return w;
        }
    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_QUADRATUREPOINTS_HPP
