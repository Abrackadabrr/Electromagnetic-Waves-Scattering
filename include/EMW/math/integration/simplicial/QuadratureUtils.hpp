//
// Created by evgen on 11.02.24.
//

#ifndef NUMEDIE_MATH_INTEGRATORS_NUMERICAL_BASE_ROUTINES_QUADRATUREPOINTS_HPP
#define NUMEDIE_MATH_INTEGRATORS_NUMERICAL_BASE_ROUTINES_QUADRATUREPOINTS_HPP

#include "EMW/types/Types.hpp"

namespace EMW::Math::Integration::Numerical::Simplicial
{
    namespace QuadratureUtils
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
            static_assert(false, "We have no barycantic coordinates for this dimension");
        };

        template <typename scalar>
        struct BarPoint<scalar, 2>
        {
            scalar u, v;

            template <typename point_t>
            constexpr point_t to_domain(const Containers::array<point_t, 3>& points) const
            {
                return u * points[0] + v * points[1] + (1 - u - v) * points[2];
            }

            template <typename point_t>
            constexpr point_t to_domain(const point_t& v1, const point_t& v2, const point_t& v3) const
            {
                return u * v1 + v * v2 + (1 - u - v) * v3;
            }

            // TODO: add the consept
            template <typename cell_t,
                      typename = std::void_t<
                          decltype(std::declval<std::remove_cvref_t<cell_t>>().barycentric(scalar{}, scalar{})),
                          typename std::remove_cvref_t<cell_t>::point_t>>
            constexpr decltype(auto) to_domain(cell_t&& cell) const
            {
                return std::forward<cell_t>(cell).barycentric(u, v);
            }
        };

        template <typename scalar>
        struct BarPoint<scalar, 1>
        {
            scalar w;

            template <typename point_t>
            constexpr point_t to_domain(const Containers::array<point_t, 2>& points) const
            {
                return w * points[0] + (1 - w) * points[1];
            }

            template <typename point_t>
            constexpr point_t to_domain(const point_t& v1, const point_t& v2) const
            {
                return w * v1 + (1 - w) * v2;
            }

            template <typename cell_t, typename = std::void_t<
                          decltype(std::declval<std::remove_cvref_t<cell_t>>().barycentric(scalar{})),
                          typename std::remove_cvref_t<cell_t>::point_t>>
            constexpr decltype(auto) to_domain(cell_t&& cell) const
            {
                return std::forward<cell_t>(cell).barycentric(w);
            }

            constexpr operator scalar()
            {
                return w;
            }
        };
    }

    template <unsigned int dimension_, unsigned int order_>
    struct GaussianPoints
    {
        static_assert("We have no quadrature of such dimension and order");
    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_QUADRATUREPOINTS_HPP
