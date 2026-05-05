//
// Created by evgen on 27.04.2026.
//

#ifndef SEGMENTQUADRATUREPOINTS_HPP
#define SEGMENTQUADRATUREPOINTS_HPP

#include "EMW/types/Types.hpp"
#include "EMW/math/integration/simplicial/QuadratureUtils.hpp"

namespace EMW::Math::Integration::Numerical::Simplicial
{
    template <>
    struct GaussianPoints<1, 1> : QuadratureUtils::BasePoints<1, 1>
    {
        using Base = QuadratureUtils::BasePoints<1, 1>;
        static constexpr size_t n_points = 1;
        static constexpr Containers::array<Types::scalar, n_points> weights = {1};
        static constexpr Containers::array<QuadratureUtils::BarPoint<Types::scalar, Base::dim>, n_points>
        nodes = {1. / 2};
    };

    template <>
    struct GaussianPoints<1, 2> : QuadratureUtils::BasePoints<1, 1>
    {
        using Base = QuadratureUtils::BasePoints<1, 1>;
        static constexpr size_t n_points = 2;

        static constexpr Containers::array<Types::scalar, n_points> weights = {
            static_cast<Types::scalar>(0.50000000000000000000),
            static_cast<Types::scalar>(0.50000000000000000000)
        };

        static constexpr Containers::array<
            QuadratureUtils::BarPoint<Types::scalar, Base::dim>, n_points>
        nodes = {
            static_cast<Types::scalar>(0.21132486540518711775),
            static_cast<Types::scalar>(0.78867513459481288225)
        };
    };

    template <>
    struct GaussianPoints<1, 5> : QuadratureUtils::BasePoints<1, 5>
    {
        using Base = QuadratureUtils::BasePoints<1, 5>;
        static constexpr size_t n_points = 3;

        static constexpr Containers::array<Types::scalar, n_points> weights = {
            static_cast<Types::scalar>(0.277777777777777777778),
            static_cast<Types::scalar>(0.444444444444444444444),
            static_cast<Types::scalar>(0.277777777777777777778)
        };

        static constexpr Containers::array<QuadratureUtils::BarPoint<Types::scalar, Base::dim>, n_points>
        nodes = {
            static_cast<Types::scalar>(0.112701665379258311482),
            static_cast<Types::scalar>(0.500000000000000000000),
            static_cast<Types::scalar>(0.887298334620741688518)
        };
    };

    template <>
    struct GaussianPoints<1, 7> : QuadratureUtils::BasePoints<1, 7>
    {
        using Base = QuadratureUtils::BasePoints<1, 7>;
        static constexpr size_t n_points = 4;

        static constexpr Containers::array<Types::scalar, n_points> weights = {
            static_cast<Types::scalar>(0.173927422568726928687),
            static_cast<Types::scalar>(0.326072577431273071313),
            static_cast<Types::scalar>(0.326072577431273071313),
            static_cast<Types::scalar>(0.173927422568726928687)
        };

        static constexpr Containers::array<QuadratureUtils::BarPoint<Types::scalar, Base::dim>, n_points>
        nodes = {
            static_cast<Types::scalar>(0.0694318442029737123880),
            static_cast<Types::scalar>(0.330009478207571867599),
            static_cast<Types::scalar>(0.669990521792428132401),
            static_cast<Types::scalar>(0.930568155797026287612)
        };
    };

    template <>
    struct GaussianPoints<1, 9> : QuadratureUtils::BasePoints<1, 9>
    {
        using Base = QuadratureUtils::BasePoints<1, 9>;
        static constexpr size_t n_points = 5;

        static constexpr Containers::array<Types::scalar, n_points> weights = {
            static_cast<Types::scalar>(0.118463442528094543757),
            static_cast<Types::scalar>(0.239314335249683234021),
            static_cast<Types::scalar>(0.284444444444444444444),
            static_cast<Types::scalar>(0.239314335249683234021),
            static_cast<Types::scalar>(0.118463442528094543757)
        };

        static constexpr Containers::array<QuadratureUtils::BarPoint<Types::scalar, Base::dim>, n_points>
        nodes = {
            static_cast<Types::scalar>(0.0469100770306680036012),
            static_cast<Types::scalar>(0.230765344947158454482),
            static_cast<Types::scalar>(0.500000000000000000000),
            static_cast<Types::scalar>(0.769234655052841545518),
            static_cast<Types::scalar>(0.953089922969331996399)
        };
    };

    template <>
    struct GaussianPoints<1, 11> : QuadratureUtils::BasePoints<1, 11>
    {
        using Base = QuadratureUtils::BasePoints<1, 11>;
        static constexpr size_t n_points = 6;

        static constexpr Containers::array<Types::scalar, n_points> weights = {
            static_cast<Types::scalar>(0.0856622461895851725201),
            static_cast<Types::scalar>(0.180380786524069303785),
            static_cast<Types::scalar>(0.233956967286345523695),
            static_cast<Types::scalar>(0.233956967286345523695),
            static_cast<Types::scalar>(0.180380786524069303785),
            static_cast<Types::scalar>(0.0856622461895851725201)
        };

        static constexpr Containers::array<QuadratureUtils::BarPoint<Types::scalar, Base::dim>, n_points>
        nodes = {
            static_cast<Types::scalar>(0.0337652428984239860938),
            static_cast<Types::scalar>(0.169395306766867743169),
            static_cast<Types::scalar>(0.380690406958401545685),
            static_cast<Types::scalar>(0.619309593041598454315),
            static_cast<Types::scalar>(0.830604693233132256831),
            static_cast<Types::scalar>(0.966234757101576013906)
        };
    };
}
#endif //SEGMENTQUADRATUREPOINTS_HPP
