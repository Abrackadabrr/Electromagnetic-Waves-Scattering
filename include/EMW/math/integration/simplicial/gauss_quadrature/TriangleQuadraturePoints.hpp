//
// Created by evgen on 28.04.2026.
//

#ifndef TRIANGLEQUADRATUREPOINTS_HPP
#define TRIANGLEQUADRATUREPOINTS_HPP

#include "EMW/types/Types.hpp"
#include "EMW/math/integration/simplicial/QuadratureUtils.hpp"

namespace EMW::Math::Integration::Numerical::Simplicial
{
    template <>
    struct GaussianPoints<2, 1> : QuadratureUtils::BasePoints<2, 1>
    {
        using Base = BasePoints<2, 1>;
        static constexpr size_t n_points = 1;
        static constexpr Containers::array<Types::scalar, n_points> weights = {1};
        static constexpr Containers::array<QuadratureUtils::BarPoint<Types::scalar, 2>, n_points> nodes = {
            {1. / 3, 1. / 3},
        };
    };

    template <>
    struct GaussianPoints<2, 2> : QuadratureUtils::BasePoints<2, 2>
    {
        using Base = BasePoints<2, 2>;
        static constexpr size_t n_points = 3;
        static constexpr Containers::array<Types::scalar, n_points> weights = {1. / 3, 1. / 3, 1. / 3};
        static constexpr Containers::array<QuadratureUtils::BarPoint<Types::scalar, 2>, n_points> nodes = {
            QuadratureUtils::BarPoint<Types::scalar, 2>{1. / 6, 1. / 6},
            {1. / 6, 2. / 3},
            {2. / 3, 1. / 6},
        };
    };

    template <>
    struct GaussianPoints<2, 5> : QuadratureUtils::BasePoints<2, 5>
    {
        using Base = BasePoints<2, 5>;
        static constexpr size_t n_points = 7;
        static constexpr Types::scalar l1 = 0.47014206410511508977; // (6 + std::sqrt(15)) / 21
        static constexpr Types::scalar w1 = 0.13239415278850618074; // (155 + std::sqrt(15)) / 1200
        static constexpr Types::scalar w2 = 0.12593918054482715260; // (155 - std::sqrt(15)) / 1200
        static constexpr Types::scalar l2 = 0.10128650732345633880; // (6 - std::sqrt(15)) / 21

        static constexpr Containers::array<Types::scalar, n_points> weights = {
            0.225, w1, w1, w1, w2, w2, w2
        };

        static constexpr Containers::array nodes = {
            QuadratureUtils::BarPoint<Types::scalar, 2>{1. / 3, 1. / 3},

            QuadratureUtils::BarPoint<Types::scalar, 2>{l1, l1},
            QuadratureUtils::BarPoint<Types::scalar, 2>{l1, 1 - 2 * l1},
            QuadratureUtils::BarPoint<Types::scalar, 2>{1 - 2 * l1, l1},

            QuadratureUtils::BarPoint<Types::scalar, 2>{l2, l2},
            QuadratureUtils::BarPoint<Types::scalar, 2>{l2, 1 - 2 * l2},
            QuadratureUtils::BarPoint<Types::scalar, 2>{1 - 2 * l2, l2}
        };
    };

    template <>
    struct GaussianPoints<2, 4> : QuadratureUtils::BasePoints<2, 4>
    {
        using Base = QuadratureUtils::BasePoints<2, 4>;
        static constexpr size_t n_points = 6;

        static constexpr Types::scalar l1 = 0.44594849091596489;
        static constexpr Types::scalar w1 = 0.22338158967801147;
        static constexpr Types::scalar l2 = 0.091576213509770854;
        static constexpr Types::scalar w2 = 0.10995174365532188;

        static constexpr Containers::array<Types::scalar, n_points> weights = {
            w1, w1, w1, w2, w2, w2
        };

        static constexpr Containers::array nodes = {
            QuadratureUtils::BarPoint<Types::scalar, 2>{l1, l1},
            QuadratureUtils::BarPoint<Types::scalar, 2>{l1, 1 - 2 * l1},
            QuadratureUtils::BarPoint<Types::scalar, 2>{1 - 2 * l1, l1},

            QuadratureUtils::BarPoint<Types::scalar, 2>{l2, l2},
            QuadratureUtils::BarPoint<Types::scalar, 2>{l2, 1 - 2 * l2},
            QuadratureUtils::BarPoint<Types::scalar, 2>{1 - 2 * l2, l2}
        };
    };
}

#endif //TRIANGLEQUADRATUREPOINTS_HPP
