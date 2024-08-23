//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_NK_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_NK_HPP

#include "types/Types.hpp"
#include "math/integration/QuadraturePoints.hpp"

namespace EMW::DefiniteIntegrals::NewtonCotess {

    template<Types::index... N>
    struct Quadrature {
        static constexpr Types::index dim = sizeof...(N);
        static constexpr Types::index size = (... * N);
        static constexpr Containers::array<Node<dim>, size> nodes = cartesian_product<Quadrature<N>...>();
    };

    template<>
    struct Quadrature<1> {
        static constexpr Types::index dim = 1;
        static constexpr Types::index size = 1;
        static constexpr Containers::array<Node<1>, size> nodes = {Node<1>{{0}, 2.}};
    };

    template<>
    struct Quadrature<2> {
        static constexpr Types::index dim = 1;
        static constexpr Types::index size = 2;
        static constexpr Containers::array<Node<1>, size> nodes = {Node<1>{{-0.5}, 1.}, Node<1>{{0.5}, 1.}};
    };

    template<>
    struct Quadrature<4> {
        static constexpr Types::index dim = 1;
        static constexpr Types::index size = 4;
        static constexpr Containers::array<Node<1>, size> nodes = {Node<1>{{-0.75}, 1. / 2}, Node<1>{{-0.25}, 1. / 2},
                                                                   Node<1>{{0.25}, 1. / 2}, Node<1>{{0.75}, 1. / 2}};
    };

    template<>
    struct Quadrature<8> {
        static constexpr Types::index dim = 1;
        static constexpr Types::index size = 8;
        static constexpr Containers::array<Node<1>, size> nodes = {Node<1>{{-7./8}, 1. / 4}, Node<1>{{-5./8}, 1. / 4},
                                                                   Node<1>{{-3./8}, 1. / 4}, Node<1>{{-1./8}, 1. / 4},
                                                                   Node<1>{{1./8}, 1. / 4}, Node<1>{{3./8}, 1. / 4},
                                                                   Node<1>{{5./8}, 1. / 4}, Node<1>{{7./8}, 1. / 4}};
    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_NK_HPP
