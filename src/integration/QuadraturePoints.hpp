//
// Created by evgen on 11.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_QUADRATUREPOINTS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_QUADRATUREPOINTS_HPP

#include "types/Types.hpp"

namespace EMW::DefiniteIntegrals {

    template<Types::index dim>
    struct Node {
        Containers::array <Types::scalar, dim> point;
        Types::scalar weight;
    };

    template<typename... Quadratures, Types::index... Is>
    constexpr Node<sizeof...(Quadratures)> getPoint(Containers::array<Types::index, sizeof...(Quadratures)> counter,
                                                    std::index_sequence<Is...>) {
        return {Containers::array < Types::scalar,
                sizeof...(Quadratures) > {Quadratures::nodes[counter[Is]].point[0]...},
                (... * Quadratures::nodes[counter[Is]].weight)};
    }

/**
 * Декартово произведение одномерных разбиений и тензорное произведение весов
 * @tparam Quadratures одномерные разбиения
 * @return
 */
    template<typename... Quadratures>
    constexpr Containers::array<Node<sizeof...(Quadratures)>, (... * Quadratures::size)> cartesian_product() {
        static_assert((... * Quadratures::dim) == 1);

        constexpr Types::index size = (... * Quadratures::size);
        Containers::array<Types::index, sizeof...(Quadratures)> dimentions = {Quadratures::size...};

        Containers::array <Node<dimentions.size()>, size> products{};
        auto counter = Containers::array < Types::index, dimentions.size() > {};  // массив нулей

        for (auto &product: products) {
            product = getPoint<Quadratures...>(counter, std::make_index_sequence<dimentions.size()>());
            ++(counter.front());
            for (Types::index i = 0; i != dimentions.size() - 1; i++) {
                if (counter[i] == dimentions[i]) {
                    counter[i] = 0;
                    counter[i + 1]++;
                } else
                    break;
            }
        }
        return products;
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_QUADRATUREPOINTS_HPP
