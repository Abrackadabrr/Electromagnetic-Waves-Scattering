//
// Created by evgen on 14.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_PRODUCTIONS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_PRODUCTIONS_HPP

#include "types/Types.hpp"

namespace EMW::Math {
    /**
     * Для действительных векторов является скалярным произвдением, для комплексных -- сверткой по одному индексу
     * @return "квазискалярное" произведение
     */
    template<typename t_1, typename t_2>
    decltype(std::declval<t_1>() * std::declval<t_2>()) quasiDot(const Types::Vector3<t_1> &first, const Types::Vector3<t_2> &second) {
        return first.x() * second.x() + first.y() * second.y() + first.z() * second.z();
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_PRODUCTIONS_HPP
