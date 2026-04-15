//
// Created by evgen on 14.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_PRODUCTIONS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_PRODUCTIONS_HPP

#include "types/Types.hpp"
#include "math/MathConstants.hpp"

namespace EMW::Math {
// Произведения векторов, аналогичные определениям для вещественных вектором

// Скалярное произведение
    template<typename t_1, typename t_2>
    decltype(std::declval<t_1>() * std::declval<t_2>()) quasiDot(const Types::Vector3<t_1> &first, const Types::Vector3<t_2> &second) {
        return first.x() * second.x() + first.y() * second.y() + first.z() * second.z();
    }

// Векторное произведение
    template<typename t_1, typename t_2>
    Types::Vector3<decltype(std::declval<t_1>() * std::declval<t_2>())> cross(const Types::Vector3<t_1> &first, const Types::Vector3<t_2> &second) {
        return { first.y() * second.z() - first.z() * second.y(),
                 first.z() * second.x() - first.x() * second.z(),
                 first.x() * second.y() - first.y() * second.x() };
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_PRODUCTIONS_HPP
