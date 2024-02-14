//
// Created by evgen on 14.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_PRODUCTIONS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_PRODUCTIONS_HPP

#include "types/Types.hpp"

namespace EMW {
    template<typename t_1, typename t_2>
    decltype(std::declval<t_1>() * std::declval<t_2>()) dot(const Types::Vector3<t_1> &first, const Types::Vector3<t_2> &second) {
        return first.x() * second.x() + first.y() * second.y() + first.z() * second.z();
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_PRODUCTIONS_HPP
