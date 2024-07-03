//
// Created by evgen on 02.07.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_PHYSICALCONDITION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_PHYSICALCONDITION_HPP

#include "types/Types.hpp"

namespace EMW::Physics {
    struct physicalConditionsCase {
        // polarization
        Types::Vector3d E0;
        // wave number
        Types::complex_d k;
        // wave unit vector
        Types::Vector3d k_vec;
    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_PHYSICALCONDITION_HPP
