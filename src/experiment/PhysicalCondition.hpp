//
// Created by evgen on 02.07.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_PHYSICALCONDITION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_PHYSICALCONDITION_HPP

#include <utility>

#include "math/Productions.hpp"
#include "types/Types.hpp"
#include "math/MathConstants.hpp"

namespace EMW::Physics {
    struct planeWaveCase {
        // polarization
        Types::Vector3d E0;
        // wave number
        Types::scalar k;
        // wave vector
        Types::Vector3d k_vec;

        planeWaveCase(Types::Vector3d polarization, Types::scalar k_fig, const Types::Vector3d& k_unit_vec) : E0(
                std::move(polarization)), k(k_fig), k_vec(k_fig * k_unit_vec) {};

        [[nodiscard]] Types::Vector3c value(const Types::Vector3d& point) const{
            return E0 * std::exp(-Math::Constants::i * Math::quasiDot(k_vec, point));
        }
    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_PHYSICALCONDITION_HPP
