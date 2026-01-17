//
// Created by evgen on 17.01.2026.
//

#ifndef UTILS_HPP
#define UTILS_HPP

#include "types/Types.hpp"

#include "math/MathConstants.hpp"

#include "../Functions.hpp"

namespace EMW::Operators::Volume::Utils {

/**
* Интеграл 1/|x - y| по кубу сначала по х, потом по у
*/
inline Types::complex_d self_interracting_gravitational_cube(Types::scalar L) {
    const Types::scalar sqrt2 = std::sqrt(2.0);
    const Types::scalar sqrt3 = std::sqrt(3.0);
    const Types::scalar C =
        std::log((1.0 + sqrt2) * (2.0 + sqrt3))
        - Math::Constants::PI<Types::scalar>() / 3.0
        + (1.0 + sqrt2 - 2.0 * sqrt3) / 5.0;

    // 2 * L^5 * C
    const Types::scalar L2 = L * L;
    return 2.0 * L2 * L2 * L * C;
}

}
#endif //UTILS_HPP
