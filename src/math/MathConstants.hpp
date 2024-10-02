//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MATHCONSTANTS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MATHCONSTANTS_HPP

#include "types/Types.hpp"

namespace EMW::Math::Constants {
    template<typename T>
    constexpr T PI() {
        return static_cast<T>(3.141592653589793238462643383279502884197169399375105820974944592307816406286208998670);
    }
    template<typename T>
    constexpr T inverse_4PI()  {
        return static_cast<T>(1) / (static_cast<T>(4) * PI<T>());
    }

    constexpr Types::complex_d i = Types::complex_d{0, 1};
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_MATHCONSTANTS_HPP
