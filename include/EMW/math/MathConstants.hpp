//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MATHCONSTANTS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MATHCONSTANTS_HPP

#include "types/Types.hpp"

namespace EMW::Math::Constants {
template <typename T> constexpr T PI() {
    return static_cast<T>(3.141592653589793238462643383279502884197169399375105820974944592307816406286208998670);
}

template <typename T> constexpr T deg_to_rad() { return PI<T>() / static_cast<T>(180); }

template <typename T> constexpr T PI_square() { return PI<T>() * PI<T>(); }

template <typename T> constexpr T inverse_4PI() { return static_cast<T>(1) / (static_cast<T>(4) * PI<T>()); }

constexpr Types::complex_d i = Types::complex_d{0, 1};

// Электродинамическое постоянные
constexpr Types::scalar c = 2.99792458 * 1e8;
constexpr Types::scalar one_div_c = (1. / 2.99792458) * 1e-8;
constexpr Types::scalar mu_0 = 1.25663706127 * 1e-6;
constexpr Types::scalar mu_0_c = 2.99792458 * 1.25663706127 * 1e2;
constexpr Types::scalar e_0_c = 1 / mu_0_c;

} // namespace EMW::Math::Constants

#endif // ELECTROMAGNETIC_WAVES_SCATTERING_MATHCONSTANTS_HPP
