//
// Created by evgen on 02.07.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_PHYSICALCONDITION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_PHYSICALCONDITION_HPP

#include <utility>

#include "math/MathConstants.hpp"
#include "math/Productions.hpp"
#include "types/Types.hpp"

namespace EMW::Physics {
struct planeWaveCase {
    // polarization
    Types::Vector3d E0;
    // wave number
    Types::complex_d k;
    // unit wave vector
    Types::Vector3d k_vec;

    planeWaveCase(Types::Vector3d polarization, Types::complex_d k_fig, Types::Vector3d k_unit_vec);
    [[nodiscard]] Types::Vector3c value(const Types::Vector3d &point) const;
};

struct HertzElectricDipole {
    // Dipole moment vector
    Types::Vector3c p_;
    // Wave number
    Types::complex_d k_;
    // Location
    Types::point_t r0_;

    HertzElectricDipole(const Types::Vector3c &dipole_moment, const Types::complex_d k_fig, const Types::point_t &r0);

    // Расчет поля диполя Герца
    Types::Vector3c operator()(const Types::point_t &r, Types::complex_d epsilon = {1.0, 0.}) const;
};

/**
 * Расчитывает волновое число по частоте в ГИГАГерцах
 */
inline constexpr Types::scalar get_k_on_frquency(const Types::scalar frequency) {
    return 2 * Math::Constants::PI<Types::scalar>() * 1e9 * frequency / Math::Constants::c;
}
} // namespace EMW::Physics
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_PHYSICALCONDITION_HPP
