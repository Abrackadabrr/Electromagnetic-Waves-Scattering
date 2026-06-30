//
// Created by evgen on 09.06.2026.
//

#include "experiment/PhysicalCondition.hpp"

#include <operators/Functions.hpp>

namespace EMW::Physics {

planeWaveCase::planeWaveCase(Types::Vector3d polarization, Types::complex_d k_fig, Types::Vector3d k_unit_vec)
        : E0(std::move(polarization)), k(k_fig), k_vec(std::move(k_unit_vec)){};

[[nodiscard]] Types::Vector3c planeWaveCase::value(const Types::Vector3d &point) const {
    return E0 * std::exp(-Math::Constants::i * k * Math::quasiDot(k_vec, point));
}

/**
 * Monochromatic point electric dipole
 */
HertzElectricDipole::HertzElectricDipole(const Types::Vector3c &dipole_moment,
                      const Types::complex_d k_fig, const Types::point_t &r0)
      : p_(dipole_moment), k_(k_fig), r0_(r0){};

// Расчет поля диполя Герца
Types::Vector3c HertzElectricDipole::operator()(const Types::point_t &r,
                             Types::complex_d epsilon) const {
    return Math::Constants::inverse_4PI<Types::scalar>() * Helmholtz::far_zone_integral_kernel(k_, r - r0_, p_);
}
};
