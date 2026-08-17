//
// Created by evgen on 09.06.2026.
//

#include "experiment/PhysicalCondition.hpp"

#include "math/integration/NumericalIntegration.hpp"

#include <operators/Functions.hpp>

#include <stdexcept>

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

LineSource::LineSource(const Types::Vector3d &direction, const Types::point_t &center,
                       Types::complex_d current_amplitude, Types::complex_d k_fig)
    : direction_(direction), center_(center), current_amplitude_(current_amplitude), k_(k_fig) {
    if (direction_.squaredNorm() == 0.0) {
        throw std::invalid_argument("LineSource direction must be non-zero");
    }
}

Types::Vector3c LineSource::operator()(const Types::point_t &r, Types::complex_d epsilon) const {
    const Types::Vector3c differential_dipole_moment = current_amplitude_ * direction_.cast<Types::complex_d>();
    const auto integrand = [this, &r, epsilon, &differential_dipole_moment](Types::scalar t) {
        const Types::point_t dipole_position = center_ + t * direction_;
        const HertzElectricDipole dipole{differential_dipole_moment, k_, dipole_position};
        return dipole(r, epsilon);
    };
    const auto stop_criterion = [](const Types::Vector3c &previous, const Types::Vector3c &current) {
        return (previous - current).norm() < 1.0e-8 * current.norm() + 1.0e-14;
    };

    return DecartIntegration::adaptive_integrate<DecartIntegration::GaussLegendre::Quadrature<8>>(
               integrand, {-0.5}, {1.0}, stop_criterion, 12)
        .first;
}
};
