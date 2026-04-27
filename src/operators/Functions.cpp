//
// Created by evgen on 08.02.24.
//

#include "operators/Functions.hpp"

#include "math/MathConstants.hpp"
#include "math/Productions.hpp"

namespace EMW::Helmholtz {
Types::complex_d F(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y) {

    // добавляю функцию cглаживания
    // тут хардкод, потому что я знаю что я сейчас считаю рупор на грубой сетке
    // чтобы выключить сглаживание надо е поставить нулем.
    // const Types::scalar e = 0; // это 2 * h
    // const auto smooth_factor = smoother(e, x, y);

    const Types::Vector3d rVec = x - y;
    const Types::scalar r = rVec.norm();
    return Math::Constants::inverse_4PI<Types::scalar>() * (std::exp(Math::Constants::i * k * r) / r);
}

// TODO: тут надо разобраться с устойчивостью формулы при r -> 0. Скорее всего в таком случае нужно переходить
// к разложению в ряд Тейлора до ~ 3 порядка (при r где-то 1е-3). Нужно посмотреть на графики устойчивости
Types::complex_d F_bounded_part(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y) {
    const Types::Vector3d rVec = x - y;
    const Types::scalar r = rVec.norm();
    if (r < 1e-6)
        return (Math::Constants::i - k * r / 2.) * k * Math::Constants::inverse_4PI<Types::scalar>();
    return Math::Constants::inverse_4PI<Types::scalar>() * ((std::exp(Math::Constants::i * k * r) - 1.) / r);
}

Types::Vector3c V(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y) {
    const Types::Vector3d r_vec = x - y;
    const Types::scalar r2 = (x - y).squaredNorm();
    const Types::scalar r = std::sqrt(r2);
    const auto ik = Math::Constants::i * k;
    return ((Math::Constants::inverse_4PI<Types::scalar>() * std::exp(ik * r) / r2) *
            (1 / r - ik)) * r_vec;
}

Types::scalar smoother(Types::scalar e, const Types::Vector3d &x, const Types::Vector3d &y) {
    const Types::scalar r2 = (x - y).squaredNorm();
    const Types::scalar r = std::sqrt(r2);
    return r < e ? (3 - 2 * (r / e)) * (r2 / (e * e)) : 1;
}

Types::Vector3c sigmaKernel(Types::complex_d k, const Types::Vector3d &tau, const Types::Vector3d &point_on_surface,
                            const Types::Vector3c &j_e, const Types::Vector3c &j_m, Types::complex_d epsilon) {
    const auto ik = Math::Constants::i * k;
    const Types::complex_d exponent = std::exp(-ik * tau.dot(point_on_surface));
    const Types::Vector3c vec_e = (j_e - tau * Math::quasiDot(j_e, tau)) / epsilon;
    const Types::Vector3c vec_m = Math::cross(tau, j_m);
    return exponent * ik * (vec_e + vec_m);
}

Types::Vector3c far_zone_integral_kernel(Types::complex_d k, const Types::point_t &r, const Types::Vector3c &j) {
    const Types::scalar inv_r_sqr_norm = 1. / r.squaredNorm();
    const Types::scalar inv_r_norm = std::sqrt(inv_r_sqr_norm);
    const Types::scalar inv_r_qube_norm = inv_r_sqr_norm * inv_r_norm;
    const Types::complex_d exponent = std::exp(Math::Constants::i * k * (1. / inv_r_norm));
    const Types::complex_d k_sqr_div_r_norm = k * k * inv_r_norm;
    const Types::complex_d ik_div_r_sqr_norm = Math::Constants::i * k * inv_r_sqr_norm;
    const Types::complex_d first_part = -inv_r_qube_norm + ik_div_r_sqr_norm + k_sqr_div_r_norm;
    const Types::complex_d second_part =
        static_cast<Types::scalar>(3) * (inv_r_qube_norm - ik_div_r_sqr_norm) - k_sqr_div_r_norm;
    const auto j_dot_r = Math::quasiDot(j, r);
    return exponent * (j * first_part + (r * inv_r_sqr_norm) * j_dot_r * second_part );
    // умножение второй части на inv_r_sqr_norm,
    // т.к.  в формулах должен быть единичный вектор направления r :-)
}

} // namespace EMW::Helmholtz

namespace EMW::Laplace {
Types::scalar F(const Types::Vector3d &x, const Types::Vector3d &y) {
    const Types::Vector3d r_vec = x - y;
    const Types::scalar r2 = (r_vec).squaredNorm();
    const Types::scalar r = std::sqrt(r2);
    return Math::Constants::inverse_4PI<Types::scalar>() / r;
}

Types::Vector3d gradF(const Types::Vector3d &x, const Types::Vector3d &y) {
    const Types::Vector3d r_vec = x - y;
    const Types::scalar r2 = (r_vec).squaredNorm();
    const Types::scalar r = std::sqrt(r2);
    return (Math::Constants::inverse_4PI<Types::scalar>() / (r2 * r)) * r_vec;
}
} // namespace EMW::Laplace
