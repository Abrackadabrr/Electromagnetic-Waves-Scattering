//
// Created by evgen on 08.02.24.
//

#include "Functions.hpp"
#include "math/MathConstants.hpp"
#include "math/Productions.hpp"

namespace EMW::Helmholtz {
Types::complex_d F(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y) {

    // добавляю функцию cглаживания
    // тут хардкод, потому что я знаю что я сейчас считаю рупор на грубой сетке
    // чтобы выключить сглаживание надо е поставить нулем.
    const Types::scalar e = 0; // это 2 * h
    const auto smooth_factor = smoother(e, x, y);

    const Types::Vector3d rVec = x - y;
    const Types::scalar r = std::sqrt(Math::quasiDot(rVec, rVec));
    return smooth_factor * Math::Constants::inverse_4PI<Types::scalar>() * (std::exp(Math::Constants::i * k * r) / r);
}

Types::complex_d F_bounded_part(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y) {
    const Types::Vector3d rVec = x - y;
    const Types::scalar r = rVec.norm();
    return Math::Constants::inverse_4PI<Types::scalar>() * ((std::exp(Math::Constants::i * k * r) - 1.) / r);
}

Types::Vector3c V(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y) {
    const Types::Vector3d r_vec = x - y;
    const Types::scalar r2 = (x - y).squaredNorm();
    const Types::scalar r = std::sqrt(r2);
    return (Math::Constants::inverse_4PI<Types::scalar>() * std::exp(Math::Constants::i * k * r) / r2) *
           (1 / r - Math::Constants::i * k) * r_vec;
}

Types::scalar smoother(Types::scalar e, const Types::Vector3d &x, const Types::Vector3d &y) {
    const Types::scalar r2 = (x - y).squaredNorm();
    const Types::scalar r = std::sqrt(r2);
    return r < e ? (3 - 2 * (r / e)) * (r2 / (e * e)) : 1;
}

Types::Vector3c sigmaKernel(Types::complex_d k, const Types::Vector3d &tau, const Types::Vector3d &point_on_surface,
                            const Types::Vector3c &j_e, const Types::Vector3c &j_m) {
    const Types::complex_d exponent = std::exp((-1.) * Math::Constants::i * k * tau.dot(point_on_surface));
    const Types::scalar epsilon = 1;
    const Types::Vector3c vec_e = Math::Constants::i * k * (j_e - tau * Math::quasiDot(j_e, tau)) / epsilon;
    const Types::Vector3c vec_m = -Math::Constants::i * k * Math::cross(tau, j_m);
    return (exponent * (vec_e - vec_m));
}

Types::Vector3c reducedK_kernel(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y,
                                const Types::Vector3c &j) {
    const Types::Vector3d rVec = x - y;
    const Types::scalar r2 = rVec.squaredNorm();
    const Types::scalar r = std::sqrt(r2);
    const Types::scalar r3 = r2 * r;
    const Types::complex_d exponent =
        Math::Constants::inverse_4PI<Types::scalar>() * (std::exp(Math::Constants::i * k * r));
    const Types::complex_d firstPart =
        (-(static_cast<Types::scalar>(1) / r3) + ((Math::Constants::i * k) / r2) + (k * k / r));
    const Types::complex_d secondPart = (static_cast<Types::scalar>(3) / r3 -
                                         (static_cast<Types::scalar>(3) * Math::Constants::i * k) / r2 - (k * k / r)) /
                                        r2;
    return exponent * (j * firstPart + rVec * (Math::quasiDot(rVec, j) * secondPart));
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
    return  (Math::Constants::inverse_4PI<Types::scalar>() / (r2 * r)) * r_vec;
    }
} // namespace EMW::Laplace

