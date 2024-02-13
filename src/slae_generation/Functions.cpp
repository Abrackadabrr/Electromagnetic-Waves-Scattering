//
// Created by evgen on 08.02.24.
//
#include "Functions.hpp"
#include "math/MathConstants.hpp"

namespace EMW::Helmholtz {
    Types::complex_d F(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y) {
        const Types::scalar r = (x - y).norm();
        return Math::Constants::inverse_4PI<Types::scalar> * (std::exp(Math::Constants::i * k * r) / r);
    }

    Types::Vector3c V(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y) {
        const Types::Vector3d r_vec = x - y;
        const Types::scalar r2 = (x - y).squaredNorm();
        const Types::scalar r = std::sqrt(r2);
        return (Math::Constants::inverse_4PI<Types::scalar> * std::exp(Math::Constants::i * k * r) / r2) *
               (1 / r - Math::Constants::i * k) * r_vec;
    }

    Types::scalar smoother(Types::scalar e, const Types::Vector3d &x, const Types::Vector3d &y) {
        const Types::scalar r2 = (x - y).squaredNorm();
        const Types::scalar r = std::sqrt(r2);
        return r < e ? (3 - 2 * (r / e)) * (r2 / (e * e)) : 1;
    }

    Types::Vector3c
    reducedK(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y, const Types::Vector3c &j) {
        const Types::Vector3d rVec = x - y;
        const Types::scalar r2 = rVec.squaredNorm();
        const Types::scalar r = std::sqrt(r2);
        const Types::scalar r3 = r2 * r;
        const auto exponent = Math::Constants::inverse_4PI<Types::scalar> * (std::exp(Math::Constants::i * k * r));
        const Types::complex_d firstPart = (-1 / r3 + (Math::Constants::i * k) / r2 + k * k / r);
        const Types::complex_d secondPart = (3 / r3 - (static_cast<Types::scalar>(3) * Math::Constants::i * k) / r2 -
                                             k * k / r);
        return exponent * (j * firstPart + rVec * (rVec.dot(j) * secondPart));
    }
}
