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
               (1/r - Math::Constants::i * k) * r_vec;
    }
}
