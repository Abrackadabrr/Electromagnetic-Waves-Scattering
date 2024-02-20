//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_FUNCTIONS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_FUNCTIONS_HPP

#include "types/Types.hpp"

namespace EMW::Helmholtz {
    Types::complex_d F(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y);

    Types::Vector3c V(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y);

    Types::Vector3c
    sigmaKernel(Types::complex_d k, const Types::Vector3d &tau, const Types::Vector3d &y, const Types::Vector3c &j);

    Types::scalar smoother(Types::scalar e, const Types::Vector3d &x, const Types::Vector3d &y);

    Types::Vector3c
    reducedK(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y, const Types::Vector3c &j);
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_FUNCTIONS_HPP
