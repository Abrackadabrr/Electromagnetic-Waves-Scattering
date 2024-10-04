//
// Created by evgen on 15.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_ESA_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_ESA_HPP

#include "math/fields/SurfaceVectorField.hpp"
#include "mesh/MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::ESA {
    Types::Vector3c sigmaOverCell(Types::complex_d k, const Types::Vector3d &tau,
                                  const Mesh::IndexedCell &cell, const Types::Vector3c& j);

    Types::scalar calculateESA(const Types::Vector3d &tau, Types::complex_d k,
                               const Math::SurfaceVectorField &j_field);
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_ESA_HPP
