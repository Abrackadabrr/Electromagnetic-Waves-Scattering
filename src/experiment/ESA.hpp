//
// Created by evgen on 15.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_ESA_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_ESA_HPP

#include "math/fields/SurfaceVectorField.hpp"
#include "mesh/MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::ESA {
Types::Vector3c sigmaOverCell(Types::complex_d k, const Types::Vector3d &tau, const Mesh::IndexedCell &cell,
                              const Types::Vector3c &j);

/**
* Расчет эффективной площади рассеяния для тела, на поверхности которого текут и электрические, и магнитные токи
*/
Types::scalar calculateESA(const Types::Vector3d &tau, Types::complex_d k, const Math::SurfaceVectorField &j_e,
                           const Math::SurfaceVectorField &j_m);

Types::scalar calculateESA(const Types::Vector3d &tau, Types::complex_d k, const Containers::vector<Math::SurfaceVectorField> &j_e_v,
                           const Containers::vector<Math::SurfaceVectorField> &j_m_v);

/**
* Расчет эффективной площади рассеяния для тела, на поверхности которого текут только электрические токи
*/
Types::scalar calculateESA(const Types::Vector3d &tau, Types::complex_d k, const Math::SurfaceVectorField &j_e);
} // namespace EMW::ESA
#endif // ELECTROMAGNETIC_WAVES_SCATTERING_ESA_HPP
