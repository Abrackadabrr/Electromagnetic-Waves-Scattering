//
// Created by evgen on 30.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP

#include "mesh/SurfaceMesh.hpp"
#include "mesh/VolumeMesh.hpp"
#include "math/fields/SurfaceField.hpp"
#include "types/Types.hpp"

namespace VTK {
    void surface_snapshot(const EMW::Mesh::SurfaceMesh &mesh, const std::string &part_to_file);

    void field_snapshot(const EMW::Math::SurfaceField & field, const std::string &part_to_file);

    void volume_snapshot(const EMW::Mesh::VolumeMesh &mesh, const std::string &part_to_file);

    void united_snapshot(const EMW::Mesh::SurfaceMesh &mesh, std::initializer_list<EMW::Math::SurfaceField> fields,
                         const std::string &path_to_file);
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP
