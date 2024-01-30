//
// Created by evgen on 30.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP

#include "mesh/Mesh.hpp"
#include "Types.hpp"

namespace VTK {
    void test_snapshot(EMW::Types::index snap_number, const EMW::Mesh::SurfaceMesh &mesh, const std::string &part_to_file);
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP
