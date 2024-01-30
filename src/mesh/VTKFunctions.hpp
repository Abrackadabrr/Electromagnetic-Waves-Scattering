//
// Created by evgen on 30.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP

#include "mesh/Mesh.hpp"
#include "Types.hpp"

namespace EMW::Mesh {
    void test_snapshot(Types::index snap_number, const Mesh::SurfaceMesh &mesh);
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP
