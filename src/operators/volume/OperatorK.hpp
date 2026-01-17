//
// Created by evgen on 17.01.2026.
//

#ifndef OPERATORK_HPP
#define OPERATORK_HPP

#include "./Utils.hpp"

#include "mesh/volume_mesh/CubeMesh.hpp"

namespace EMW::Operators::Volume {

class operator_K_over_cube_mesh {
    const Mesh::VolumeMesh::CubeMesh &mesh;
    Types::complex_d k_;

    Types::Matrix3c matrix_2_coef(Types::index k, Types::index p);
    Types::Matrix3c matrix_3_coef(Types::index k, Types::index p);

  public:
    explicit operator_K_over_cube_mesh(const Mesh::VolumeMesh::CubeMesh &mesh) : mesh(mesh) {};
};

} // namespace EMW::Operators::Volume

#endif // OPERATORK_HPP
