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
    Types::complex_d wave_number;

    Types::Matrix3c matrix_2_coef(Types::index k, Types::index p);
    Types::Matrix3c matrix_3_coef(Types::index k, Types::index p);

    Types::scalar newton_potential_over_cube(Types::index k, Types::point_t r);

  public:
    explicit operator_K_over_cube_mesh(const Mesh::VolumeMesh::CubeMesh &mesh) : mesh(mesh) {};
};

} // namespace EMW::Operators::Volume

#endif // OPERATORK_HPP
