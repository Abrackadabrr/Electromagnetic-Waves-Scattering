//
// Created by evgen on 18.01.2025.
//

#include "HomogeneousStructure.hpp"

#include "mesh/Utils.hpp"

namespace EMW::Geometry {

HomogeneousStructure::HomogeneousStructure(const Containers::vector<Types::Vector3d> origins,
                                           const Mesh::SurfaceMesh &base_mesh) {
    basic_meshes.resize(origins.size());
    for (Types::index i = 0; i != origins.size(); ++i) {
        basic_meshes[i] =
            std::move(Mesh::Utils::move_by_vector(base_mesh, origins[i]));
    }
}

} // namespace EMW::Geometry
