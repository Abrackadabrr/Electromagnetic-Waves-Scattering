#include "geometry/HomogeneousStructure.hpp"

#include "mesh/Utils.hpp"

namespace EMW::Geometry {

HomogeneousStructure::HomogeneousStructure(const Containers::vector<Types::Vector3d> origins,
                                           const Mesh::SurfaceMesh &base_mesh)
    : origins(origins) {
    basic_meshes.reserve(origins.size());
    for (const auto &origin : origins)
        basic_meshes.emplace_back(Mesh::Utils::move_by_vector(base_mesh, origin));
}

void HomogeneousStructure::add_submesh(Types::index i, Mesh::IndexedCell::Tag tag) {
    additional_meshes.emplace_back(basic_meshes[i].getSubmesh(tag));
}

} // namespace EMW::Geometry
