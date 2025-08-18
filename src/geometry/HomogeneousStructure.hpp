//
// Created by evgen on 18.01.2025.
//

#ifndef HOMOGENEOUSSTRUCTURE_HPP
#define HOMOGENEOUSSTRUCTURE_HPP

#include "types/Types.hpp"

#include "mesh/SurfaceMesh.hpp"
#include "mesh/MeshTypes.hpp"

namespace EMW::Geometry {
class HomogeneousStructure {
public:
  using mesh_t = Mesh::SurfaceMesh;
private:

    Containers::vector<Types::Vector3d> origins;
    Containers::vector<mesh_t> basic_meshes;
    Containers::vector<mesh_t> additional_meshes;

protected:
    void add_submesh(Types::index i, Mesh::IndexedCell::Tag tag);

public:
    HomogeneousStructure(const Containers::vector<Types::Vector3d> origins,
                     const Mesh::SurfaceMesh &base_mesh);

    const Containers::vector<Types::Vector3d> &get_origins() const noexcept { return origins; }
    const Types::Vector3d &get_origin(Types::index i) const noexcept { return origins[i]; }
    const Containers::vector<mesh_t> &get_meshes() const noexcept { return basic_meshes; }
    [[nodiscard]] const mesh_t &get(Types::index index) const noexcept { return basic_meshes[index]; }
    [[nodiscard]] constexpr Types::index size() const noexcept { return basic_meshes.size(); }
};

} // namespace EMW::Geometry

#endif // HOMOGENEOUSSTRUCTURE_HPP
