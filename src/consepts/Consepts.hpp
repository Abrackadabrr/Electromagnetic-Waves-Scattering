//
// Created by evgen on 05.02.2025.
//

#ifndef CONSEPTS_HPP
#define CONSEPTS_HPP

#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"

namespace EMW::Concepts {

/** Концепт, генерализирующий сетку */
template <typename Topology>
concept manifold_like = requires(Topology topology) {
    { topology.getCells() } -> std::same_as<Containers::vector<typename Topology::CellType>>;
};
/** Концепт, генерализирующий множество сеток (полную геометрию) */
template <typename TopologicalStructure>
concept GeomtricalStructure = requires(TopologicalStructure structure, Types::index index) {
    { structure.get(index) } -> std::same_as<const Mesh::SurfaceMesh &>;
    { structure.size() } -> std::same_as<Types::index>;
};

} // namespace EMW::Consepts

#endif // CONSEPTS_HPP
