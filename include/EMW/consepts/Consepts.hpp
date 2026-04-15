//
// Created by evgen on 05.02.2025.
//

#ifndef CONSEPTS_HPP
#define CONSEPTS_HPP

#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"

namespace EMW::Concepts {

// hepler consepts
template <typename T, typename... U>
concept IsAnyOf = (std::same_as<T, U> || ...);

/** Концепт, генерализирующий сетку */
template <typename Topology>
concept manifold_like = requires(Topology topology) {
    { topology.getCells() } -> std::same_as<Containers::vector<typename Topology::CellType>>;
};

/** Концепт, генерализирующий множество сеток (полную геометрию) */
template <typename TopologicalStructure>
concept GeomtricalStructure = requires(TopologicalStructure structure, Types::index index) {
    { structure.get(index) } -> IsAnyOf<const Mesh::SurfaceMesh &, Mesh::SurfaceMesh &> ;
    { structure.size() } -> std::same_as<Types::index>;
};

/** Концепт, описывающий, что такое блок в больших структурных матрицах */
template<typename vector_t, typename block_t>
concept GenericMatrixBlock = requires(block_t block, Types::index i, Types::index j, vector_t vector)
{
    // Наличие оператора с круглыми скобочками
    block(i, j);
    std::is_scalar_v<decltype(block(i, j))> || std::is_same_v<decltype(block(i, j)), Types::complex_d>;
    // Наличие методов rows() и cols()
    { block.rows() } -> std::convertible_to<Types::index>;
    { block.cols() } -> std::convertible_to<Types::index>;
    block * vector;
};

} // namespace EMW::Consepts

#endif // CONSEPTS_HPP
