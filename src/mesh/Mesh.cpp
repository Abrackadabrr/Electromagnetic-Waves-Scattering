//
// Created by evgen on 30.01.24.
//

#include "Mesh.hpp"
#include <ranges>

namespace EMW::Mesh {
    SurfaceMesh::SurfaceMesh(Containers::vector<Point> nodes,
                             Containers::vector<Containers::array<Types::index, 4>> cells) : nodes_(nodes) {
        const auto cellsConstructed = cells | std::views::transform(
        [const &nodes] (const Containers::array<Types::index, 4> &indexes) -> IndexedCell{
                return IndexedCell(indexes, nodes);
        }
        );
        cells_ = Containers::vector{std::ranges::begin(cellsConstructed), std::ranges::end(cellsConstructed)};
    };
}
