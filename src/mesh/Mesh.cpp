//
// Created by evgen on 30.01.24.
//

#include "Mesh.hpp"
#include <ranges>

namespace EMW::Mesh {
    SurfaceMesh::SurfaceMesh(Containers::vector<Point> nodes,
                             Containers::vector<Containers::array<Types::index, 4>> cells) : nodes_(nodes) {
        const auto cellsConstructed = cells | std::views::transform(
                [&nodes](const Containers::array<Types::index, 4> &indexes) -> IndexedCell {
                    auto cell = IndexedCell(indexes, nodes);
                    cell.collPoint_.E_ = {0, 0, 2};
                    cell.collPoint_.H_ = {0, 0, 1};
                    return cell;
                }
        );
        cells_ = Containers::vector<IndexedCell>{std::ranges::begin(cellsConstructed),
                                                 std::ranges::end(cellsConstructed)};
    };

//    SurfaceMesh::SurfaceMesh(Containers::vector<Point> nodes,
//                             Containers::vector<Containers::array<Types::index, 4>> cells,
//                             Containers::vector<Types::Vector3d> E_field, Containers::vector<Types::Vector3d> H_field)
//            : nodes_(nodes) {
//        const auto cellsConstructed = cells | std::views::enumerate | std::views::transform(
//                [&nodes, &E_field, &H_field](const Types::index i, const Containers::array<Types::index, 4> &indexes) -> IndexedCell {
//                    auto cell = IndexedCell(indexes, nodes);
//                    cell.collPoint_.E_ = E_field[i];
//                    cell.collPoint_.H_ = H_field[i];
//                    return cell;
//                }
//        );
//        cells_ = Containers::vector<IndexedCell>{std::ranges::begin(cellsConstructed),
//                                                 std::ranges::end(cellsConstructed)};
//    };
}
