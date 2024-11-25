//
// Created by evgen on 30.01.24.
//

#include "SurfaceMesh.hpp"
#include <ranges>

namespace EMW::Mesh {
    SurfaceMesh::SurfaceMesh(Containers::vector<point_t> nodes,
                             Containers::vector<Containers::array<Types::index, 4>> cells) : nodes_(nodes) {
        const auto cellsConstructed = cells | std::views::transform(
                [&nodes](const Containers::array<Types::index, 4> &indexes) -> IndexedCell {
                    auto cell = IndexedCell(indexes, nodes);
                    return cell;
                }
        );
        cells_ = Containers::vector<IndexedCell>{std::ranges::begin(cellsConstructed),
                                                 std::ranges::end(cellsConstructed)};
    };

    SurfaceMesh::SurfaceMesh(Containers::vector<point_t> nodes, Containers::vector<Containers::array<Types::index, 4>> cells,
        const std::function<point_t(const Containers::array<Types::index, 4> &,
        const Containers::vector<point_t> &)> &getPoint) : nodes_(nodes) {
        const auto cellsConstructed = cells | std::views::transform(
                [&nodes, &getPoint](const Containers::array<Types::index, 4> &indexes) -> IndexedCell {
                    auto cell = IndexedCell(indexes, nodes, getPoint);
                    return cell;
                }
        );
        cells_ = Containers::vector<IndexedCell>{std::ranges::begin(cellsConstructed),
                                                 std::ranges::end(cellsConstructed)};
    };

    void SurfaceMesh::customLocalBasis(const std::function<std::array<Types::Vector3d, 3>(const Mesh::IndexedCell&)>& func) {
        for (auto &cell: cells_) {
            const auto J = func(cell);
            cell.normal = J[2];
            cell.tau[0] = J[0];
            cell.tau[1] = J[1];
        }
    }

    void SurfaceMesh::customCollocationPpoints(const std::function<Types::Vector3d(const Mesh::IndexedCell &)>& func) {
        for (auto &cell: cells_) {
            cell.collPoint_ = func(cell);
        }
    }
}
