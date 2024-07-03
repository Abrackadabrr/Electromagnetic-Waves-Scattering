//
// Created by evgen on 30.01.24.
//

#include "Mesh.hpp"
#include <ranges>

namespace EMW::Mesh {
    SurfaceMesh::SurfaceMesh(Containers::vector<Point> nodes,
                             Containers::vector<Containers::array<Types::index, 4>> cells) : nodes_(nodes),
                                                                                             jFilled_(false) {
        const auto cellsConstructed = cells | std::views::transform(
                [&nodes](const Containers::array<Types::index, 4> &indexes) -> IndexedCell {
                    auto cell = IndexedCell(indexes, nodes);
                    return cell;
                }
        );
        cells_ = Containers::vector<IndexedCell>{std::ranges::begin(cellsConstructed),
                                                 std::ranges::end(cellsConstructed)};
    };

    void SurfaceMesh::fillJ(const Types::VectorXc &j) {
        const long N = static_cast<long>(cells_.size());
        for (auto [i, cell]: cells_ | std::views::enumerate) {
            cell.collPoint_.J_ = j(i) * cell.tau[0] + j(i + N) * cell.tau[1];
        }
        jFilled_ = true;
    }
#if 0
    SurfaceMesh::SurfaceMesh(Containers::vector<Point> nodes,
                             Containers::vector<Containers::array<Types::index, 4>> cells,
                             Containers::vector<Types::Vector3c> E_field, Containers::vector<Types::Vector3c> H_field)
            : nodes_(nodes), jFilled_(false) {
        const auto cellsConstructed = cells | std::views::enumerate | std::views::transform(
                [&nodes, &E_field, &H_field](const auto &zip) -> IndexedCell {
                    const auto [i, indexes] = zip;
                    auto cell = IndexedCell(indexes, nodes);
                    cell.collPoint_.E_ = E_field[i];
                    cell.collPoint_.H_ = H_field[i];
                    return cell;
                }
        );
        cells_ = Containers::vector<IndexedCell>{std::ranges::begin(cellsConstructed),
                                                 std::ranges::end(cellsConstructed)};
    }
#endif
}
