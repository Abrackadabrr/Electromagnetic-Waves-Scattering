//
// Created by evgen on 30.01.24.
//

#include "SurfaceMesh.hpp"

#include <iostream>
#include <ranges>

namespace EMW::Mesh {
SurfaceMesh::SurfaceMesh(Containers::vector<point_t> nodes,
                         Containers::vector<Containers::array<Types::index, 4>> cells)
    : nodes_(nodes) {
    std::transform(cells.begin(), cells.end(), std::back_inserter(cells_),
                   [&nodes](const Containers::array<Types::index, 4> &indexes) -> IndexedCell {
                       return IndexedCell(indexes, nodes);
                   });
#define WAVEGUIDE_CALCULATION 1
   // std::cout << "Waveguide submesh creation:" << WAVEGUIDE_CALCULATION << std::endl;
#if WAVEGUIDE_CALCULATION
    // БОЛЬШУЩИЙ КОСТЫЛЬ ДЛЯ РАСЧЕТА ВОЛНОВОДА
    // этот параметр показывает количество ячеек в плоскости,
    // где есть матгнитный ток и задается импедансное условие
    // эти точки помечаются отдельно, чтобы в дальнейшем их вынуть
    // в своем формате сетки я точно знаю, что эти точки лежат в конце
    int amount_of_cells_in_active_surface = 200;
   // std::cout << "amount_of_cells_in_active_surface " << amount_of_cells_in_active_surface << std::endl;
    for (int i = 1; i <= amount_of_cells_in_active_surface; i++) {
        cells_[cells_.size() - i].tag = IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION;
    }
#endif
};

SurfaceMesh::SurfaceMesh(Containers::vector<point_t> nodes,
                         Containers::vector<Containers::array<Types::index, 4>> cells,
                         Containers::vector<std::string> tags)
    : nodes_(nodes) {
    std::transform(cells.begin(), cells.end(), std::back_inserter(cells_),
                   [&nodes](const Containers::array<Types::index, 4> &indexes) -> IndexedCell {
                       return IndexedCell(indexes, nodes);
                   });
    for (auto &&[index, cell] : std::ranges::enumerate_view(cells_))
        cell.determineTag(tags[index]);

    std::cout << "Size of waveguide cross-section is "
              << this->getSubmeshInfo(IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION).cells_size << std::endl;
};

SurfaceMesh::SurfaceMesh(Containers::vector<point_t> nodes,
                         Containers::vector<Containers::array<Types::index, 4>> cells,
                         const std::function<point_t(const Containers::array<Types::index, 4> &,
                                                     const Containers::vector<point_t> &)> &getPoint)
    : nodes_(nodes) {
    const auto cellsConstructed =
        cells |
        std::views::transform([&nodes, &getPoint](const Containers::array<Types::index, 4> &indexes) -> IndexedCell {
            auto cell = IndexedCell(indexes, nodes, getPoint);
            return cell;
        });
    cells_ = Containers::vector<IndexedCell>{std::ranges::begin(cellsConstructed), std::ranges::end(cellsConstructed)};
};

void SurfaceMesh::customLocalBasis(
    const std::function<std::array<Types::Vector3d, 3>(const Mesh::IndexedCell &)> &func) {
    for (auto &cell : cells_) {
        const auto J = func(cell);
        cell.normal = J[2];
        cell.tau[0] = J[0];
        cell.tau[1] = J[1];
    }
}

void SurfaceMesh::customCollocationPpoints(const std::function<Types::Vector3d(const Mesh::IndexedCell &)> &func) {
    for (auto &cell : cells_) {
        cell.collPoint_ = func(cell);
    }
}

SurfaceMesh SurfaceMesh::getSubmesh(IndexedCell::Tag tag) const {
    const auto predicate = [tag](const IndexedCell &cell) { return cell.tag == tag; };
    std::vector<IndexedCell> subcells;
    Types::index size = std::count_if(cells_.begin(), cells_.end(), predicate);
    subcells.reserve(size);
    std::copy_if(cells_.begin(), cells_.end(), std::back_inserter(subcells), predicate);
    return {nodes_, std::move(subcells)};
}

SurfaceMesh::mesh_info_t SurfaceMesh::getSubmeshInfo(IndexedCell::Tag tag) const {
    const auto predicate = [&tag](const IndexedCell &cell) { return cell.tag == tag; };
    Types::index size = std::count_if(cells_.begin(), cells_.end(), predicate);
    return {nodes_.size(), size, std::string{"submesh_of_"} + name};
}
}

