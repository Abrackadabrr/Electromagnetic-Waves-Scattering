//
// Created by evgen on 17.01.2026.
//

#ifndef VOLUMECELLS_HPP
#define VOLUMECELLS_HPP

#include "../MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::Mesh::VolumeCells {

/**
* Этот класс должен быть приватным классом кубической сетки, конечно же, потому что его можно спокойно создать в
* неконсистентном состоянии. Но пока что это просто класс какой-то тут. Нужен для быстро посмотрения расчетов.
*/
struct IndexedCube {
    enum Axis { X, Y, Z };
    enum Direction { Minus, Plus };

    using nodes_t = Containers::array<Types::index, 8>;
    using full_points_t = Containers::vector<point_t>;

    nodes_t nodes_;
    Types::scalar volume_{};
    point_t center_;

    IndexedCube(const full_points_t &full_points, const nodes_t &full_indices);

    Mesh::IndexedCell getFace(Axis ax, Direction dir, const full_points_t& fp) const;

    Mesh::IndexedCell getXface(Direction dir, const full_points_t& fp) const;

    Mesh::IndexedCell getYface(Direction dir, const full_points_t& fp) const;

    Mesh::IndexedCell getZface(Direction dir, const full_points_t& fp) const;
};

} // namespace EMW::Mesh::VolumeCells

#endif // VOLUMECELLS_HPP
