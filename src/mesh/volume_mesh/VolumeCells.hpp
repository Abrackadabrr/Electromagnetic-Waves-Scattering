//
// Created by evgen on 17.01.2026.
//

#ifndef VOLUMECELLS_HPP
#define VOLUMECELLS_HPP

#include "../MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::Mesh::VolumeCells {
/**
 *
 */
struct CubeAlongAxis {
    Types::point_t initial_point;
    Containers::array<Types::scalar, 3> lengths;

    enum Axis { X, Y, Z };
    enum Direction { Minus, Plus };

    Mesh::Cell getXface(Direction dir) const {
        const Types::Vector3d ex = Types::Vector3d{1, 0, 0};
        const Types::Vector3d ey = Types::Vector3d{0, 1, 0};
        const Types::Vector3d ez = Types::Vector3d{0, 0, 1};
        const Types::Vector3d start_point = initial_point + static_cast<Types::index>(dir) * ex;
        const Types::Vector3d point2;
        const Types::Vector3d point3;
        const Types::Vector3d point4;
    }

    Mesh::Cell getYface(Direction dir) const {}

    Mesh::Cell getZface(Direction dir) const {}
};

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
