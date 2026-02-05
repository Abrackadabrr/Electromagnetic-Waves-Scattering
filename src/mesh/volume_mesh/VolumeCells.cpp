//
// Created by evgen on 17.01.2026.
//

#include "VolumeCells.hpp"

#include <bitset>

namespace EMW::Mesh::VolumeCells {
IndexedCube::IndexedCube(const Containers::vector<point_t> &full_points, const nodes_t &full_indices)
    : nodes_(full_indices), volume_((full_points[full_indices[0]] - full_points[full_indices[1]]).norm() *
                                    (full_points[full_indices[0]] - full_points[full_indices[2]]).norm() *
                                    (full_points[full_indices[0]] - full_points[full_indices[4]]).norm()) {
    center_ = Types::Vector3d::Zero();
    for (int i = 0; i < full_indices.size(); i++) {
        center_ += full_points[full_indices[i]];
    }
    center_ /= full_indices.size();
};

Mesh::IndexedCell IndexedCube::getFace(Axis ax, Direction dir, const full_points_t &fp) const {
    switch (ax) {
    case Axis::X:
        return getXface(dir, fp);
        break;
    case Axis::Y:
        return getYface(dir, fp);
        break;
    case Axis::Z:
        return getZface(dir, fp);
        break;
    }
}

Types::point_t IndexedCube::getLeftDownCornerOfFace(Axis ax, Direction dir, const full_points_t &fp) const {
    const auto res = getFace(ax, dir, fp);
    return fp[res.points_[0]];
}

Mesh::IndexedCell IndexedCube::getXface(Direction dir, const full_points_t &fp) const {
    Containers::array<Types::index, 4> xface{
        nodes_[dir],
        nodes_[static_cast<Types::index>(4 - dir)],
        nodes_[6 + dir],
        nodes_[static_cast<Types::index>(2 + 3 * dir)],
    };
    return IndexedCell{xface, fp};
};

Mesh::IndexedCell IndexedCube::getYface(Direction dir, const full_points_t &fp) const {
    Containers::array<Types::index, 4> xface{
        nodes_[2 * dir],
        nodes_[1 + 5 * dir],
        nodes_[5 + 2 * dir],
        nodes_[4 - dir],
    };
    return IndexedCell{xface, fp};
};

Mesh::IndexedCell IndexedCube::getZface(Direction dir, const full_points_t &fp) const {
    Containers::array<Types::index, 4> xface{
        nodes_[4 * dir],
        nodes_[2 + 3 * dir],
        nodes_[3 + 4 * dir],
        nodes_[1 + 5 * dir],
    };
    return IndexedCell{xface, fp};
};

} // namespace EMW::Mesh::VolumeCells
