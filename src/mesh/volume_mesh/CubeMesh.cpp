//
// Created by evgen on 17.01.2026.
//

#include "CubeMesh.hpp"

namespace EMW::Mesh::VolumeMesh {
CubeMesh::CubeMesh(Types::Vector3d minCorner, Types::scalar xs, Types::scalar ys, Types::scalar zs, std::size_t nx,
                   std::size_t ny, std::size_t nz) {
    //assert(xs > 0.0 && ys > 0.0 && zs > 0.0);
    //assert(nx >= 2 && ny >= 2 && nz >= 2);

    nx_ = nx; ny_ = ny; nz_ = nz;

    // 1) Точки (уникальные)
    nodes_.resize(nx * ny * nz);

    const double dx = xs / static_cast<double>(nx - 1);
    const double dy = ys / static_cast<double>(ny - 1);
    const double dz = zs / static_cast<double>(nz - 1);
    dx_ = dx; dy_ = dy; dz_ = dz;

    for (std::size_t k = 0; k < nz; ++k) {
        const double z = minCorner.z() + dz * static_cast<double>(k);
        for (std::size_t j = 0; j < ny; ++j) {
            const double y = minCorner.y() + dy * static_cast<double>(j);
            for (std::size_t i = 0; i < nx; ++i) {
                const double x = minCorner.x() + dx * static_cast<double>(i);
                nodes_[point_idx(i, j, k)] = Types::Vector3d{x, y, z};
            }
        }
    }

    // 2) Ячейки-кубы (гексаэдры)
    // Каждая ячейка (i,j,k) соединяет "соседние" точки:
    // нижний слой z=k:      v000 v100 v010 v110
    // верхний слой z=k+1:   v001 v101 v011 v111
    cells_.reserve((nx - 1) * (ny - 1) * (nz - 1));

    for (std::size_t k = 0; k + 1 < nz; ++k) {
        for (std::size_t j = 0; j + 1 < ny; ++j) {
            for (std::size_t i = 0; i + 1 < nx; ++i) {
                const auto v000 = static_cast<Types::index>(point_idx(i,   j,   k));
                const auto v100 = static_cast<Types::index>(point_idx(i+1, j,   k));
                const auto v010 = static_cast<Types::index>(point_idx(i,   j+1, k));
                const auto v110 = static_cast<Types::index>(point_idx(i+1, j+1, k));

                const auto v001 = static_cast<Types::index>(point_idx(i,   j,   k+1));
                const auto v101 = static_cast<Types::index>(point_idx(i+1, j,   k+1));
                const auto v011 = static_cast<Types::index>(point_idx(i,   j+1, k+1));
                const auto v111 = static_cast<Types::index>(point_idx(i+1, j+1, k+1));

                // Порядок вершин в ячейке:
                // [0..3] нижняя грань, затем [4..7] верхняя грань
                cells_.push_back(VolumeCells::IndexedCube(nodes_, {v000, v100, v010, v110,
                                   v001, v101, v011, v111}));
            }
        }
    }
}



}; // namespace EMW::Mesh::VolumeMesh
