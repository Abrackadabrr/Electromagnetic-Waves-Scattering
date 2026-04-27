//
// Created by evgen on 17.01.2026.
//

#include "mesh/volume_mesh/CubeMesh.hpp"

namespace EMW::Mesh::VolumeMesh {
CubeMesh::CubeMesh(Types::Vector3d minCorner, Types::scalar xs, Types::scalar ys, Types::scalar zs, std::size_t nx,
                   std::size_t ny, std::size_t nz) {
    //assert(xs > 0.0 && ys > 0.0 && zs > 0.0);
    //assert(nx >= 2 && ny >= 2 && nz >= 2);

    nx_ = nx;
    ny_ = ny;
    nz_ = nz;

    // 1) Точки (уникальные)
    nodes_.resize(nx * ny * nz);

    const double dx = xs / static_cast<double>(nx - 1);
    const double dy = ys / static_cast<double>(ny - 1);
    const double dz = zs / static_cast<double>(nz - 1);
    dx_ = dx;
    dy_ = dy;
    dz_ = dz;

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
                const auto v000 = static_cast<Types::index>(point_idx(i, j, k));
                const auto v100 = static_cast<Types::index>(point_idx(i + 1, j, k));
                const auto v010 = static_cast<Types::index>(point_idx(i, j + 1, k));
                const auto v110 = static_cast<Types::index>(point_idx(i + 1, j + 1, k));

                const auto v001 = static_cast<Types::index>(point_idx(i, j, k + 1));
                const auto v101 = static_cast<Types::index>(point_idx(i + 1, j, k + 1));
                const auto v011 = static_cast<Types::index>(point_idx(i, j + 1, k + 1));
                const auto v111 = static_cast<Types::index>(point_idx(i + 1, j + 1, k + 1));

                // Порядок вершин в ячейке:
                // [0..3] нижняя грань, затем [4..7] верхняя грань
                cells_.push_back(VolumeCells::IndexedCube(nodes_, {v000, v100, v010, v110,
                                                                   v001, v101, v011, v111}));
            }
        }
    }
}

Eigen::PermutationMatrix<Eigen::Dynamic> CubeMesh::getPermutation(size_t Nx, size_t Ny, size_t Nz) const {
    // Создаем матрицу перестановки
    Eigen::PermutationMatrix<Eigen::Dynamic> p_mat(3 * cells_.size());
    decltype(auto) p_mat_cubes = getPermutationForCubes(Nx, Ny, Nz);
    for (size_t idx = 0; idx < p_mat_cubes.size(); ++idx) {
        p_mat.indices()[3 * idx] = 3 * p_mat_cubes.indices()[idx];
        p_mat.indices()[3 * idx + 1] = 3 * p_mat_cubes.indices()[idx] + 1;
        p_mat.indices()[3 * idx + 2] = 3 * p_mat_cubes.indices()[idx] + 2;
    }
    return p_mat.transpose();
}

Eigen::PermutationMatrix<Eigen::Dynamic> CubeMesh::getPermutationForCubes(size_t Nx, size_t Ny, size_t Nz) const {
    // Создаем матрицу перестановки
    Eigen::PermutationMatrix<Eigen::Dynamic> p_mat(cells_.size());

    // Ходим по кубами и последовательно записываем в матрицу правильные индексы
    size_t new_x_size = nCubesX() / Nx;
    size_t new_y_size = nCubesY() / Ny;
    size_t new_z_size = nCubesZ() / Nz;
    size_t values_per_cube = Nx * Ny * Nz;

    for (size_t zdx = 0; zdx < new_z_size; ++zdx) {
        for (size_t ydx = 0; ydx < new_y_size; ++ydx) {
            for (size_t xdx = 0; xdx < new_x_size; ++xdx) {
                // Ищем нижний левый куб для текущего большого куба
                const size_t lx = xdx * Nx;
                size_t ly = ydx * Ny;
                size_t lz = zdx * Nz;
                // Заполняем матрицу перестановки всеми кубами из большого куба подряд
                for (size_t k = 0; k < Nz; ++k) {
                    for (size_t j = 0; j < Ny; ++j) {
                        for (size_t i = 0; i < Nx; ++i) {
                            // Индексы текущего куба в общей сетке
                            const size_t cur_cube_x = i + lx;
                            const size_t cur_cube_y = j + ly;
                            const size_t cur_cube_z = k + lz;
                            // Линейный индекс куба в сетке
                            const size_t full_mesh_idx = cube_idx(cur_cube_x, cur_cube_y, cur_cube_z);
                            // Сколько индексов уже прошло до этого момента
                            // = 3 * кол-во больших кубов + количество маленьких кубов
                            const size_t index_to_begin =
                                values_per_cube * (xdx + new_x_size * (ydx + new_y_size * zdx)) +  // кол-во бол. кубов
                                (i + Nx * (j + Ny * k)); // количество маленьких кубов в текущем кубе
                            // В векторе соотвествующие компоненты лежат в позициях full_mesh_idx.
                            p_mat.indices()[index_to_begin] = full_mesh_idx;
                        }
                    }
                }
            }
        }
    }
    return p_mat;
}


}; // namespace EMW::Mesh::VolumeMesh
