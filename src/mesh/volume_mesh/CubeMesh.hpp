//
// Created by evgen on 17.01.2026.
//

#ifndef CUBEMESH_HPP
#define CUBEMESH_HPP

#include "VolumeCells.hpp"

namespace EMW::Mesh::VolumeMesh {
class CubeMesh {

    struct mesh_info_t {
        Types::index nodes_size;
        Types::index cells_size;
        std::string name;
    };

    using cells_container_t = Containers::vector<VolumeCells::IndexedCube>;
    using nodes_container_t = Containers::vector<Types::point_t>;

    nodes_container_t nodes_;
    cells_container_t cells_;
    Types::index nx_, ny_, nz_;
    Types::scalar dx_, dy_, dz_;

    std::string name = "default_mesh_name";

  public:
    CubeMesh(Types::Vector3d minCorner, Types::scalar xs, Types::scalar ys, Types::scalar zs, std::size_t nx,
             std::size_t ny, std::size_t nz);

    // --- 3D Indecies to 1D --- //
    Types::index inline point_idx(std::size_t i, std::size_t j, std::size_t k) const {
        assert(i < nx_ && j < ny_ && k < nz_);
        return i + nx_ * (j + ny_ * k);
    };

    inline std::size_t cube_idx(std::size_t cx, std::size_t cy, std::size_t cz) const {
        assert(nx_ >= 2 && ny_ >= 2 && nz_ >= 2);
        const std::size_t mx = nx_ - 1;
        const std::size_t my = ny_ - 1;
        const std::size_t mz = nz_ - 1;
        assert(cx < mx && cy < my && cz < mz);

        return cx + mx * (cy + my * cz);
    }

    // --- Getters --- //
    const cells_container_t &getCells() const { return cells_; }
    const nodes_container_t &getNodes() const { return nodes_; }
    const std::string getName() const { return name; }
    const Types::scalar h() const { return std::sqrt(dx_*dx_ + dy_*dy_ + dz_*dz_); }
};
} // namespace EMW::Mesh::VolumeMesh

#endif // CUBEMESH_HPP
