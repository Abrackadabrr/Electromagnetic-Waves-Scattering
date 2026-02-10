//
// Created by evgen on 17.01.2026.
//

#ifndef CUBEMESH_HPP
#define CUBEMESH_HPP

#include "VolumeCells.hpp"

namespace EMW::Mesh::VolumeMesh {
class CubeMesh {
protected:
    struct mesh_info_t {
        Types::index nodes_size;
        Types::index cells_size;
        std::string name;
    };
    using CellsType = VolumeCells::IndexedCube;
    using CellsContainer_t = Containers::vector<CellsType>;
    using NodesContainer_t = Containers::vector<Types::point_t>;

    NodesContainer_t nodes_;
    CellsContainer_t cells_;
    Types::index nx_, ny_, nz_;
    Types::scalar dx_, dy_, dz_;

    std::string name = "default_mesh_name";

    /**
     * Параллелепипедная сетка в пространстве
     */
    CubeMesh(Types::Vector3d minCorner, Types::scalar xs, Types::scalar ys, Types::scalar zs, std::size_t nx,
             std::size_t ny, std::size_t nz);

  public:
    /**
     * Кубическая сетка в пространстве (пока что оператор умеет работать только с ней)
     *
     * @param minCorner угол куба с наименьшими координатами
     * @param xs длина стороны куба
     * @param nx количество ячеек разбиения по стороне куба
     */
    CubeMesh(Types::Vector3d minCorner, Types::scalar xs, std::size_t nx) : CubeMesh(minCorner, xs, xs, xs, nx, nx, nx) {};

    // --- Конвертация 3D индекс в 1D --- //
    // индекс точки
    [[nodiscard]] Types::index inline point_idx(std::size_t i, std::size_t j, std::size_t k) const {
        // assert(i < nx_ && j < ny_ && k < nz_);
        return i + nx_ * (j + ny_ * k);
    };

    // индекс куба
    [[nodiscard]] inline std::size_t cube_idx(std::size_t cx, std::size_t cy, std::size_t cz) const {
        //assert(nx_ >= 2 && ny_ >= 2 && nz_ >= 2);
        const std::size_t mx = nx_ - 1;
        const std::size_t my = ny_ - 1;
        const std::size_t mz = nz_ - 1;
        //assert(cx < mx && cy < my && cz < mz);

        return cx + mx * (cy + my * cz);
    }

    // --- Getters --- //
    [[nodiscard]] const CellsContainer_t &getCells() const { return cells_; }
    [[nodiscard]] const NodesContainer_t &getNodes() const { return nodes_; }
    [[nodiscard]] const std::string& getName() const { return name; }
    [[nodiscard]] Types::scalar h() const { return std::sqrt(dx_*dx_ + dy_*dy_ + dz_*dz_); }
    [[nodiscard]] const Types::point_t& leftDownCorner(Types::index k) const { return nodes_[cells_[k].nodes_[0]]; };
    [[nodiscard]] Types::scalar dx() const { return dx_; };
    [[nodiscard]] Types::scalar dy() const { return dy_; };
    [[nodiscard]] Types::scalar dz() const { return dz_; };
};
} // namespace EMW::Mesh::VolumeMesh

#endif // CUBEMESH_HPP
