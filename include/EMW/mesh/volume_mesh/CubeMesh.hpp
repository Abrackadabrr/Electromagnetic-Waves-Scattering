//
// Created by evgen on 17.01.2026.
//

#ifndef CUBEMESH_HPP
#define CUBEMESH_HPP

#include "VolumeCells.hpp"

namespace EMW::Mesh::VolumeMesh {
class CubeMesh {
  public:
    using CellsType = VolumeCells::IndexedCube;
    using CellsContainer_t = Containers::vector<CellsType>;
    using NodesContainer_t = Containers::vector<Types::point_t>;

  protected:
    struct mesh_info_t {
        Types::index nodes_size;
        Types::index cells_size;
        std::string name;
    };

    Types::scalar dx_, dy_, dz_;
    Types::point_t rel_center_pos_;
    NodesContainer_t nodes_;
    NodesContainer_t left_down_corners_;
    CellsContainer_t cells_;
    Types::index nx_, ny_, nz_;

    std::string name_ = "default_mesh_name";

  public:
    CubeMesh() = default;
    /**
     * Cетка на прямоугольнике из прямоугольников
     *
     * @param minCorner угол области с наименьшими координатами
     * @param xs, ys, zs длины сторон
     * @param nx, ny, nx количества ячеек разбиения по стороне области
     */
    CubeMesh(Types::Vector3d minCorner, Types::scalar xs, Types::scalar ys, Types::scalar zs, std::size_t nx,
             std::size_t ny, std::size_t nz);

    /**
     * Кубическая сетка в пространстве (пока что оператор умеет работать только с ней)
     *
     * @param minCorner угол куба с наименьшими координатами
     * @param xs длина стороны куба
     * @param nx количество ячеек разбиения по стороне куба
     */
    CubeMesh(Types::Vector3d minCorner, Types::scalar xs, std::size_t nx)
        : CubeMesh(minCorner, xs, xs, xs, nx, nx, nx){};

    // --- Конвертация 3D индекс в 1D --- //
    // индекс точки
    [[nodiscard]] Types::index inline point_idx(std::size_t i, std::size_t j, std::size_t k) const {
        // assert(i < nx_ && j < ny_ && k < nz_);
        return i + nx_ * (j + ny_ * k);
    };

    // индекс куба
    [[nodiscard]] inline std::size_t cube_idx(std::size_t cx, std::size_t cy, std::size_t cz) const {
        // assert(nx_ >= 2 && ny_ >= 2 && nz_ >= 2);
        const std::size_t mx = nx_ - 1;
        const std::size_t my = ny_ - 1;
        const std::size_t mz = nz_ - 1;
        // assert(cx < mx && cy < my && cz < mz);

        return cx + mx * (cy + my * cz);
    }

    // --- Getters --- //
    [[nodiscard]] inline const CellsContainer_t &getCells() const noexcept { return cells_; }
    [[nodiscard]] inline const NodesContainer_t &getNodes() const noexcept { return nodes_; }
    [[nodiscard]] inline const NodesContainer_t &getLDC() const noexcept { return left_down_corners_; }
    [[nodiscard]] inline const std::string &getName() const noexcept { return name_; }
    [[nodiscard]] inline Types::scalar h() const noexcept { return std::sqrt(dx_ * dx_ + dy_ * dy_ + dz_ * dz_); }
    [[nodiscard]] inline const Types::point_t &leftDownCorner(Types::index k) const noexcept {
        return left_down_corners_[k];
    }
    [[nodiscard]] inline Types::scalar dx() const noexcept { return dx_; };
    [[nodiscard]] inline Types::scalar dy() const noexcept { return dy_; };
    [[nodiscard]] inline Types::scalar dz() const noexcept { return dz_; };
    [[nodiscard]] inline size_t nx() const noexcept { return nx_; };
    [[nodiscard]] inline size_t ny() const noexcept { return ny_; };
    [[nodiscard]] inline size_t nz() const noexcept { return nz_; };
    [[nodiscard]] inline size_t nCubesX() const noexcept { return nx_ - 1; };
    [[nodiscard]] inline size_t nCubesY() const noexcept { return ny_ - 1; };
    [[nodiscard]] inline size_t nCubesZ() const noexcept { return nz_ - 1; };
    [[nodiscard]] inline size_t nCells() const noexcept { return (nx_ - 1) * (ny_ - 1) * (nz_ - 1); }
    [[nodiscard]] inline Types::scalar distance(size_t k, size_t p) const noexcept {
        return (left_down_corners_[k] - left_down_corners_[p]).norm();
    };

    // --- Setters --- //
    void setName(const std::string &name) noexcept { name_ = name; };

    // --- Specials --- //
    [[nodiscard]] Eigen::PermutationMatrix<Eigen::Dynamic> getPermutation(size_t Nx, size_t Ny,
                                                                          size_t Nz) const noexcept;
    [[nodiscard]] Eigen::PermutationMatrix<Eigen::Dynamic> getPermutationForCubes(
            size_t Nx, size_t Ny, size_t Nz) const noexcept;

    // --- Calculations --- //
    [[nodiscard]] Containers::array<Mesh::IndexedCell, 6> getFacesOfCube(Types::index k) const;
    [[nodiscard]] Containers::array<Mesh::ParallelogramFace, 6> newGetFacesOfCube(Types::index k) const;
    [[nodiscard]] Containers::array<Mesh::ParallelogramFace, 6> verynewGetFacesOfCube(Types::index k) const;
};
} // namespace EMW::Mesh::VolumeMesh

#endif // CUBEMESH_HPP
