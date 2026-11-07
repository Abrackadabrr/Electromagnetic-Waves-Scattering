#ifndef SHIFTEDPERIODICSTRUCTURE_HPP
#define SHIFTEDPERIODICSTRUCTURE_HPP

#include "mesh/SurfaceMesh.hpp"
#include "mesh/Utils.hpp"
#include "types/Types.hpp"

#include <stdexcept>
#include <utility>

namespace EMW::Geometry {
/** Periodic meshes centered around the origin in the plane spanned by dir1 and dir2. */
class PeriodicStructure {
  public:
    template <typename T> using structure_matrix_t = Containers::vector<Containers::vector<T>>;
    using mesh_t = Mesh::SurfaceMesh;

  private:
    Types::index n1_, n2_;
    Types::scalar h1_, h2_;
    structure_matrix_t<Types::Vector3d> origin_matrix;
    Containers::array<Types::Vector3d, 2> directions;
    Containers::vector<mesh_t> basic_meshes;
    Containers::vector<mesh_t> additional_meshes;

    structure_matrix_t<Types::Vector3d> calculate_origin_matrix() const {
        structure_matrix_t<Types::Vector3d> result(n1_, Containers::vector<Types::Vector3d>(n2_));
        const Types::scalar bias_1 = (h1_ * (n1_ - 1)) / 2;
        const Types::scalar bias_2 = (h2_ * (n2_ - 1)) / 2;
        for (Types::index j = 0; j < n2_; ++j)
            for (Types::index i = 0; i < n1_; ++i)
                result[i][j] = directions[0] * (h1_ * static_cast<Types::scalar>(i) - bias_1)
                             + directions[1] * (h2_ * static_cast<Types::scalar>(j) - bias_2);
        return result;
    }

  protected:
    void add_submesh(Types::index i, Mesh::IndexedCell::Tag tag) {
        additional_meshes.emplace_back(basic_meshes[i].getSubmesh(tag));
    }

  public:
    PeriodicStructure(Types::index n1, Types::index n2, Types::Vector3d dir1, Types::Vector3d dir2,
                      Types::scalar h1, Types::scalar h2, const Mesh::SurfaceMesh &base_mesh)
        : n1_(n1), n2_(n2), h1_(h1), h2_(h2), directions{dir1, dir2} {
        if (n1_ == 0 || n2_ == 0) throw std::invalid_argument("PeriodicStructure dimensions must be positive");
        origin_matrix = calculate_origin_matrix();
        basic_meshes.reserve(n1_ * n2_);
        for (Types::index i = 0; i < n1_ * n2_; ++i)
            basic_meshes.emplace_back(Mesh::Utils::move_by_vector(base_mesh, get_origin(i)));
    }

    PeriodicStructure(Types::index n1, Types::index n2, Types::scalar h1, Types::scalar h2,
                      const Mesh::SurfaceMesh &base_mesh)
        : PeriodicStructure(n1, n2, Types::Vector3d{1, 0, 0}, Types::Vector3d{0, 1, 0},
                            h1, h2, base_mesh) {}

    [[nodiscard]] std::pair<Types::index, Types::index> linear_to_double(Types::index index) const noexcept {
        return {index / n2_, index % n2_};
    }
    [[nodiscard]] Types::index double_to_linear(Types::index i, Types::index j) const noexcept {
        return i * n2_ + j;
    }
    [[nodiscard]] const structure_matrix_t<Types::Vector3d> &get_origin_matrix() const noexcept {
        return origin_matrix;
    }
    [[nodiscard]] structure_matrix_t<Types::Vector3d> &get_origin_matrix() noexcept { return origin_matrix; }
    [[nodiscard]] const Types::Vector3d &get_origin(Types::index i, Types::index j) const noexcept {
        return origin_matrix[i][j];
    }
    [[nodiscard]] const Types::Vector3d &get_origin(Types::index i) const noexcept {
        const auto [row, col] = linear_to_double(i);
        return origin_matrix[row][col];
    }
    [[nodiscard]] Types::Vector3d &get_origin(Types::index i) noexcept {
        const auto [row, col] = linear_to_double(i);
        return origin_matrix[row][col];
    }
    [[nodiscard]] const Containers::vector<mesh_t> &get_meshes() const noexcept { return basic_meshes; }
    [[nodiscard]] const mesh_t &get(Types::index i) const noexcept { return basic_meshes[i]; }
    [[nodiscard]] mesh_t &get(Types::index i) noexcept { return basic_meshes[i]; }
    [[nodiscard]] const mesh_t &get(Types::index i, Types::index j) const noexcept {
        return basic_meshes[double_to_linear(i, j)];
    }
    [[nodiscard]] Types::index rows() const noexcept { return n1_; }
    [[nodiscard]] Types::index cols() const noexcept { return n2_; }
    [[nodiscard]] Types::index size() const noexcept { return n1_ * n2_; }

    [[nodiscard]] PeriodicStructure expand_without_saving_nice_origin() const {
        return {2 * n1_ - 1, 2 * n2_ - 1, directions[0], directions[1], h1_, h2_, basic_meshes[0]};
    }
};

} // namespace EMW::Geometry

#endif
