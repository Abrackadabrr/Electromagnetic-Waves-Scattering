//
// Created by evgen on 24.06.2025.
//

#ifndef SHIFTEDPERIODICSTRUCTURE_HPP
#define SHIFTEDPERIODICSTRUCTURE_HPP

#include "mesh/SurfaceMesh.hpp"
#include "mesh/Utils.hpp"
#include "types/Types.hpp"

#include <iostream>

namespace EMW::Geometry {
/**
 * Периодическая структура с отсчетами в плоскости Oxy.
 * Позволяет работать с периодическими структурами
 */
template <Types::index N1, Types::index N2> class ShiftedStructure {
  public:
    template <typename T> using structure_matrix_t = Containers::array<Containers::array<T, N2>, N1>;
    using mesh_t = Mesh::SurfaceMesh;
    static constexpr Types::index N1_ = N1;
    static constexpr Types::index N2_ = N2;

  private:
    Types::scalar h1_, h2_;
    structure_matrix_t<Types::Vector3d> origin_matrix;
    Containers::array<Types::Vector3d, 2> directions;
    Containers::vector<mesh_t> basic_meshes;
    Containers::vector<mesh_t> additional_meshes;

  protected:
    void add_submesh(Types::index i, Mesh::IndexedCell::Tag tag);

    structure_matrix_t<Types::Vector3d> calculate_origin_matrix(Types::Vector3d dir1, Types::Vector3d dir2,
                                                                Types::scalar h1, Types::scalar h2) const;

    std::pair<Types::index, Types::index> static constexpr linear_to_double(Types::index index) noexcept {
        return {index / N2, index % N2};
    }
    Types::index static constexpr double_to_linear(Types::index i, Types::index j) noexcept { return i * N2 + j; }

  public:
    ShiftedStructure(Types::Vector3d dir1, Types::Vector3d dir2, Types::scalar h1, Types::scalar h2,
                     const Mesh::SurfaceMesh &base_mesh);

    // Selectors
    const structure_matrix_t<Types::Vector3d> &get_origin_matrix() const noexcept { return origin_matrix; }
    const Types::Vector3d &get_origin(Types::index i, Types::index j) const noexcept { return origin_matrix[i][j]; }
    const structure_matrix_t<mesh_t> &get_meshes() const noexcept { return basic_meshes; }
    [[nodiscard]] const mesh_t &get(Types::index index) const noexcept { return basic_meshes[index]; }
    [[nodiscard]] const mesh_t &get(Types::index i, Types::index j) const noexcept {
        return basic_meshes[double_to_linear(i, j)];
    }
    [[nodiscard]] static constexpr Types::index rows() noexcept { return N1; }
    [[nodiscard]] static constexpr Types::index cols() noexcept { return N2; }
    [[nodiscard]] static constexpr Types::index size() noexcept { return N1 * N2; }

    ShiftedStructure<2 * N1 - 1, 2 * N2 - 1> expand_without_saving_nice_origin() const;
};

template <Types::index N1, Types::index N2>
typename ShiftedStructure<N1, N2>::template structure_matrix_t<Types::Vector3d>
ShiftedStructure<N1, N2>::calculate_origin_matrix(Types::Vector3d dir1, Types::Vector3d dir2,
                                                  Types::scalar h1, Types::scalar h2) const {
    structure_matrix_t<Types::Vector3d> result{};

    // Предполагаем, что точка (0, 0) -- это середина периодической структуры
    // Строчки == изменение координаты по dir1
    // Столбца == изменение координаты по dir2
    const Types::scalar bias_1 = h1 * ((N1 - 1) / 2);  // сколько нужно сдвинуть по первому направлению
    const Types::scalar bias_2 = h2 * ((N2 - 1) / 2);  // сколько нужно сдвинуть по второму направлению
    for (Types::index j = 0; j != N1; ++j) {
        for (Types::index k = 0; k != N2; ++k) {
            const Types::scalar l1 = h1 * static_cast<Types::scalar>(k) - bias_1;
            const Types::scalar l2 = h2 * static_cast<Types::scalar>(j) - bias_2;
            result[j][k] = dir1 * l1 + dir2 * l2;
        }
    }
    return result;
    // Соответственно будет построена периоддическая структура около переданной сетки
    // Если количество строк и/или столбцов четное, то сетка переедет не дальше, чем [h1/2, h2/2]
    // Если количество строк и стоблцов нечетное, то базовая сетка никуда не переедет
}

template <Types::index N1, Types::index N2>
ShiftedStructure<N1, N2>::ShiftedStructure(Types::Vector3d dir1, Types::Vector3d dir2, Types::scalar h1,
                                           Types::scalar h2, const Mesh::SurfaceMesh &base_mesh)
    : h1_(h1), h2_(h2), directions{dir1, dir2} {
    basic_meshes.resize(N1 * N2);
    origin_matrix = calculate_origin_matrix(dir1, dir2, h1, h2);
    for (Types::index i = 0; i != N1 * N2; ++i) {
        const auto double_index = linear_to_double(i);
        basic_meshes[i] =
            std::move(Mesh::Utils::move_by_vector(base_mesh, origin_matrix[double_index.first][double_index.second]));
    }
}

template <Types::index N1, Types::index N2>
void ShiftedStructure<N1, N2>::add_submesh(Types::index i, Mesh::IndexedCell::Tag tag) {
    additional_meshes.emplace_back(basic_meshes[i].getSubmesh(tag));
}

template<Types::index N1, Types::index N2>
ShiftedStructure<2*N1 - 1, 2*N2 - 1> ShiftedStructure<N1,N2>::expand_without_saving_nice_origin() const {
    return {directions[0], directions[1], h1_, h2_, basic_meshes[0]};
}

} // namespace EMW::Geometry

#endif //SHIFTEDPERIODICSTRUCTURE_HPP
