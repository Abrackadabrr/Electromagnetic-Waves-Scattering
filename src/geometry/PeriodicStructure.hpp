//
// Created by evgen on 18.01.2025.
//

#ifndef PERIODICSTRUCTURE_HPP
#define PERIODICSTRUCTURE_HPP

#include "mesh/SurfaceMesh.hpp"
#include "mesh/Utils.hpp"
#include "types/Types.hpp"

#include <iostream>

namespace EMW::Geometry {
/**
 * Периодическая структура с отсчетами в плоскости Oxy.
 * Позволяет работать с периодическими структурами
 */
template <Types::index N1, Types::index N2> class PeriodicStructure {
  public:
    template <typename T> using structure_matrix_t = Containers::array<Containers::array<T, N2>, N1>;
    using mesh_t = Mesh::SurfaceMesh;

  private:
    Types::scalar h1_, h2_;
    structure_matrix_t<Types::Vector3d> origin_matrix;
    Containers::vector<mesh_t> basic_meshes;
    Containers::vector<mesh_t> additional_meshes;

  protected:
    void add_submesh(Types::index i, Mesh::IndexedCell::Tag tag);

    structure_matrix_t<Types::Vector3d> calculate_origin_matrix(Types::scalar h1, Types::scalar h2) const;

    std::pair<Types::index, Types::index> static constexpr linear_to_double(Types::index index) noexcept {
        return {index / N2, index % N2};
    }
    Types::index static constexpr double_to_linear(Types::index i, Types::index j) noexcept { return i * N2 + j; }

  public:
    PeriodicStructure(Types::scalar h1, Types::scalar h2, const Mesh::SurfaceMesh &base_mesh);

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

    PeriodicStructure<2*N1 - 1, 2*N2 - 1> expand_without_saving_nice_origin() const;
};

template <Types::index N1, Types::index N2>
typename PeriodicStructure<N1, N2>::template structure_matrix_t<Types::Vector3d>
PeriodicStructure<N1, N2>::calculate_origin_matrix(Types::scalar h1, Types::scalar h2) const {
    structure_matrix_t<Types::Vector3d> result{};

    // Предполагаем, что точка (0, 0) -- это середина периодической структуры
    // Строчки == изменение координаты по x
    // Столбца == изменение координаты по y
    for (Types::index j = 0; j != N1; ++j) {
        for (Types::index k = 0; k != N2; ++k) {
            result[j][k] =
                Types::Vector3d{h1 * (static_cast<Types::scalar>(k) - static_cast<Types::scalar>(N2 - 1) / 2),
                                h2 * (static_cast<Types::scalar>(j) - static_cast<Types::scalar>(N1 - 1) / 2), 0};
        }
    }
    return result;
    // Соответственно будет построена периоддическая структура около переданной сетки
    // Если количество строк и/или столбцов четное, то сетка переедет не дальше, чем [h1/2, h2/2]
    // Если количество строк и стоблцов нечетное, то базовая сетка никуда не переедет
}

template <Types::index N1, Types::index N2>
PeriodicStructure<N1, N2>::PeriodicStructure(Types::scalar h1, Types::scalar h2, const Mesh::SurfaceMesh &base_mesh)
    : h1_(h1), h2_(h2) {
//    std::cout << N1 << ' ' << N2 << std::endl;
    basic_meshes.resize(N1 * N2);
    origin_matrix = calculate_origin_matrix(h1, h2);
    for (Types::index i = 0; i != N1 * N2; ++i) {
        const auto double_index = linear_to_double(i);
//        std::cout << std::boolalpha << (i == double_to_linear(linear_to_double(i).first, linear_to_double(i).second))
//                                                << std::endl;
        // std::cout << origin_matrix[double_index.first][double_index.second] << std::endl;
        basic_meshes[i] =
            std::move(Mesh::Utils::move_by_vector(base_mesh, origin_matrix[double_index.first][double_index.second]));
    }
}

template <Types::index N1, Types::index N2>
void PeriodicStructure<N1, N2>::add_submesh(Types::index i, Mesh::IndexedCell::Tag tag) {
    additional_meshes.emplace_back(basic_meshes[i].getSubmesh(tag));
}

template<Types::index N1, Types::index N2>
PeriodicStructure<2*N1 - 1, 2*N2 - 1> PeriodicStructure<N1,N2>::expand_without_saving_nice_origin() const {
    return {h1_, h2_, basic_meshes[0]};
}
} // namespace EMW::Geometry

#endif // PERIODICSTRUCTURE_HPP
