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
    template <typename T> using structure_matrix_t = Containers::array<Containers::array<T, N2>, N1>;
    using mesh_t = Mesh::SurfaceMesh;

    structure_matrix_t<Types::Vector3d> origin_matrix;
    Containers::vector<mesh_t> basic_meshes;

    structure_matrix_t<Types::Vector3d> calculate_origin_matrix(Types::scalar h1, Types::scalar h2) const;

    std::pair<Types::index, Types::index> static constexpr linear_to_double(Types::index index) noexcept {
        return {index / N2, index % N2};
    }

  public:
    PeriodicStructure(Types::scalar h1, Types::scalar h2, const Mesh::SurfaceMesh &base_mesh);

    // Selectors
    const structure_matrix_t<Types::Vector3d> &get_origin_matrix() const noexcept { return origin_matrix; }
    const structure_matrix_t<mesh_t> &get_meshes() const noexcept { return basic_meshes; }
    [[nodiscard]] const mesh_t &get(Types::index index) const noexcept { return basic_meshes[index]; }
    [[nodiscard]] const mesh_t &get(Types::index i, Types::index j) const noexcept {
        return basic_meshes[i * N2 + j];
    }

    [[nodiscard]] static constexpr Types::index rows() noexcept { return N1; }
    [[nodiscard]] static constexpr Types::index cols() noexcept { return N2; }
    [[nodiscard]] static constexpr Types::index size() noexcept { return N1 * N2; }
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
}

template <Types::index N1, Types::index N2>
PeriodicStructure<N1, N2>::PeriodicStructure(Types::scalar h1, Types::scalar h2, const Mesh::SurfaceMesh &base_mesh)  {
    basic_meshes.resize(N1 * N2);
    origin_matrix = calculate_origin_matrix(h1, h2);
    for (Types::index i = 0; i != N1 * N2; ++i) {
        const auto double_index = linear_to_double(i);
        // std::cout << origin_matrix[double_index.first][double_index.second] << std::endl;
        basic_meshes[i] = std::move(Mesh::Utils::move_by_vector(base_mesh, origin_matrix[double_index.first][double_index.second]));
    }
}

} // namespace EMW::Geometry

#endif // PERIODICSTRUCTURE_HPP
