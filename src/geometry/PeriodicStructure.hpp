//
// Created by evgen on 18.01.2025.
//

#ifndef PERIODICSTRUCTURE_HPP
#define PERIODICSTRUCTURE_HPP

#include "mesh/SurfaceMesh.hpp"
#include "mesh/Utils.hpp"
#include "types/Types.hpp"

namespace EMW::Geometry {
/**
 * Периодическая структура с отсчетами в плоскости Oxy.
 * Позволяет работать с периодическими структурами
 */
template <Types::index N1, Types::index N2> class PeriodicStructure {
    template <typename T> using structure_matrix_t = Containers::array<Containers::array<T, N2>, N1>;
    using mesh_t = Mesh::SurfaceMesh;

    structure_matrix_t<Types::Vector3d> origin_matrix;
    structure_matrix_t<mesh_t> meshes_matrix;

    structure_matrix_t<Types::Vector3d> calculate_origin_matrix(Types::scalar h1, Types::scalar h2) const;

    std::pair<Types::index, Types::index> static constexpr linear_to_double(Types::index index) noexcept {
        return {index / N2, index % N2};
    }

  public:
    PeriodicStructure(Types::scalar h1, Types::scalar h2, const Mesh::SurfaceMesh &base_mesh);

    // Selectors
    const structure_matrix_t<Types::Vector3d> &get_origin_matrix() const noexcept { return origin_matrix; }
    const structure_matrix_t<mesh_t> &get_mesh_matrix() const noexcept { return meshes_matrix; }
    const mesh_t &get(Types::index index) const noexcept {
        const auto double_index = linear_to_double(index);
        return meshes_matrix[double_index.first][double_index.second];
    }

    [[nodiscard]] static constexpr Types::index rows() noexcept { return N1; }
    [[nodiscard]] static constexpr Types::index cols() noexcept { return N2; }
    [[nodiscard]] static constexpr Types::index size() noexcept { return N1 * N2; }

    // Подумать, нужны ли сюда методы, которые поддерживают тёплицевость
};

template <Types::index N1, Types::index N2>
typename PeriodicStructure<N1, N2>::template structure_matrix_t<Types::Vector3d>
PeriodicStructure<N1, N2>::calculate_origin_matrix(Types::scalar h1, Types::scalar h2) const {

    structure_matrix_t<Types::Vector3d> result{};

    // Предполагаем, что точка (0, 0) -- это середина периодической структуры
    // Строчки == изменение координаты по x
    // Столбца == изменение координаты по y
    for (Types::index j = 0; j != N1; ++j) {
        for (Types::index k = N2; k != N2; ++k) {
            result[j][k] =
                Types::Vector3d{h1 * (static_cast<Types::scalar>(k) - static_cast<Types::scalar>(N1 / 2)),
                                h2 * (static_cast<Types::scalar>(j) - static_cast<Types::scalar>(N2 / 2)), 0};
        }
    }
    return result;
}

template <Types::index N1, Types::index N2>
PeriodicStructure<N1, N2>::PeriodicStructure(Types::scalar h1, Types::scalar h2, const Mesh::SurfaceMesh &base_mesh)
    : origin_matrix(calculate_origin_matrix(h1, h2)) {

    for (Types::index i = 0; i != N1; ++i) {
        for (Types::index j = 0; j != N2; ++j) {
            meshes_matrix[i][j] = Mesh::Utils::move_by_vector(base_mesh, origin_matrix[i][j]);
        }
    }
}

} // namespace EMW::Geometry

#endif // PERIODICSTRUCTURE_HPP
