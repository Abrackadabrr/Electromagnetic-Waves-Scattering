//
// Created by evgen on 04.10.2024.
//

#ifndef SURFACEFIELDBASE_H
#define SURFACEFIELDBASE_H

#include <ranges>

#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "types/Types.hpp"

namespace EMW::Math {
template <typename field> class SurfaceFieldBase {
  protected:
    using manifold_t = Mesh::SurfaceMesh;
    using field_t = field;

    const manifold_t &manifold_;
    Containers::vector<field_t> field_data_;
    bool initialized;
    std::string name_ = "default_field_name";

  public:
    /**
     * Консткруктор "никакого" поля на поверхности
     * @param manifold_ref
     */
    SurfaceFieldBase(const manifold_t &manifold_ref)
        : manifold_(manifold_ref), field_data_(Containers::vector<field_t>{}),
          initialized(false) {}

    /**
     * Конструктор поля по форматированным данным
     * @param man_ref
     * @param field_data
     */
    SurfaceFieldBase(const manifold_t &man_ref,
                     const Containers::vector<field_t> &field_data)
        : manifold_(man_ref), field_data_(field_data), initialized(true) {}

    /**
     * Конструктор поля по аналитическому заданию (сеточная функция на узалх)
     * @param man_ref
     * @param function
     */
    SurfaceFieldBase(
        const manifold_t &man_ref,
        const std::function<field_t(const Mesh::point_t &)> &function);

    /**
     * Конструктор поля по ячейкам (сеточная яункция по ячейкам)
     * @param man_ref
     * @param function
     */
    SurfaceFieldBase(
        const manifold_t &man_ref,
        const std::function<field_t(const Mesh::IndexedCell &)> &function);

    // --- Getters & Setters --- //
    [[nodiscard]] const manifold_t &getManifold() const { return manifold_; };

    [[nodiscard]] const Containers::vector<field_t> &getField() const {
        return field_data_;
    };

    [[nodiscard]] std::string getName() const { return name_; };

    void setName(const std::string &name) { name_ = name; };
};

template <typename field_t>
SurfaceFieldBase<field_t>::SurfaceFieldBase(
    const SurfaceFieldBase<field_t>::manifold_t &man_ref,
    const std::function<field_t(const Mesh::point_t &)> &function)
    : manifold_(man_ref), initialized(true) {
    field_data_.reserve(man_ref.getCells().size());
    const auto coll_points_data_view =
        man_ref.getCells() |
        std::views::transform(
            [](const Mesh::IndexedCell &cell) -> Mesh::point_t {
                return cell.collPoint_.point_;
            });
    // странно, поскольку не могу применить подряд несколько трансформов
    auto points = Containers::vector<Mesh::point_t>{
        coll_points_data_view.begin(), coll_points_data_view.end()};

    const auto filed_data_view = points | std::views::transform(function);
    field_data_ =
        std::vector<field_t>{filed_data_view.begin(), filed_data_view.end()};
    field_data_.shrink_to_fit();
}

template <typename field_t>
SurfaceFieldBase<field_t>::SurfaceFieldBase(
    const manifold_t &man_ref,
    const std::function<field_t(const Mesh::IndexedCell &)> &function)
    : manifold_(man_ref), initialized(true) {
    field_data_.reserve(man_ref.getCells().size());
 const auto filed_data_view = man_ref.getCells() | std::views::transform(function);
 field_data_ = std::vector<field_t>{filed_data_view.begin(), filed_data_view.end()};
}
}

#endif //SURFACEFIELDBASE_H
