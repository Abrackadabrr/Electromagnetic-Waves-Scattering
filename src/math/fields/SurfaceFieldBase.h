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
    virtual ~SurfaceFieldBase() = default;
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
    field_data_.resize(man_ref.getCells().size());
#pragma omp parallel for num_threads(14) schedule(dynamic)
    for (int i = 0; i < man_ref.getCells().size(); ++i) {
        field_data_[i] = function(man_ref.getCells()[i].collPoint_.point_);
    }
    field_data_.shrink_to_fit();
}

template <typename field_t>
SurfaceFieldBase<field_t>::SurfaceFieldBase(
    const manifold_t &man_ref,
    const std::function<field_t(const Mesh::IndexedCell &)> &function)
    : manifold_(man_ref), initialized(true) {
    field_data_.reserve(man_ref.getCells().size());
#pragma omp parallel for num_threads(14) schedule(dynamic)
    for (int i = 0; i < man_ref.getCells().size(); ++i) {
        field_data_[i] = function(man_ref.getCells()[i]);
    }
    field_data_.shrink_to_fit();
}
}

#endif //SURFACEFIELDBASE_H
