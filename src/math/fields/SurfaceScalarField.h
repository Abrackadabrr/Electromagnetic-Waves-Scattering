//
// Created by evgen on 04.10.2024.
//

#ifndef SURFACESCALARFIELD_H
#define SURFACESCALARFIELD_H

#include "SurfaceFieldBase.h"
#include <iostream>

namespace EMW::Math {
template <typename Scalar> class SurfaceScalarField : public SurfaceFieldBase<Scalar> {
    static_assert(std::is_same_v<Scalar, Types::scalar> || std::is_same_v<Scalar, Types::complex_d>);
    using Base = SurfaceFieldBase<Scalar>;
    using manifold_t = Base::manifold_t;
    using field_t = Base::field_t;

  public:
    SurfaceScalarField(const manifold_t &manifold_ref) : Base(manifold_ref){};
    SurfaceScalarField(const manifold_t &man_ref, const Containers::vector<field_t> &field_data)
        : Base(man_ref, field_data) {}
    SurfaceScalarField(const manifold_t &manifold_ref, const std::function<field_t(const Mesh::point_t &)> &function)
        : Base(manifold_ref, function) {}
    SurfaceScalarField(const manifold_t &manifold_ref,
                       const std::function<field_t(const Mesh::IndexedCell &)> &function)
        : Base(manifold_ref, function) {}

    static constexpr int getNumberOfComponents() { return 1; }

    void multiply(const std::function<Types::scalar(const Types::Vector3d &)> &function);

    void multiply(const std::function<Types::scalar(const Mesh::IndexedCell &)> &function);

    void modify(const std::function<field_t(const Mesh::point_t &, const field_t &)> &function);

    // ------- Операции для преобразования поля в вектрное представление -------- //
    Types::VectorX<Scalar> formVector() const;

    Types::VectorX<Scalar> formVector(Mesh::IndexedCell::Tag tag) const;

    // ------- Factories ------ //

    static SurfaceScalarField tagField(const manifold_t &manifold);

    static SurfaceScalarField sequenceNumberField(const manifold_t &manifold);

    static SurfaceScalarField fromSLAESolution(const manifold_t &manifold, const Types::VectorX<Scalar> &values);
};
} // namespace EMW::Math

template <typename Scalar>
void EMW::Math::SurfaceScalarField<Scalar>::multiply(
    const std::function<Types::scalar(const EMW::Types::Vector3d &)> &function) {
    for (int i = 0; i < static_cast<int>(this->manifold_.getCells().size()); i++) {
        const auto cell = this->manifold_.getCells()[i];
        this->field_data_[i] = this->field_data_[i] * function(cell.collPoint_);
    }
}

template <typename Scalar>
void EMW::Math::SurfaceScalarField<Scalar>::multiply(
    const std::function<Types::scalar(const EMW::Mesh::IndexedCell &)> &function) {
    for (int i = 0; i < static_cast<int>(this->manifold_.getCells().size()); i++) {
        const auto cell = this->manifold_.getCells()[i];
        this->field_data_[i] = this->field_data_[i] * function(cell);
    }
}

template <typename Scalar>
void EMW::Math::SurfaceScalarField<Scalar>::modify(
    const std::function<field_t(const Mesh::point_t &, const field_t &)> &function) {
    for (int i = 0; i < static_cast<int>(this->manifold_.getCells().size()); i++) {
        const auto cell = this->manifold_.getCells()[i];
        const field_t new_field = function(cell.collPoint_, this->field_data_[i]);
        this->field_data_[i] = new_field;
    }
}

template <typename Scalar>
EMW::Math::SurfaceScalarField<Scalar> EMW::Math::SurfaceScalarField<Scalar>::tagField(const manifold_t &manifold) {
    Containers::vector<Scalar> field_data{};
    const auto &cells = manifold.getCells();
    field_data.resize(manifold.getCells().size());
#pragma omp parallel for num_threads(14)
    for (int i = 0; i != manifold.getCells().size(); i++) {
        field_data[i] = static_cast<Scalar>(cells[i].tag);
    }
    return SurfaceScalarField(manifold, field_data);
}

template <typename Scalar>
EMW::Math::SurfaceScalarField<Scalar>
EMW::Math::SurfaceScalarField<Scalar>::sequenceNumberField(const manifold_t &manifold) {
    Containers::vector<Scalar> field_data{};
    const auto &cells = manifold.getCells();
    field_data.resize(manifold.getCells().size());
#pragma omp parallel for num_threads(14)
    for (int i = 0; i != manifold.getCells().size(); i++) {
        field_data[i] = static_cast<Scalar>(i);
    }
    return {manifold, field_data};
}

template <typename Scalar>
EMW::Math::SurfaceScalarField<Scalar>
EMW::Math::SurfaceScalarField<Scalar>::fromSLAESolution(const manifold_t &manifold,
                                                        const Types::VectorX<Scalar> &values) {
    if (manifold.getCells().size() != values.rows())
        throw std::invalid_argument("SurfaceScalarField::fromSLAESolution: size of vector and mesh_cells mismatch");

    Containers::vector<Scalar> field_data{}; field_data.reserve(manifold.getCells().size());

    for (int i = 0; i < static_cast<int>(values.rows()); i++) {
        field_data.push_back(values[i]);
    }

    auto res = SurfaceScalarField{manifold, field_data};
    res.setName("number_field");
    return res;
}

template <typename Scalar>
EMW::Types::VectorX<Scalar> EMW::Math::SurfaceScalarField<Scalar>::formVector() const{
    Types::VectorX<Scalar> res = Types::VectorX<Scalar>::Zero(static_cast<int>(this->field_data_.size()));
    for (int i = 0; i < static_cast<int>(this->field_data_.size()); i++) {
        res[i] = this->field_data_[i];
    }
    return res;
}

template <typename Scalar>
EMW::Types::VectorX<Scalar> EMW::Math::SurfaceScalarField<Scalar>::formVector(Mesh::IndexedCell::Tag tag) const {
    const auto n_cols = this->manifold_.getSubmeshInfo(tag).cells_size;
    Types::VectorX<Scalar> res = Types::VectorX<Scalar>::Zero(n_cols);
    int result_index = 0;
    for (int i = 0; i < static_cast<int>(this->field_data_.size()); i++) {
        if (this->manifold_.getCells()[i].tag == tag) {
            res[result_index] = this->field_data_[i];
            result_index++;
        }
    }
    if (result_index != n_cols)
        throw std::runtime_error("error in form vector by tag");
    return res;
}

#endif // SURFACESCALARFIELD_H
