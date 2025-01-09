//
// Created by evgen on 04.10.2024.
//

#include "SurfaceScalarField.h"


void EMW::Math::SurfaceScalarField::multiply(
    const std::function<Types::scalar(const EMW::Types::Vector3d &)> &function) {
    for (int i = 0; i < static_cast<int>(manifold_.getCells().size()); i++) {
        const auto cell = manifold_.getCells()[i];
        field_data_[i] = field_data_[i] * function(cell.collPoint_);
    }
}

void EMW::Math::SurfaceScalarField::multiply(
    const std::function<Types::scalar(const EMW::Mesh::IndexedCell &)> &function) {
    for (int i = 0; i < static_cast<int>(manifold_.getCells().size()); i++) {
        const auto cell = manifold_.getCells()[i];
        field_data_[i] = field_data_[i] * function(cell);
    }
}
void EMW::Math::SurfaceScalarField::modify(const std::function<field_t(const Mesh::point_t &, const field_t &)> &function) {
    for (int i = 0; i < static_cast<int>(manifold_.getCells().size()); i++) {
        const auto cell = manifold_.getCells()[i];
        const field_t new_field = function(cell.collPoint_, field_data_[i]);
        field_data_[i] = new_field;
    }
}

EMW::Math::SurfaceScalarField EMW::Math::SurfaceScalarField::tagField(const manifold_t& manifold) {
    Containers::vector<Types::complex_d> field_data{};
    const auto& cells = manifold.getCells();
    field_data.resize(manifold.getCells().size());
#pragma omp parallel for num_threads(14)
    for (int i = 0; i != manifold.getCells().size(); i++) {
        field_data[i] = {static_cast<Types::scalar>(cells[i].tag), 0};
    }
    return SurfaceScalarField(manifold, field_data);
}

EMW::Math::SurfaceScalarField EMW::Math::SurfaceScalarField::sequenceNumberField(const manifold_t& manifold) {
    Containers::vector<Types::complex_d> field_data{};
    const auto& cells = manifold.getCells();
    field_data.resize(manifold.getCells().size());
#pragma omp parallel for num_threads(14)
    for (int i = 0; i != manifold.getCells().size(); i++) {
        field_data[i] = {static_cast<Types::scalar>(i), 0};
    }
    return SurfaceScalarField(manifold, field_data);
}
