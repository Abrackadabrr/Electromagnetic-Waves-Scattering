//
// Created by evgen on 04.10.2024.
//

#include "SurfaceScalarField.h"


void EMW::Math::SurfaceScalarField::multiply(
    const std::function<Types::scalar(const EMW::Types::Vector3d &)> &function) {
    for (int i = 0; i < static_cast<int>(manifold_.getCells().size()); i++) {
        const auto cell = manifold_.getCells()[i];
        field_data_[i] = field_data_[i] * function(cell.collPoint_.point_);
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
        const field_t new_field = function(cell.collPoint_.point_, field_data_[i]);
        field_data_[i] = new_field;
    }
}
