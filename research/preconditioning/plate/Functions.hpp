//
// Created by evgen on 19.08.24.
//

#ifndef PRECONDITIONING_FUNCTIONS_HPP
#define PRECONDITIONING_FUNCTIONS_HPP

#include "types/Types.hpp"
#include "math/MathConstants.hpp"
#include "math/Productions.hpp"
#include "mesh/Algorithms.hpp"

using namespace EMW;
using namespace EMW::Types;

complex_d phi(const Vector3d &point, const scalar k) {
    return Math::Constants::inverse_4PI<scalar>() * std::exp(Math::Constants::i * k * point.norm()) *
           (Math::Constants::i * k - 1 / point.norm()) / point.squaredNorm();
}

complex_d xi(const Vector3d &point, const scalar k) {
    const scalar r_2 = point.squaredNorm();
    const scalar r = std::sqrt(r_2);
    const scalar r_3 = r_2 * r;
    return Math::Constants::inverse_4PI<scalar>() * std::exp(Math::Constants::i * k * r) *
           (2 / r_2 - k * k + (1 / r) * (1. - 3. * Math::Constants::i * k)) / r_3;
}

Vector3c analyticalInfinitePlate(const Vector3d &e, const Mesh::IndexedCell &cell, const scalar &k) {
    const Vector3d support{0, 0, 0};
    if (Mesh::Algorithm::PointInTriangle(support, cell.getVertex())) {
        return {0, 0, 0};
    } else {
        const Mesh::point_t point = cell.collPoint_.point_;
        return point * xi(point, k) * Math::quasiDot(point, e) - e * xi(point, k) * point.squaredNorm() -
               2. * phi(point, k) * e;
    }
}

#endif //PRECONDITIONING_FUNCTIONS_HPP
