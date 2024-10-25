//
// Created by evgen on 21.10.2024.
//

#ifndef PRECONDITIONING_HPP
#define PRECONDITIONING_HPP

#include "mesh/SurfaceMesh.hpp"
#include "types/Types.hpp"

namespace EMW::Matrix::Preconditioning {

namespace detail {
Types::VectorXc getLine(const Mesh::SurfaceMesh &mesh, const Mesh::point_t &coll_point, const Types::scalar radius);
} // namespace detail
/*
 * Производит расчет предобуславливателя для задачи n x K[u] = m
 */
Types::MatrixXc getPreconditiotner(const Mesh::SurfaceMesh &mesh, const Types::scalar radius,
                                   const EMW::Types::scalar k);

}

#endif //PRECONDITIONING_HPP
