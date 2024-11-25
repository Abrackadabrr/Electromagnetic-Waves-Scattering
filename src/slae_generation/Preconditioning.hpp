//
// Created by evgen on 21.10.2024.
//

#ifndef PRECONDITIONING_HPP
#define PRECONDITIONING_HPP

#include "mesh/SurfaceMesh.hpp"
#include "types/Types.hpp"

namespace EMW::Matrix::Preconditioning {
/*
 * Производит расчет предобуславливателя для задачи n x K[u] = m
 * Возвращает матрицу P, умноженную на k^2 / 4.
 * В таком случае предполагается, что P * K ~=~ (k^2 / 4) * I
 */
Types::MatrixXc getPreconditiotner(const Mesh::SurfaceMesh &mesh, const Types::scalar radius,
                                   const EMW::Types::complex_d k);

Types::MatrixXc getInverseBasedPreconditioner(const Mesh::SurfaceMesh &mesh, const Types::MatrixXc &metrix,
                                              const Types::scalar radius,
                                   const EMW::Types::complex_d k);

}

#endif //PRECONDITIONING_HPP
