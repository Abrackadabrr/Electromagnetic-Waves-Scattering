//
// Created by evgen on 21.10.2024.
//

#include "Preconditioning.hpp"
#include "slae_generation/MatrixGeneration.hpp"

#include <operators/OperatorK.hpp>

EMW::Types::MatrixXc EMW::Matrix::Preconditioning::getPreconditiotner(const Mesh::SurfaceMesh &mesh,
                                                                      const Types::scalar radius,
                                                                      const EMW::Types::complex_d k) {
    const auto &cells = mesh.getCells();
    const long N = static_cast<long>(cells.size());
    Types::MatrixXc result = Types::MatrixXc::Zero(2 * N, 2 * N);
    const EMW::Types::complex_d multiplier = 4. / (k * k);
// #pragma omp parallel for schedule(dynamic) num_threads(14) collapse(2)
    for (long i = 0; i < N; ++i) {
        for (long j = 0; j < N; ++j) {
            if ((cells[i].collPoint_ - cells[j].collPoint_).norm() < radius) {
                const auto coefs = Matrix::DiscreteK::getMatrixCoefs(i, j, k, cells);
                result(i, j) = multiplier * coefs.a22;
                result(i + N, j) = -(multiplier * coefs.a12);
                result(i, j + N) = -(multiplier * coefs.a21);
                result(i + N, j + N) = multiplier * coefs.a11;
            }
        }
    }
    return result;
}

EMW::Types::MatrixXc EMW::Matrix::Preconditioning::getInverseBasedPreconditioner(const Mesh::SurfaceMesh &mesh,
                                                                                 const Types::MatrixXc &matrix,
                                                                                 const Types::scalar radius,
                                                                                 const EMW::Types::complex_d k) {
    Types::MatrixXc result = matrix;
    const auto &cells = mesh.getCells();
    const long N = static_cast<long>(cells.size());
#pragma omp parallel for schedule(dynamic) num_threads(14) collapse(2)
    for (long i = 0; i < N; ++i) {
        for (long j = 0; j < N; ++j) {
            if ((cells[i].collPoint_ - cells[j].collPoint_).norm() > radius) {
                result(i, j) = 0;
            }
        }
    }
    return result.inverse();
}
