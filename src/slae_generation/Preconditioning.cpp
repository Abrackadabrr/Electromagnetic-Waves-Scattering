//
// Created by evgen on 21.10.2024.
//


#include "Preconditioning.hpp"
#include "slae_generation/MatrixGeneration.hpp"

EMW::Types::MatrixXc EMW::Matrix::Preconditioning::getPreconditiotner(const Mesh::SurfaceMesh &mesh,
                                                                      const Types::scalar radius, const EMW::Types::scalar k) {
    const auto & cells = mesh.getCells();
    const long N = static_cast<long>(cells.size());
    Types::MatrixXc result = Types::MatrixXc::Zero(2 * N, 2 * N);

#pragma omp parallel for schedule(dynamic) num_threads(14) collapse(2)
    for (long i = 0; i < N; ++i) {
        for (long j = 0; j < N; ++j) {
            if ((cells[i].collPoint_.point_ - cells[j].collPoint_.point_).norm() < radius) {
                const auto coefs = getMatrixCoefs(i, j, k, cells);
                result(i, j) = coefs.a22;
                result(i + N, j) = -coefs.a12;
                result(i, j + N) = -coefs.a21;
                result(i + N, j + N) = coefs.a11;
            }
        }
    }
    return result;
}
