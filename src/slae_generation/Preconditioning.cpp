//
// Created by evgen on 21.10.2024.
//


#include "Preconditioning.hpp"
#include "slae_generation/MatrixGeneration.hpp"

EMW::Types::VectorXc EMW::Matrix::Preconditioning::detail::getLine(const Mesh::SurfaceMesh &mesh,
                                                                   const Mesh::point_t &coll_point,
                                                                   const Types::scalar radius) {
        const auto &cell_r = mesh.getCells();
        Types::VectorXc result = Types::VectorXc::Zero(2 * cell_r.size());
        for (int i = 0; i < mesh.getCells().size(); i++) {
            if ((cell_r[i].collPoint_.point_ - coll_point).norm() < radius) {
            }
        }
}

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
                result(i, j) = coefs.a11;
                result(i + N, j) = coefs.a21;
                result(i, j + N) = coefs.a12;
                result(i + N, j + N) = coefs.a22;
            }
        }
    }
    return result;
}
