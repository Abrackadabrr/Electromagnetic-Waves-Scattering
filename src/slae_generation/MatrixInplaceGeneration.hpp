//
// Created by evgen on 04.08.2025.
//

#ifndef MATRIXINPLACEGENERATION_HPP
#define MATRIXINPLACEGENERATION_HPP

#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "types/Types.hpp"
#include "MatrixGeneration.hpp"

namespace EMW::Matrix {

// ---------------- Методы для расчета матриц операторов К и R inplace --------------//
template<typename MatrixLike>
void getMatrixK_inplace(MatrixLike& result, Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh) {
    const auto &cells = surface_mesh.getCells();
    const long N = static_cast<long>(cells.size());
    if (result.rows() != 2 * N || result.cols() != 2 * N)
        throw std::invalid_argument("Matrix size for result is not compatible with the surface mesh in calculation "
                                    "of operator R matrix with one mesh");

#pragma omp parallel for num_threads(14) collapse(2)
    for (long i = 0; i < N; ++i) {
        for (long j = 0; j < N; ++j) {
            const auto coefs = DiscreteK::getMatrixCoefs(i, j, k, cells);
            result(i, j) = coefs.a11;
            result(i + N, j) = coefs.a21;
            result(i, j + N) = coefs.a12;
            result(i + N, j + N) = coefs.a22;
        }
    }
};

template<typename MatrixLike>
void getMatrixK_inplace(MatrixLike& result, Types::complex_d k, const Mesh::SurfaceMesh &integration_mesh,
                           const Mesh::SurfaceMesh &mesh_with_coll_points) {
    const auto &cells = mesh_with_coll_points.getCells();
    const auto &cells_to_integrate = integration_mesh.getCells();
    const long N = static_cast<long>(cells.size());
    const long M = static_cast<long>(cells_to_integrate.size());
    if (result.rows() != 2 * N || result.cols() != 2 * M)
        throw std::invalid_argument("Matrix size for result is not compatible with the surface meshes in calculation "
                                    "of operator K matrix with two different meshes");

#pragma omp parallel for num_threads(14) collapse(2)
    for (long i = 0; i < N; ++i) {
        for (long j = 0; j < M; ++j) {
            const auto coefs = DiscreteK::getMatrixCoefs(cells[i], cells_to_integrate[j], k);
            result(i, j) = coefs.a11;
            result(i + N, j) = coefs.a21;
            result(i, j + M) = coefs.a12;
            result(i + N, j + M) = coefs.a22;
        }
    }
};

template<typename MatrixLike>
void getMatrixR_inplace(MatrixLike& result, Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh) {
    const auto &cells = surface_mesh.getCells();
    const long N = static_cast<long>(cells.size());
    if (result.rows() != 2 * N || result.cols() != 2 * N)
        throw std::invalid_argument("Matrix size for result is not compatible with the surface mesh in calculation "
                                    "of operator R matrix with one mesh");

#pragma omp parallel for num_threads(14) collapse(2)
    for (long i = 0; i < N; ++i) {
        for (long j = 0; j < N; ++j) {
            const auto coefs = DiscreteR::getMatrixCoefs(i, j, k, cells);
            result(i, j) = coefs.a11;
            result(i + N, j) = coefs.a21;
            result(i, j + N) = coefs.a12;
            result(i + N, j + N) = coefs.a22;
        }
    }
};

template<typename MatrixLike>
void getMatrixR_inplace(MatrixLike& result, Types::complex_d k, const Mesh::SurfaceMesh &integration_mesh,
                           const Mesh::SurfaceMesh &mesh_with_coll_points) {
    const auto &cells = mesh_with_coll_points.getCells();
    const auto &cells_to_integrate = integration_mesh.getCells();
    const long N = static_cast<long>(cells.size());
    const long M = static_cast<long>(cells_to_integrate.size());
    if (result.rows() != 2 * N || result.cols() != 2 * M)
        throw std::invalid_argument("Matrix size for result is not compatible with the surface meshes in calculation "
                                    "of operator R matrix with two different meshes");

#pragma omp parallel for num_threads(14) collapse(2)
    for (long i = 0; i < N; ++i) {
        for (long j = 0; j < M; ++j) {
            const auto coefs = DiscreteR::getMatrixCoefs(cells[i], cells_to_integrate[j], k);
            // if (i == j) std::cout << coefs.a11 << " " << coefs.a12 << " " << coefs.a21 << ' ' << coefs.a22 <<
            // std::endl;
            result(i, j) = coefs.a11;
            result(i + N, j) = coefs.a21;
            result(i, j + M) = coefs.a12;
            result(i + N, j + M) = coefs.a22;
        }
    }
}

};

#endif //MATRIXINPLACEGENERATION_HPP
