//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP

#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "types/Types.hpp"

namespace EMW::Matrix {

enum class operator_t {
  R = 0,
  K = 1
};

/**
 * Contain parts of matrix coefficients
 * first figure means number of tau_i (i - number of an equation)
 * second figure means number of tau_j (j - number of a part in sum)
 */
struct MatrixCoefs {
    Types::complex_d a11;
    Types::complex_d a12;
    Types::complex_d a21;
    Types::complex_d a22;
};

namespace DiscreteK {
struct ContourIntegralParts {
    Types::Vector3c ab;
    Types::Vector3c bc;
    Types::Vector3c cd;
    Types::Vector3c da;
};

/***
 * Returning the first part of matrix coefficient (surface integral) divided by k^2
 * (aka dot product of K_1{ tau[m]_j, \sigma_j} and tau[m]_i)
 * @param i - index of an equation
 * @param j - index inside the equation
 * @param k - wave number
 * @param cells - cells in the mesh (with collocation nodes)
 * @return
 */
Types::complex_d getFirstPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                                      const Containers::vector<Mesh::IndexedCell> &cells);

/**
 * Returning the zero part of matrix coefficients (contour integral)
 * @param i - index of an equation
 * @param j - index inside the equation
 * @param k - wave number
 * @param cells - cells in the mesh (with collocation nodes)
 * @return
 */
Types::Matrix3c getZeroPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                                    const Containers::vector<Mesh::IndexedCell> &cells);

MatrixCoefs getMatrixCoefs(Types::index i, Types::index j, Types::complex_d k,
                           const Containers::vector<Mesh::IndexedCell> &cells);
} // namespace DiscreteK

namespace DiscreteR {
MatrixCoefs getMatrixCoefs(Types::index i, Types::index j, Types::complex_d k,
                           const Containers::vector<Mesh::IndexedCell> &cells);
}

Types::MatrixXc getMatrixK(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh);

Types::MatrixXc getMatrixR(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh);

} // namespace EMW::Matrix

#endif // ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP
