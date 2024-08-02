//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP

#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/MeshTypes.hpp"

namespace EMW::Matrix {

    /**
     * Contain K0 or K1 parts of matrix coefficients
     * first figure means number of tau_i (i - number of an equation)
     * second figure means number of tau_j (j - number of a part in sum)
     */
    struct MatrixCoefs {
        Types::complex_d a11;
        Types::complex_d a12;
        Types::complex_d a21;
        Types::complex_d a22;
    };

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
    Types::complex_d
    getFirstPartIntegral(Types::index i, Types::index j, Types::scalar k,
                         const Containers::vector<Mesh::IndexedCell> &cells);


    /**
     * Returning the zero part of matrix coefficients (contour integral)
     * @param i - index of an equation
     * @param j - index inside the equation
     * @param k - wave number
     * @param cells - cells in the mesh (with collocation nodes)
     * @return
     */
    Types::Matrix3c
    getZeroPartIntegral(Types::index i, Types::index j, Types::scalar k,
                        const Containers::vector<Mesh::IndexedCell> &cells);

    MatrixCoefs
    getMatrixCoefs(Types::index i, Types::index j, Types::scalar k,
                   const Containers::vector<Mesh::IndexedCell> &cells);

    Types::MatrixXc getMatrix(Types::scalar k, const Mesh::SurfaceMesh& surface_mesh);

    Types::VectorXc getRHS(const Types::Vector3d &pol, const Types::Vector3d &k_vec,
                           const Containers::vector<Mesh::IndexedCell> &cells);
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP
