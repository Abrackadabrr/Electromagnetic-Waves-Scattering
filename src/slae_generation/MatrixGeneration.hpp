//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP

#include "types/Types.hpp"
#include "mesh/MeshTypes.hpp"
#include "slae_generation/Functions.hpp"
#include "integration/gauss_quadrature/GaussQuadrature.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"

namespace EMW::Matrix {

    /**
     * Contain K1 parts of matrix coeffitients
     * first figure means number of tau_i (i - number of an equation)
     * second figure means number of tau_j (j - number of a part in sum)
     */
    struct K1_coefs {
        Types::complex_d a11;
        Types::complex_d a12;
        Types::complex_d a21;
        Types::complex_d a22;
    };

    /***
     * Returning the first part of matrix coefficient divided by k^2
     * (aka dot product of K_1{ tau[m]_j, \sigma_j} and tau[m]_i)
     * @param i - index of an equation
     * @param j - index inside the equation
     * @param k - wave number
     * @param cells - cells in the mesh (with collocation nodes)
     * @return K1_coefs -- a_{ij}^{(1, 1), (1, 2) (2, 1) (2, 2)} / k^2 -- four coeffitients in matrix
     */
    K1_coefs
    getFirstPart(Types::index i, Types::index j, Types::complex_d k,
                 const Containers::vector<Mesh::IndexedCell> &cells) {
        const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
            const Types::Vector3d y = cells[j].parametrization(p, q);
            const Types::scalar mul = cells[j].multiplier(p, q);
            return Helmholtz::F(k, cells[i].collPoint_.point_, y) * mul;
        };
        const Types::scalar c11 = cells[i].tau[0].dot(cells[j].tau[0]);
        const Types::scalar c12 = cells[i].tau[0].dot(cells[j].tau[1]);
        const Types::scalar c21 = cells[i].tau[1].dot(cells[j].tau[0]);
        const Types::scalar c22 = cells[i].tau[1].dot(cells[j].tau[1]);
        const Types::complex_d integral = DefiniteIntegrals::integrate<DefiniteIntegrals::Quadrature<6, 6>>(phi,
                                                                                                            {-1, -1},
                                                                                                            {2, 2});
        return {c11 * integral, c12 * integral, c21 * integral, c22 * integral};
    }

    /**
     *
     * @param i
     * @param j
     * @param p
     * @param m
     * @param k
     * @param cells
     * @return
     */
    Types::complex_d
    getZeroPart(Types::index i, Types::index j, Types::index p, Types::index m, Types::complex_d k,
                const Containers::vector<Mesh::IndexedCell> &cells) {
        return {0, 0}
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP
