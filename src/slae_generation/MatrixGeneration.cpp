//
// Created by evgen on 31.01.24.
//

#include "MatrixGeneration.hpp"
#include "slae_generation/Functions.hpp"
#include "integration/gauss_quadrature/GaussQuadrature.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"

namespace EMW::Matrix {
    Types::complex_d
    getFirstPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                         const Containers::vector<Mesh::IndexedCell> &cells) {
        const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
            const Types::Vector3d y = cells[j].parametrization(p, q);
            const Types::scalar mul = cells[j].multiplier(p, q);
            return Helmholtz::F(k, cells[i].collPoint_.point_, y) * mul;
        };
        return DefiniteIntegrals::integrateImproper<DefiniteIntegrals::Quadrature<20, 20>>(phi, {0, 0}, {1./2, 1./2}) +
                DefiniteIntegrals::integrateImproper<DefiniteIntegrals::Quadrature<20, 20>>(phi, {0, 1./2}, {1./2, 1./2}) +
                DefiniteIntegrals::integrateImproper<DefiniteIntegrals::Quadrature<20, 20>>(phi, {1./2, 0}, {1./2, 1./2}) +
                DefiniteIntegrals::integrateImproper<DefiniteIntegrals::Quadrature<20, 20>>(phi, {1./2, 1./2}, {1./2, 1./2});
    }

    Types::Matrix3c
    getZeroPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                        const Containers::vector<Mesh::IndexedCell> &cells) {
        const auto AB = [&](Types::scalar t) -> Types::Matrix3c {
            const Types::Vector3d y = cells[j].parametrization(t, 0);
            return Helmholtz::V(k, cells[i].collPoint_.point_, y) * cells[j].integrationParameters.mul[0].transpose();
        };
        const auto BC = [&](Types::scalar t) -> Types::Matrix3c {
            const Types::Vector3d y = cells[j].parametrization(1, t);
            return Helmholtz::V(k, cells[i].collPoint_.point_, y) * cells[j].integrationParameters.mul[1].transpose();
        };
        const auto CD = [&](Types::scalar t) -> Types::Matrix3c {
            const Types::Vector3d y = cells[j].parametrization(1 - t, 1);
            return Helmholtz::V(k, cells[i].collPoint_.point_, y) * cells[j].integrationParameters.mul[2].transpose();
        };
        const auto DA = [&](Types::scalar t) -> Types::Matrix3c {
            const Types::Vector3d y = cells[j].parametrization(0, 1 - t);
            return Helmholtz::V(k, cells[i].collPoint_.point_, y) * cells[j].integrationParameters.mul[3].transpose();
        };

        return DefiniteIntegrals::integrate<DefiniteIntegrals::Quadrature<8>>(AB, {0}, {1}) +
               DefiniteIntegrals::integrate<DefiniteIntegrals::Quadrature<8>>(BC, {0}, {1}) +
               DefiniteIntegrals::integrate<DefiniteIntegrals::Quadrature<8>>(CD, {0}, {1}) +
               DefiniteIntegrals::integrate<DefiniteIntegrals::Quadrature<8>>(DA, {0}, {1});
    }

    MatrixCoefs
    getMatrixCoefs(Types::index i, Types::index j, Types::complex_d k,
                   const Containers::vector<Mesh::IndexedCell> &cells) {
        const Types::Matrix3c int0 = getZeroPartIntegral(i, j, k, cells);
        const Types::complex_d int1_k2 = k * k * getFirstPartIntegral(i, j, k, cells);

        const Types::complex_d a11_0 = cells[i].tau[0].transpose() * int0 * cells[j].tau[0];
        const Types::complex_d a12_0 = cells[i].tau[0].transpose() * int0 * cells[j].tau[1];
        const Types::complex_d a21_0 = cells[i].tau[1].transpose() * int0 * cells[j].tau[0];
        const Types::complex_d a22_0 = cells[i].tau[1].transpose() * int0 * cells[j].tau[1];

        const Types::complex_d a11_1 = cells[i].tau[0].dot(cells[j].tau[0]) * int1_k2;
        const Types::complex_d a12_1 = cells[i].tau[0].dot(cells[j].tau[1]) * int1_k2;
        const Types::complex_d a21_1 = cells[i].tau[1].dot(cells[j].tau[0]) * int1_k2;
        const Types::complex_d a22_1 = cells[i].tau[1].dot(cells[j].tau[1]) * int1_k2;

        return {a11_0 + a11_1, a12_0 + a12_1, a21_0 + a21_1, a22_0 + a22_1};
    }

    Types::MatrixXc getMatrix(Types::complex_d k, const Containers::vector<Mesh::IndexedCell> &cells) {
        const long N = static_cast<long>(cells.size());
        Types::MatrixXc result = Types::MatrixXc::Zero(2 * N, 2 * N);

        for (long i = 0; i < N; ++i) {
            for (long j = 0; j < N; ++j) {
                const auto coefs = getMatrixCoefs(i, j, k, cells);
                result(i, j) = coefs.a11;
                result(i + N, j) = coefs.a21;
                result(i, j + N) = coefs.a12;
                result(i + N, j + N) = coefs.a22;
            }
        }
        return result;
    }
}