//
// Created by evgen on 31.01.24.
//

#include "MatrixGeneration.hpp"
#include "operators/Functions.hpp"
#include "operators/Operators.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/Productions.hpp"

namespace EMW::Matrix {
    Types::complex_d
    getFirstPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                         const Containers::vector<Mesh::IndexedCell> &cells) {
        const Types::scalar critical_rho = 0.1;
        return (cells[i].collPoint_.point_ - cells[i].collPoint_.point_).norm() > critical_rho ?
               EMW::Operators::detail::K1OverSingularCellRnDWithSingularityExtraction<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(
                       cells[i].collPoint_.point_, cells[j], k)
                      :
               EMW::Operators::detail::K1OverSingularCellReducedAndDivided<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(
                       cells[i].collPoint_.point_, cells[j], k);
    }

    Types::Matrix3c
    getZeroPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                        const Containers::vector<Mesh::IndexedCell> &cells) {
        return
//        i == j ? Types::Matrix3c::Zero() :
        EMW::Operators::detail::K0TensorOverSingularCell<DefiniteIntegrals::GaussLegendre::Quadrature<8>>(
                cells[i].collPoint_.point_, cells[j], k);
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

        const Types::complex_d a11_1 = Math::quasiDot(cells[i].tau[0], cells[j].tau[0]) * int1_k2;
        const Types::complex_d a12_1 = Math::quasiDot(cells[i].tau[0], cells[j].tau[1]) * int1_k2;
        const Types::complex_d a21_1 = Math::quasiDot(cells[i].tau[1], cells[j].tau[0]) * int1_k2;
        const Types::complex_d a22_1 = Math::quasiDot(cells[i].tau[1], cells[j].tau[1]) * int1_k2;

        return {a11_0 + a11_1, a12_0 + a12_1, a21_0 + a21_1, a22_0 + a22_1};
    }

    Types::MatrixXc getMatrix(Types::complex_d k, const Mesh::SurfaceMesh& surface_mesh) {
        const auto & cells = surface_mesh.getCells();
        const long N = static_cast<long>(cells.size());
        Types::MatrixXc result = Types::MatrixXc::Zero(2 * N, 2 * N);

#pragma omp parallel for num_threads(14) collapse(2)
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
