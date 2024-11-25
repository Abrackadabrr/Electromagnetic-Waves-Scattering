//
// Created by evgen on 31.01.24.
//

#include "MatrixGeneration.hpp"
#include "math/Productions.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/integration/newton_cotess/Rectangular.hpp"
#include "operators/OperatorK.hpp"
#include "operators/OperatorR.hpp"

#include <iostream>

namespace EMW::Matrix {

namespace DiscreteK {
Types::complex_d getFirstPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                                      const Containers::vector<Mesh::IndexedCell> &cells) {
    return EMW::OperatorK::detail::K1OverSingularCellRnDWithSingularityExtraction<
        DefiniteIntegrals::NewtonCotess::Quadrature<8, 8>>(cells[i].collPoint_, cells[j], k);
}

Types::Matrix3c getZeroPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                                    const Containers::vector<Mesh::IndexedCell> &cells) {
    return EMW::OperatorK::detail::K0TensorOverSingularCell<DefiniteIntegrals::NewtonCotess::Quadrature<8>>(
        cells[i].collPoint_, cells[j], k);
}

MatrixCoefs getMatrixCoefs(Types::index i, Types::index j, Types::complex_d k,
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
} // namespace DiscreteK

namespace DiscreteR {
MatrixCoefs getMatrixCoefs(Types::index i, Types::index j, Types::complex_d k,
                           const Containers::vector<Mesh::IndexedCell> &cells) {
// это мы знаем из свойств оператора R
    if (j == i)
        return {Types::complex_d{0, 0}, Types::complex_d{0, 0}, Types::complex_d{0, 0}, Types::complex_d{0, 0}};
// иначе честно расчитываем
    const Types::Vector3c integral =
        OperatorR::detail::forMatrix::commonIntegralPart<DefiniteIntegrals::NewtonCotess::Quadrature<8, 8>>(
            cells[j], cells[i].collPoint_, k);
    const Types::complex_d a11 = Math::quasiDot(integral, cells[j].tau[0].cross(cells[i].tau[0]));
    const Types::complex_d a12 = Math::quasiDot(integral, cells[j].tau[1].cross(cells[i].tau[0]));
    const Types::complex_d a21 = Math::quasiDot(integral, cells[j].tau[0].cross(cells[i].tau[1]));
    const Types::complex_d a22 = Math::quasiDot(integral, cells[j].tau[1].cross(cells[i].tau[1]));
    return {a11, a12, a21, a22};
}
} // namespace DiscreteR

Types::MatrixXc getMatrixK(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh) {
    const auto &cells = surface_mesh.getCells();
    const long N = static_cast<long>(cells.size());
    Types::MatrixXc result = Types::MatrixXc::Zero(2 * N, 2 * N);

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
    return result;
}

Types::MatrixXc getMatrixR(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh) {
    const auto &cells = surface_mesh.getCells();
    const long N = static_cast<long>(cells.size());
    Types::MatrixXc result = Types::MatrixXc::Zero(2 * N, 2 * N);

#pragma omp parallel for num_threads(1) collapse(2)
    for (long i = 0; i < N; ++i) {
        for (long j = 0; j < N; ++j) {
            const auto coefs = DiscreteR::getMatrixCoefs(i, j, k, cells);
            result(i, j) = coefs.a11;
            result(i + N, j) = coefs.a21;
            result(i, j + N) = coefs.a12;
            result(i + N, j + N) = coefs.a22;
        }
    }
    return result;
}
} // namespace EMW::Matrix
