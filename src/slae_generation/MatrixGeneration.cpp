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

Types::scalar compare(const MatrixCoefs &a, const MatrixCoefs &b) {
    Types::complex_d a11 = a.a11 - b.a11;
    Types::complex_d a12 = a.a12 - b.a12;
    Types::complex_d a21 = a.a21 - b.a21;
    Types::complex_d a22 = a.a22 - b.a22;
    return std::abs(a11) + std::abs(a12) + std::abs(a21) + std::abs(a22);
}

Types::complex_d getFirstPartIntegral(const Mesh::IndexedCell &cell_i, const Mesh::IndexedCell &cell_j,
                                      Types::complex_d k) {
    if ((cell_i.collPoint_ - cell_j.collPoint_).norm() < 1e-15) {
        // std::cout << "Self affecting" << std::endl;
        return EMW::OperatorK::detail::K1OverSingularCellRnDWithSingularityExtraction<
            DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>>(cell_i.collPoint_, cell_j, k);
    }
    return EMW::OperatorK::detail::K1OverSingularCellReducedAndDivided<
        DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>>(cell_i.collPoint_, cell_j, k);
}

Types::complex_d getFirstPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                                      const Containers::vector<Mesh::IndexedCell> &cells) {
    return getFirstPartIntegral(cells[i], cells[j], k);
}

Types::Matrix3c getZeroPartIntegral(const Mesh::IndexedCell &cell_i, const Mesh::IndexedCell &cell_j,
                                    Types::complex_d k) {
    return EMW::OperatorK::detail::K0TensorOverSingularCell<DefiniteIntegrals::GaussLegendre::Quadrature<8>>(
        cell_i.collPoint_, cell_j, k);
}

Types::Matrix3c getZeroPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                                    const Containers::vector<Mesh::IndexedCell> &cells) {
    return getZeroPartIntegral(cells[i], cells[j], k);
}

MatrixCoefs getMatrixCoefs(const Mesh::IndexedCell &cell_i, const Mesh::IndexedCell &cell_j, Types::complex_d k) {
    const Types::Matrix3c int0 = getZeroPartIntegral(cell_i, cell_j, k);
    const Types::complex_d int1_k2 = k * k * getFirstPartIntegral(cell_i, cell_j, k);

    const Types::complex_d a11_0 = cell_i.tau[0].transpose() * int0 * cell_j.tau[0];
    const Types::complex_d a12_0 = cell_i.tau[0].transpose() * int0 * cell_j.tau[1];
    const Types::complex_d a21_0 = cell_i.tau[1].transpose() * int0 * cell_j.tau[0];
    const Types::complex_d a22_0 = cell_i.tau[1].transpose() * int0 * cell_j.tau[1];

    const Types::complex_d a11_1 = Math::quasiDot(cell_i.tau[0], cell_j.tau[0]) * int1_k2;
    const Types::complex_d a12_1 = Math::quasiDot(cell_i.tau[0], cell_j.tau[1]) * int1_k2;
    const Types::complex_d a21_1 = Math::quasiDot(cell_i.tau[1], cell_j.tau[0]) * int1_k2;
    const Types::complex_d a22_1 = Math::quasiDot(cell_i.tau[1], cell_j.tau[1]) * int1_k2;

    return {a11_0 + a11_1, a12_0 + a12_1, a21_0 + a21_1, a22_0 + a22_1};
}

Containers::array<Types::complex_d, 4> getMatrixCoefsInArray(const Mesh::IndexedCell &cell_i,
                                                             const Mesh::IndexedCell &cell_j, Types::complex_d k) {
    const auto res = getMatrixCoefs(cell_i, cell_j, k);
    return {res.a11, res.a12, res.a21, res.a22};
}

MatrixCoefs getMatrixCoefs(Types::index i, Types::index j, Types::complex_d k,
                           const Containers::vector<Mesh::IndexedCell> &cells) {
    return getMatrixCoefs(cells[i], cells[j], k);
}
} // namespace DiscreteK

namespace DiscreteR {

MatrixCoefs getMatrixCoefs(const Mesh::IndexedCell &cell_i, const Mesh::IndexedCell &cell_j, Types::complex_d k) {
    // честно рассчитываем
    const Types::Vector3c integral =
        OperatorR::detail::forMatrix::commonIntegralPart<DefiniteIntegrals::NewtonCotess::Quadrature<2, 2>>(
            cell_j, cell_i.collPoint_, k);
#if 1
    const Types::complex_d a11 = Math::quasiDot(integral, cell_j.tau[0].cross(cell_i.tau[0]));
    const Types::complex_d a12 = Math::quasiDot(integral, cell_j.tau[1].cross(cell_i.tau[0]));
    const Types::complex_d a21 = Math::quasiDot(integral, cell_j.tau[0].cross(cell_i.tau[1]));
    const Types::complex_d a22 = Math::quasiDot(integral, cell_j.tau[1].cross(cell_i.tau[1]));
#endif
#if 0
    const Types::complex_d a11 = Math::quasiDot(cell_i.tau[0], Math::cross(integral, cell_j.tau[0]));
    const Types::complex_d a12 = Math::quasiDot(cell_i.tau[0], Math::cross(integral, cell_j.tau[1]));
    const Types::complex_d a21 = Math::quasiDot(cell_i.tau[1], Math::cross(integral, cell_j.tau[0]));
    const Types::complex_d a22 = Math::quasiDot(cell_i.tau[1], Math::cross(integral, cell_j.tau[1]));
#endif
    return {a11, a12, a21, a22};
}

Containers::array<Types::complex_d, 4> getMatrixCoefsInArray(const Mesh::IndexedCell &cell_i,
                                                             const Mesh::IndexedCell &cell_j, Types::complex_d k) {
    const auto res = getMatrixCoefs(cell_i, cell_j, k);
    return {res.a11, res.a12, res.a21, res.a22};
}

MatrixCoefs getMatrixCoefs(Types::index i, Types::index j, Types::complex_d k,
                           const Containers::vector<Mesh::IndexedCell> &cells) {
    // это мы знаем из свойств оператора R
    if (j == i)
        return {Types::complex_d{0, 0}, Types::complex_d{0, 0}, Types::complex_d{0, 0}, Types::complex_d{0, 0}};
    // если i == j, то имеем коэффициент, который описывает воздействие на себя
    // это означает совпадение не только точек коллокации, но и плоскостей ячеек
    // именно из-за совпадения плоскостей мы получаем нуль в результате интегрирования

    // иначе честно рассчитываем
    return getMatrixCoefs(cells[i], cells[j], k);
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

Types::MatrixXc getMatrixK(Types::complex_d k, const Mesh::SurfaceMesh &integration_mesh,
                           const Mesh::SurfaceMesh &mesh_with_coll_points) {
    const auto &cells = mesh_with_coll_points.getCells();
    const auto &cells_to_integrate = integration_mesh.getCells();
    const long N = static_cast<long>(cells.size());
    const long M = static_cast<long>(cells_to_integrate.size());
    Types::MatrixXc result = Types::MatrixXc::Zero(2 * N, 2 * M);

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
    return result;
}

Types::MatrixXc getMatrixR(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh) {
    const auto &cells = surface_mesh.getCells();
    const long N = static_cast<long>(cells.size());
    Types::MatrixXc result = Types::MatrixXc::Zero(2 * N, 2 * N);

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
    return result;
}

Types::MatrixXc getMatrixR(Types::complex_d k, const Mesh::SurfaceMesh &integration_mesh,
                           const Mesh::SurfaceMesh &mesh_with_coll_points) {
    const auto &cells = mesh_with_coll_points.getCells();
    const auto &cells_to_integrate = integration_mesh.getCells();
    const long N = static_cast<long>(cells.size());
    const long M = static_cast<long>(cells_to_integrate.size());
    Types::MatrixXc result = Types::MatrixXc::Zero(2 * N, 2 * M);

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
    return result;
}

Types::MatrixXc getMatrixCrossNormal(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh) {
    const long N = static_cast<long>(surface_mesh.getCells().size());
    Types::MatrixXc result = Types::MatrixXc::Zero(2 * N, 2 * N);

    result.block(0, N, N, N) = Types::MatrixXc::Identity(N, N);
    result.block(N, 0, N, N) = -Types::MatrixXc::Identity(N, N);

    return result;
}

Types::MatrixXc getMatrixIdentity(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh) {
    return Types::MatrixXc::Identity(2 * surface_mesh.getCells().size(), 2 * surface_mesh.getCells().size());
}

} // namespace EMW::Matrix
