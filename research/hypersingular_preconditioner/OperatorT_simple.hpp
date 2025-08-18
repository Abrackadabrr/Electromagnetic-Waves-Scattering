//
// Created by evgen on 06.08.2025.
//

#ifndef OPERATORT_SIMPLE_HPP
#define OPERATORT_SIMPLE_HPP

#include "math/Productions.hpp"
#include "math/integration/CurveIntegrals.hpp"
#include "math/integration/Quadrature.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"

#include "operators/Functions.hpp"

#include "types/Types.hpp"

#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"

/**
 * Дискретизация скалярного оператора T методом кусочно-постоянных аппроксимаций на некоторой сетке.
 * В данном подоходе оператор T заменяется матрицей, действующей на вектор из констант,
 * соответствующим кусочно-постоянным значениям функций на ячейках сетки.
 *
 * Таким образом схема: функция кусочно-постоянная, оператор гиперсингулярный
 */
namespace EMW::helper_fucntions {
inline Types::Vector3d field_to_circulate(const Types::Vector3d &x, const Types::Vector3d &y) {
    const auto gradF = Laplace::gradF(x, y);
    return {gradF[1], -gradF[0], 0};
}
} // namespace EMW::helper_fucntions

namespace EMW::OperatorT::LaplaceEquation {
template <typename Quadrature>
Types::scalar T_over_single_cell(const Types::Vector3d point, const Mesh::IndexedCell &cell) {
    // Задаем функции одной переменной: векторное поле на краю ячейки
    const auto AB = [&](Types::scalar t) -> Types::Vector3d {
        const Types::Vector3d y = cell.parametrization(t, 0);
        return Laplace::gradF(point, y);
    };
    const auto BC = [&](Types::scalar t) -> Types::Vector3d {
        const Types::Vector3d y = cell.parametrization(1, t);
        return Laplace::gradF(point, y);
    };
    const auto CD = [&](Types::scalar t) -> Types::Vector3d {
        const Types::Vector3d y = cell.parametrization(1 - t, 1);
        return Laplace::gradF(point, y);
    };
    const auto DA = [&](Types::scalar t) -> Types::Vector3d {
        const Types::Vector3d y = cell.parametrization(0, 1 - t);
        return Laplace::gradF(point, y);
    };
    // Интегрируем вектор-функции по краю, а зачем домножаем на касательный вектор (так как он всюду одинаков)
    return DefiniteIntegrals::integrate<Quadrature>(AB, {0}, {1}).dot(cell.integrationParameters.mul[0]) +
           DefiniteIntegrals::integrate<Quadrature>(BC, {0}, {1}).dot(cell.integrationParameters.mul[1]) +
           DefiniteIntegrals::integrate<Quadrature>(CD, {0}, {1}).dot(cell.integrationParameters.mul[2]) +
           DefiniteIntegrals::integrate<Quadrature>(DA, {0}, {1}).dot(cell.integrationParameters.mul[3]);
}

Types::MatrixXd T_over_mesh(const Mesh::SurfaceMesh &mesh_with_collocation_point,
                            const Mesh::SurfaceMesh &mesh_to_integrate) {
    const Types::index N = mesh_with_collocation_point.getCells().size();
    const Types::index M = mesh_to_integrate.getCells().size();
    Types::MatrixXd result = Types::MatrixXd::Zero(N, M);

    const auto &cells_int = mesh_to_integrate.getCells();
    const auto &cells_col = mesh_with_collocation_point.getCells();

#pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            result(i, j) = T_over_single_cell<EMW::DefiniteIntegrals::GaussLegendre::Quadrature<5>>(
                cells_col[i].collPoint_, cells_int[j]);
        }
    }
    return result;
}
} // namespace EMW::OperatorT::LaplaceEquation

namespace EMW::OperatorT::HelmholtzEquation {
template <typename Quadrature>
Types::complex_d T_curved_over_single_cell(Types::complex_d k, const Types::Vector3d point,
                                           const Mesh::IndexedCell &cell) {
    // Задаем функции одной переменной: векторное поле на краю ячейки
    const auto AB = [&](Types::scalar t) -> Types::Vector3c {
        const Types::Vector3d y = cell.parametrization(t, 0);
        return -Helmholtz::V(k, point, y);
    };
    const auto BC = [&](Types::scalar t) -> Types::Vector3c {
        const Types::Vector3d y = cell.parametrization(1, t);
        return -Helmholtz::V(k, point, y);
    };
    const auto CD = [&](Types::scalar t) -> Types::Vector3c {
        const Types::Vector3d y = cell.parametrization(1 - t, 1);
        return -Helmholtz::V(k, point, y);
    };
    const auto DA = [&](Types::scalar t) -> Types::Vector3c {
        const Types::Vector3d y = cell.parametrization(0, 1 - t);
        return -Helmholtz::V(k, point, y);
    };
    // Интегрируем вектор-функции по краю, а зачем домножаем на касательный вектор (так как он всюду одинаков)
    return Math::quasiDot(DefiniteIntegrals::integrate<Quadrature>(AB, {0}, {1}), cell.integrationParameters.mul[0]) +
           Math::quasiDot(DefiniteIntegrals::integrate<Quadrature>(BC, {0}, {1}), cell.integrationParameters.mul[1]) +
           Math::quasiDot(DefiniteIntegrals::integrate<Quadrature>(CD, {0}, {1}), cell.integrationParameters.mul[2]) +
           Math::quasiDot(DefiniteIntegrals::integrate<Quadrature>(DA, {0}, {1}), cell.integrationParameters.mul[3]);
}

template <typename Quadrature>
Types::complex_d T_surfaced_over_single_cell(Types::complex_d k, const Types::Vector3d point,
                                             const Mesh::IndexedCell &cell) {
    const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
        const Types::Vector3d y = cell.parametrization(p, q);
        const Types::scalar mul = cell.multiplier(p, q);
        return Helmholtz::F_bounded_part(k, point, y) * mul;
    };
    return k * k *
           (DefiniteIntegrals::integrate<Quadrature>(phi, {0, 0}, {1., 1.}) +
            Math::Constants::inverse_4PI<Types::scalar>() *
                Math::AnalyticalIntegration::integrate_1_div_r(point, cell));
}

Types::MatrixXc T_over_mesh(Types::complex_d k, const Mesh::SurfaceMesh &mesh_with_collocation_point,
                            const Mesh::SurfaceMesh &mesh_to_integrate) {
    using quadrature_1d = DefiniteIntegrals::GaussLegendre::Quadrature<5>;
    using quadrature_2d = DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>;
    Types::index N = mesh_with_collocation_point.getCells().size();
    const Types::index M = mesh_to_integrate.getCells().size();
    Types::MatrixXc result = Types::MatrixXc::Zero(N, M);

    const auto& cells_int = mesh_to_integrate.getCells();
    const auto& cells_col = mesh_with_collocation_point.getCells();

#pragma omp parallel for collapse(2)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                result(i,j) = T_curved_over_single_cell<quadrature_1d>(k, cells_col[i].collPoint_, cells_int[j]) -
                                     T_surfaced_over_single_cell<quadrature_2d>(k, cells_col[i].collPoint_, cells_int[j]);
            }
        }
        return result;
    }
}


#endif //OPERATORT_SIMPLE_HPP
