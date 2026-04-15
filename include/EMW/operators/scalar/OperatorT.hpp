//
// Created by evgen on 06.08.2025.
//

#ifndef OPERATORT_HPP
#define OPERATORT_HPP

#include "math/Productions.hpp"
#include "math/integration/Quadrature.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"

#include "operators/Functions.hpp"

#include "types/Types.hpp"

#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"


/**
 * Дискретизация скалярного оператора T методом кусочно-постоянных аппроксимаций на некоторой сетке.
 * В данном подоходе оператор T заменяется матрицей, действующей на вектор,
 * соответствующий кусочно-постоянным значениям функций на ячейках сетки.
 *
 * Таким образом схема: функция кусочно-постоянная, метод коллокаций
 */
namespace EMW::Operators {

enum OperatorType { LAPLACE, HELMHOLTZ };

template <OperatorType OP = HELMHOLTZ> class T_operator {
    using mesh_t = Mesh::SurfaceMesh;
    using cell_t = Mesh::IndexedCell;
    using point_t = Mesh::point_t;
    using vector_t = Types::VectorXc;
    using matrix_t = Types::MatrixXc;

    Types::complex_d k_;
    Types::scalar h_;

    /**
     *  Функция для расчета криволинейной части для оператора T на одной ячейке
     */
    template <typename Quadrature>
    [[nodiscard]] Types::complex_d T_contour_part_over_cell(const point_t &point, const cell_t &cell, const Types::Vector3d &n_x,
                                             const Types::Vector3d &n_y) const {
        // 1. Задаем функции одной переменной: векторное поле на краю ячейки
        const auto AB = [&](Types::scalar t) -> Types::Vector3c {
            const Types::Vector3d y = cell.parametrization(t, 0);
            return Helmholtz::V(k_, point, y);
        };
        const auto BC = [&](Types::scalar t) -> Types::Vector3c {
            const Types::Vector3d y = cell.parametrization(1, t);
            return Helmholtz::V(k_, point, y);
        };
        const auto CD = [&](Types::scalar t) -> Types::Vector3c {
            const Types::Vector3d y = cell.parametrization(1 - t, 1);
            return Helmholtz::V(k_, point, y);
        };
        const auto DA = [&](Types::scalar t) -> Types::Vector3c {
            const Types::Vector3d y = cell.parametrization(0, 1 - t);
            return Helmholtz::V(k_, point, y);
        };

        // 2. Интегрируем -\nabla_x(F) по краю

        const Types::Vector3c v1 = DefiniteIntegrals::integrate<Quadrature>(AB, {0}, {1});
        const Types::Vector3c v2 = DefiniteIntegrals::integrate<Quadrature>(BC, {0}, {1});
        const Types::Vector3c v3 = DefiniteIntegrals::integrate<Quadrature>(CD, {0}, {1});
        const Types::Vector3c v4 = DefiniteIntegrals::integrate<Quadrature>(DA, {0}, {1});

        // 3. Домножаем на нормальные вектора к контуру в соотвествующем интеграле
        // и расчитываем два слагаемых, как в формулке для
        const auto &nu1 = cell.integrationParameters.mul[0];
        const auto &nu2 = cell.integrationParameters.mul[1];
        const auto &nu3 = cell.integrationParameters.mul[2];
        const auto &nu4 = cell.integrationParameters.mul[3];

        const Types::complex_d part_1 = n_y.dot(n_x) * (Math::quasiDot(v1, nu1) + Math::quasiDot(v2, nu2) +
                                                        Math::quasiDot(v3, nu3) + Math::quasiDot(v4, nu4));
        const Types::complex_d part_2 = Math::quasiDot(v1, n_y) * nu1.dot(n_x) +
                                        Math::quasiDot(v2, n_y) * nu2.dot(n_x) +
                                        Math::quasiDot(v3, n_y) * nu3.dot(n_x) + Math::quasiDot(v4, n_y) * nu4.dot(n_x);

        // 4. Итоговый ответ для контурного интеграла
        return part_1 - part_2;
    }

    /**
     *  Функция для расчета поверхностной части для оператора T на одной ячейке
     */
    template <typename Quadrature>
    [[nodiscard]] Types::complex_d T_surface_part_over_cell(const point_t point, const cell_t &cell) const {
        if ((cell.collPoint_ - point).norm() < h_) {

            const auto phi_bounded = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
                const Types::Vector3d y = cell.parametrization(p, q);
                const Types::scalar mul = cell.multiplier(p, q);
                return Helmholtz::F_bounded_part(k_, point, y) * mul;
            };

            const auto numeric_part = DefiniteIntegrals::integrate<Quadrature>(phi_bounded, {0, 0}, {1., 1.});
            const auto analytic_part = Math::Constants::inverse_4PI<Types::scalar>() *
                                       Math::AnalyticalIntegration::integrate_1_div_r(point, cell);

            return k_ * k_ * (numeric_part + analytic_part);
        }

        const auto phi_full = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
            const Types::Vector3d y = cell.parametrization(p, q);
            const Types::scalar mul = cell.multiplier(p, q);
            return Helmholtz::F(k_, point, y) * mul;
        };

        return k_ * k_ * DefiniteIntegrals::integrate<Quadrature>(phi_full, {0, 0}, {1., 1.});
    }

  public:
    T_operator(Types::complex_d k, Types::scalar h) : k_(k), h_(h){};

    [[nodiscard]] matrix_t matrix(const mesh_t &mesh_with_point, const mesh_t &mesh_to_integrate) const {
        using quadrature_1d = DefiniteIntegrals::GaussLegendre::Quadrature<5>;
        using quadrature_2d = DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>;
        Types::index N = mesh_with_point.getCells().size();
        const Types::index M = mesh_to_integrate.getCells().size();
        Types::MatrixXc result = Types::MatrixXc::Zero(N, M);

        const auto &cells_int = mesh_to_integrate.getCells();
        const auto &cells_col = mesh_with_point.getCells();

#pragma omp parallel for collapse(2)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                const auto& n_x = cells_col[i].normal;
                const auto& n_y = cells_int[j].normal;
                const auto fp = T_contour_part_over_cell<quadrature_1d>(cells_col[i].collPoint_, cells_int[j], n_x, n_y);
                const auto sp = n_x.dot(n_y) * T_surface_part_over_cell<quadrature_2d>(cells_col[i].collPoint_, cells_int[j]);
                result(i, j) = 2. * (fp + sp);
            }
        }
        return result;
    }
};

} // namespace EMW::Operators
#endif //OPERATORT_HPP
