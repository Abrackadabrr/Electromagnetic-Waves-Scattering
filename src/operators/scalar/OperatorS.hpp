//
// Created by evgen on 29.09.2025.
//

#ifndef OPERATORS_HPP
#define OPERATORS_HPP

#include "math/Productions.hpp"
#include "math/integration/Quadrature.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"

#include "operators/Functions.hpp"

#include "types/Types.hpp"

#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"

/**
 * Дискретизация скалярного оператора S методом кусочно-постоянных аппроксимаций на некоторой сетке.
 * В данном подоходе оператор S заменяется матрицей, действующей на вектор,
 * соответствующий кусочно-постоянным значениям функций на ячейках сетки.
 *
 * Таким образом схема: функция кусочно-постоянная, оператор слабосингулярный
 */

namespace EMW::Operators {

class S_operator {
    using mesh_t = Mesh::SurfaceMesh;
    using cell_t = Mesh::IndexedCell;
    using point_t = Mesh::point_t;
    using vector_t = Types::VectorXc;
    using matrix_t = Types::MatrixXc;

    Types::complex_d k_;
    Types::scalar h_;

    /**
     * Расчет поверхностного интеграла оператора S по ячейке
     */
    template <typename Quadrature> Types::complex_d S_over_cell(const point_t point, const cell_t &cell) const {
        if ((cell.collPoint_ - point).norm() < h_) {

            const auto phi_bounded = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
                const Types::Vector3d y = cell.parametrization(p, q);
                const Types::scalar mul = cell.multiplier(p, q);
                return Helmholtz::F_bounded_part(k_, point, y) * mul;
            };

            const auto numeric_part = DefiniteIntegrals::integrate<Quadrature>(phi_bounded, {0, 0}, {1., 1.});
            const auto analytic_part = Math::Constants::inverse_4PI<Types::scalar>() *
                                       Math::AnalyticalIntegration::integrate_1_div_r(point, cell);

            return 2. * k_ * k_ * (numeric_part + analytic_part);
        }

        const auto phi_full = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
            const Types::Vector3d y = cell.parametrization(p, q);
            const Types::scalar mul = cell.multiplier(p, q);
            return Helmholtz::F(k_, point, y) * mul;
        };

        return 2. * k_ * k_ * DefiniteIntegrals::integrate<Quadrature>(phi_full, {0, 0}, {1., 1.});
    }

  public:
    S_operator(Types::complex_d k, Types::scalar h): k_(k), h_(h) {};

    matrix_t matrix(const mesh_t &mesh_with_points, const mesh_t &mesh_to_integrate) const {
        const Types::index N = mesh_with_points.getCells().size();
        const Types::index M = mesh_to_integrate.getCells().size();
        Types::MatrixXc result = Types::MatrixXd::Zero(N, M);

        const auto &cells_int = mesh_to_integrate.getCells();
        const auto &cells_col = mesh_with_points.getCells();

        using quadrature_2d = DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>;

#pragma omp parallel for collapse(2)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                result(i, j) = S_over_cell<quadrature_2d>(cells_col[i].collPoint_, cells_int[j]);
            }
        }
        return result;
    };
};

} // namespace EMW::Operators


#endif //OPERATORS_HPP
