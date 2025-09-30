//
// Created by evgen on 29.09.2025.
//

#ifndef OPERATORK_HPP
#define OPERATORK_HPP

#include "math/Productions.hpp"
#include "math/integration/Quadrature.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"

#include "operators/Functions.hpp"

#include "types/Types.hpp"

#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"

/**
 * Дискретизация скалярного сингулярного оператора K методом кусочно-постоянных аппроксимаций на некоторой сетке.
 * В данном подоходе оператор K заменяется матрицей, действующей на вектор,
 * соответствующий кусочно-постоянным значениям функций на ячейках сетки.
 *
 * Таким образом схема: функция кусочно-постоянная, оператор сингулярный
 */

namespace EMW::Operators {

class K_operator {
    using mesh_t = Mesh::SurfaceMesh;
    using cell_t = Mesh::IndexedCell;
    using point_t = Mesh::point_t;
    using vector_t = Types::VectorXc;
    using matrix_t = Types::MatrixXc;

    Types::complex_d k_;

    template <typename Quadrature> Types::complex_d K_over_cell(const point_t &point, const cell_t &cell) const {
        bool point_inside_cell =
            (point - cell.collPoint_).norm() < 1e-6; // тут должна быть нормальная функция, пока косыль
        if (point_inside_cell)
            return {0, 0};

        const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
            const Types::Vector3d y = cell.parametrization(p, q);
            const Types::scalar mul = cell.multiplier(p, q);
            return Math::quasiDot(Helmholtz::V(k_, point, y), cell.normal) * mul;
        };

        return -2. * DefiniteIntegrals::integrate<Quadrature>(phi, {0, 0}, {1, 1});
    }

  public:
    K_operator(Types::complex_d k) : k_(k){};

    matrix_t matrix(const mesh_t &mesh_with_points, const mesh_t &mesh_to_integrate) const {
        using quadrature_2d = DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>;
        Types::index N = mesh_with_points.getCells().size();
        const Types::index M = mesh_to_integrate.getCells().size();
        Types::MatrixXc result = Types::MatrixXc::Zero(N, M);

        const auto &cells_int = mesh_to_integrate.getCells();
        const auto &cells_col = mesh_with_points.getCells();

#pragma omp parallel for collapse(2)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                result(i, j) = K_over_cell<quadrature_2d>(cells_col[i].collPoint_, cells_int[j]);
            }
        }
        return result;
    }

};
}

#endif //OPERATORK_HPP
