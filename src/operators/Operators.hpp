//
// Created by evgen on 13.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_OPERATORS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_OPERATORS_HPP

#include <types/Types.hpp>
#include "mesh/MeshTypes.hpp"
#include "Functions.hpp"
#include "integration/Quadrature.hpp"
#include "mesh/Mesh.hpp"

namespace EMW::Operators {
    namespace detail {
        /**
         * Расчет интергального оператора К1/{k^2} на одной ячейке без умножения на вектор
         * @tparam Quadrature (тип) точки квадратуры
         * @tparam cell_t тип ячейки (Cell, IndexedCell)
         * @param point точка для расчета
         * @param cell ячейка, на которой считаем
         * @param k волновое число
         * @return интеграл от F(x - y) по ячейке cell
         */
        template<typename Quadrature, typename cell_t>
        Types::complex_d
        K1OverSingularReducedAndDivided(const Mesh::Point &point, const cell_t &cell, const Types::scalar k) {
            const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
                const Types::Vector3d y = cell.parametrization(p, q);
                const Types::scalar mul = cell.multiplier(p, q);
                return Helmholtz::F(k, point, y) * mul;
            };
            return DefiniteIntegrals::integrate<Quadrature>(phi, {0, 0}, {1., 1.});
        }

        /**
         * Расчет множеителя - интергала по одной ячейке в операторе К0 в виде тензора
         * @tparam Quadrature тип квадратуры
         * @tparam cell_t тип ячейки
         * @param point точка для расчета
         * @param cell ячейка
         * @param k волновое число
         * @return тензор, умножение на вектор происходит справа
         */
        template<typename Quadrature, typename cell_t>
        Types::Matrix3c
        K0TensorOverSingularCell(const Mesh::Point &point, const cell_t &cell, Types::scalar k) {
            const auto AB = [&](Types::scalar t) -> Types::Vector3c {
                const Types::Vector3d y = cell.parametrization(t, 0);
                return Helmholtz::V(k, point, y);
            };
            const auto BC = [&](Types::scalar t) -> Types::Vector3c {
                const Types::Vector3d y = cell.parametrization(1, t);
                return Helmholtz::V(k, point, y);
            };
            const auto CD = [&](Types::scalar t) -> Types::Vector3c {
                const Types::Vector3d y = cell.parametrization(1 - t, 1);
                return Helmholtz::V(k, point, y);
            };
            const auto DA = [&](Types::scalar t) -> Types::Vector3c {
                const Types::Vector3d y = cell.parametrization(0, 1 - t);
                return Helmholtz::V(k, point, y);
            };

            return DefiniteIntegrals::integrate<Quadrature>(AB, {0}, {1}) *
                   cell.integrationParameters.mul[0].transpose() +
                   DefiniteIntegrals::integrate<Quadrature>(BC, {0}, {1}) *
                   cell.integrationParameters.mul[1].transpose() +
                   DefiniteIntegrals::integrate<Quadrature>(CD, {0}, {1}) *
                   cell.integrationParameters.mul[2].transpose() +
                   DefiniteIntegrals::integrate<Quadrature>(DA, {0}, {1}) *
                   cell.integrationParameters.mul[3].transpose();
        }
    }

    /**
     * Расчет интергального оператора К1/{k^2} на одной ячейке
     * @tparam Quadrature (тип) точки квадратуры
     * @tparam cell_t тип ячейки (Cell, IndexedCell)
     * @param point точка счёта
     * @param j постоянное векторное поле на ячейке
     * @param cell ячейка
     * @param k волновое число
     * @return значение K1 в точке point деленое на k^2
     */
    template<typename Quadrature, typename cell_t>
    Types::Vector3c K1OverSingularCellDivided(const Mesh::Point &point, const Types::Vector3c &j, const cell_t &cell,
                                              const Types::scalar k) {
        return j * detail::K1OverSingularReducedAndDivided<Quadrature>(point, cell, k);
    }

    /** Полный расчет оператора K1 в точке point */
    template<typename Quadrature, typename cell_t>
    Types::Vector3c K1(const Mesh::Point &point, const Containers::vector<cell_t> &cells, const Types::scalar k) {
        Types::Vector3c result = Types::Vector3c::Zero();
        for (const auto &cell: cells) {
            result += K1OverSingularCellDivided<Quadrature>(point, cell.collPoint_.J_, cell, k);
        }
        return result * k * k;
    }

    template<typename Quadrature, typename cell_t>
    Types::Vector3c K0OverSingularCell(const Mesh::Point &point, const Types::Vector3c &j, const cell_t &cell, Types::scalar k) {
        return detail::K0TensorOverSingularCell<Quadrature>(point, cell, k) * j;
    }

    /** Полный расчет оператора K0 в точке point */
    template<typename Quadrature, typename cell_t>
    Types::Vector3c K0(const Mesh::Point &point, const Containers::vector<cell_t> &cells, const Types::scalar k) {
        Types::Vector3c result = Types::Vector3c::Zero();
        for (const auto &cell: cells) {
            result += K0OverSingularCell<Quadrature>(point, cell.collPoint_.J_, cell, k);
        }
        return result;
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_OPERATORS_HPP
