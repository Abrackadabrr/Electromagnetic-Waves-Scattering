//
// Created by evgen on 13.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_OPERATORS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_OPERATORS_HPP

#include <types/Types.hpp>
#include "mesh/MeshTypes.hpp"
#include "Functions.hpp"
#include "math/fields/SurfaceField.hpp"
#include "math/integration/Quadrature.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/MathConstants.hpp"

namespace EMW::Operators {
    namespace detail {
        /**
         * Расчет интергального оператора К1/{k^2} на одной ячейке без умножения на вектор
         * Значение интеграла считается численно при помощи квадратуры, подаваемой в шаблонный параметр
         * Происзодит преобразование ячейки в стандрный квадрат, на котором проивходит расчет с помощью
         * квадратурной формулы.
         * Данный метод не подходит для расчета self-terms в матрице, так как не позволяет с достаточной точностью
         * рассчитывать несобственные интегралы.
         * @tparam Quadrature (тип) точки квадратуры
         * @tparam cell_t тип ячейки (Cell, IndexedCell)
         * @param point точка для расчета
         * @param cell ячейка, на которой считаем
         * @param k волновое число
         * @return интеграл от F(x - y) dy по ячейке cell
         */
        template<typename Quadrature, typename cell_t>
        Types::complex_d
        K1OverSingularCellReducedAndDivided(const Mesh::point_t &point, const cell_t &cell, const Types::scalar k) {
            const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
                const Types::Vector3d y = cell.parametrization(p, q);
                const Types::scalar mul = cell.multiplier(p, q);
                return Helmholtz::F(k, point, y) * mul;
            };
            return DefiniteIntegrals::integrate<Quadrature>(phi, {0, 0}, {1., 1.});
        }

        /**
         * Расчет интергального оператора К1/{k^2} на одной ячейке без умножения на вектор
         * Расчет интеграла от функции Грина оператора Гельмгольца с регуляризацией
         * Квадратура используется для интегрирования непрерывно-дифференцируемого слагаемого, интегрирование просходит
         * по определению поверхностного интеграла с параметризацией поверхности
         *
         * @tparam Quadrature (тип) точки квадратуры
         * @tparam cell_t тип ячейки (Cell, IndexedCell)
         * @param point точка для расчета
         * @param cell ячейка, на которой считаем
         * @param k волновое число
         * @return интеграл от F(x - y) dy по ячейке cell
         */
        template<typename Quadrature, typename cell_t>
        Types::complex_d
        K1OverSingularCellRnDWithSingularityExtraction(const Mesh::point_t &point, const cell_t &cell,
                                                       const Types::scalar k) {
            const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
                const Types::Vector3d y = cell.parametrization(p, q);
                const Types::scalar mul = cell.multiplier(p, q);
                return Helmholtz::F_bounded_part(k, point, y) * mul;
            };
            return DefiniteIntegrals::integrate<Quadrature>(phi, {0, 0}, {1., 1.}) +
                   Math::Constants::inverse_4PI<Types::scalar>() *
                   Math::AnalyticalIntegration::integrate_1_div_r(point, cell);
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
        K0TensorOverSingularCell(const Mesh::point_t &point, const cell_t &cell, Types::scalar k) {
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
     * Численный расчет интергального оператора К1/{k^2} на одной ячейке
     * @tparam Quadrature (тип) точки квадратуры
     * @tparam cell_t тип ячейки (Cell, IndexedCell)
     * @param point точка счёта
     * @param j постоянное векторное поле на ячейке
     * @param cell ячейка
     * @param k волновое число
     * @return значение K1 в точке point деленое на k^2
     */
    template<typename Quadrature, typename cell_t>
    Types::Vector3c K1OverSingularCellDivided(const Mesh::point_t &point, const Types::Vector3c &j, const cell_t &cell,
                                              const Types::scalar k) {
        return j * detail::K1OverSingularCellReducedAndDivided<Quadrature>(point, cell, k);
    }

    /**
     * Полуаналитический расчет интергального оператора К1/{k^2} на одной ячейке
     * @tparam Quadrature (тип) точки квадратуры
     * @tparam cell_t тип ячейки (Cell, IndexedCell)
     * @param point точка счёта
     * @param j постоянное векторное поле на ячейке
     * @param cell ячейка
     * @param k волновое число
     * @return значение K1 в точке point деленое на k^2
     */
    template<typename Quadrature, typename cell_t>
    Types::Vector3c K1OverSingularCellDividedSingularityExtraction(const Mesh::point_t &point, const Types::Vector3c &j, const cell_t &cell,
                                              const Types::scalar k) {
        return j * detail::K1OverSingularCellRnDWithSingularityExtraction<Quadrature>(point, cell, k);
    }

    /**
     * Численный интергального оператора К1 на одной ячейке
     * @tparam Quadrature (тип) точки квадратуры
     * @tparam cell_t тип ячейки (Cell, IndexedCell)
     * @param point точка счёта
     * @param j постоянное векторное поле на ячейке
     * @param cell ячейка
     * @param k волновое число
     * @return значение K1 в точке point деленое на k^2
     */
    template<typename Quadrature, typename cell_t>
    Types::Vector3c K1OverSingularCell(const Mesh::point_t &point, const Types::Vector3c &j, const cell_t &cell,
                                       const Types::scalar k) {
        return k * k * j * detail::K1OverSingularCellReducedAndDivided<Quadrature>(point, cell, k);
    }

    /**
     * Полуаналитический интергального оператора К1 на одной ячейке
     * @tparam Quadrature (тип) точки квадратуры
     * @tparam cell_t тип ячейки (Cell, IndexedCell)
     * @param point точка счёта
     * @param j постоянное векторное поле на ячейке
     * @param cell ячейка
     * @param k волновое число
     * @return значение K1 в точке point деленое на k^2
     */
    template<typename Quadrature, typename cell_t>
    Types::Vector3c K1OverSingularCellSingularityExtraction(const Mesh::point_t &point, const Types::Vector3c &j, const cell_t &cell,
                                       const Types::scalar k) {
        return k * k * j * detail::K1OverSingularCellRnDWithSingularityExtraction<Quadrature>(point, cell, k);
    }

    template<typename Quadrature, typename cell_t>
    Types::Vector3c
    K0OverSingularCell(const Mesh::point_t &point, const Types::Vector3c &j, const cell_t &cell, Types::scalar k) {
        return detail::K0TensorOverSingularCell<Quadrature>(point, cell, k) * j;
    }


    /** Численный расчет оператора K1 в точке point */
    template<typename Quadrature>
    Types::Vector3c K1(const Mesh::point_t &point, const Types::scalar k, const Math::SurfaceField &field) {
        const auto &cells = field.getManifold().getCells();
        const auto &f = field.getField();
        Types::Vector3c result = Types::Vector3c::Zero();
        for (int i = 0; i != cells.size(); i++) {
            result += K1OverSingularCellDivided<Quadrature>(point, f[i], cells[i], k);
        }
        return result * k * k;
    }

    /** Полуаналитический расчет оператора K1 в точке point */
    template<typename Quadrature>
    Types::Vector3c K1_singularityExtraction(const Mesh::point_t &point, const Types::scalar k, const Math::SurfaceField &field) {
        const auto &cells = field.getManifold().getCells();
        const auto &f = field.getField();
        Types::Vector3c result = Types::Vector3c::Zero();
        for (int i = 0; i != cells.size(); i++) {
            result += K1OverSingularCellDividedSingularityExtraction<Quadrature>(point, f[i], cells[i], k);
        }
        return result * k * k;
    }

    /** Полный расчет оператора K0 в точке point */
    template<typename Quadrature>
    Types::Vector3c K0(const Mesh::point_t &point, const Types::scalar k, const Math::SurfaceField &field) {
        const auto &cells = field.getManifold().getCells();
        const auto &f = field.getField();
        Types::Vector3c result = Types::Vector3c::Zero();
        for (int i = 0; i != cells.size(); i++) {
            result += K0OverSingularCell<Quadrature>(point, f[i], cells[i], k);
        }
        return result;
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_OPERATORS_HPP
