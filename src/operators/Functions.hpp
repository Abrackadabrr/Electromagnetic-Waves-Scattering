//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_FUNCTIONS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_FUNCTIONS_HPP

#include "types/Types.hpp"

namespace EMW::Helmholtz
{
    // фундаментальное решение уравнения Гельмгольца
    Types::complex_d F(Types::complex_d k, const Types::Vector3d& x, const Types::Vector3d& y);

    // ограниченный кусок от фундаментального решения, вычли 1 / (4 pi r)
    Types::complex_d F_bounded_part(Types::complex_d k, const Types::Vector3d& x, const Types::Vector3d& y);

    // МИНУС градиент по x фундаментального решения уравнения Гельмгольца
    Types::Vector3c V(Types::complex_d k, const Types::Vector3d& x, const Types::Vector3d& y);

    // Подинтегральное выражение для диаграммы направленности
    Types::Vector3c
    sigmaKernel(Types::complex_d k, const Types::Vector3d& tau, const Types::Vector3d& point_on_surface,
                const Types::Vector3c& j_e, const Types::Vector3c& j_m, Types::complex_d epsilon = {1., 0.});

    // функция для сглаживания сингулярности
    Types::scalar smoother(Types::scalar e, const Types::Vector3d& x, const Types::Vector3d& y);

    /**
     * @brief Ядро интегрального оператора при внесении производных под интеграл
     *
     * @note Е. В. Захаров, Г. В. Рыжаков, А. В. Сетуха,
     * "ЧИСЛЕННОЕ РЕШЕНИЕ ТРЕХМЕРНЫХ ЗАДАЧ ДИФРАКЦИИ ЭЛЕКТРОМАГНИТНЫХ ВОЛН
     * НА СИСТЕМЕ ИДЕАЛЬНО ПРОВОДЯЩИХ ПОВЕРХНОСТЕЙ МЕТОДОМ ГИПЕРСИНГУЛЯРНЫХ
     * ИНТЕГРАЛЬНЫХ УРАВНЕНИЙ", "ДИФФЕРЕНЦИАЛЬНЫЕ УРАВНЕНИЯ", 2014, том 50, № 9, с. 1253–1263
     *
     * @param r = x - y
     * @param j ток на ячейке в точке (кусочно-постоянная функция на ячейке)
     *
     * @warning Формула верна для кусочно постоянной функции j на ячейке
     */
    [[nodiscard]] Types::Vector3c far_zone_integral_kernel(Types::complex_d k, const Types::point_t &r, const Types::Vector3c &j);
}

namespace EMW::Laplace
{
    // фундаментальное решение уравнения Лапласа
    Types::scalar F(const Types::Vector3d& x, const Types::Vector3d& y);

    // градиент фундаментального решения уравнения Лапласа
    Types::Vector3d gradF(const Types::Vector3d& x, const Types::Vector3d& y);
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_FUNCTIONS_HPP
