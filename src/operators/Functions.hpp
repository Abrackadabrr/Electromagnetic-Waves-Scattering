//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_FUNCTIONS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_FUNCTIONS_HPP

#include "types/Types.hpp"

namespace EMW::Helmholtz {

    // фундаментальное решение уравнения Гельмгольца
    Types::complex_d F(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y);

    // ограниченный кусок от фундаментального решения, вычли 1 / r
    Types::complex_d F_bounded_part(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y);

    // МИНУС градиент по x фундаментального решения уравнения Гельмгольца
    Types::Vector3c V(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y);

    // функция, которая считает подинтегральное выражение для интеграла диаграммы направленности
    Types::Vector3c
    sigmaKernel(Types::complex_d k, const Types::Vector3d &tau, const Types::Vector3d &point_on_surface, const Types::Vector3c &j_e, const Types::Vector3c& j_m);

    // функция для сглаживания сингулярности
    Types::scalar smoother(Types::scalar e, const Types::Vector3d &x, const Types::Vector3d &y);

    Types::Vector3c
    reducedK_kernel(Types::complex_d k, const Types::Vector3d &x, const Types::Vector3d &y, const Types::Vector3c &j);
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_FUNCTIONS_HPP
