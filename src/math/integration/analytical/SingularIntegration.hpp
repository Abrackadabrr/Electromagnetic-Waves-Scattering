//
// Created by evgen on 23.08.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_SINGULARINTEGRATION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_SINGULARINTEGRATION_HPP

#include "mesh/MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::Math::AnalyticalIntegration{
    /**
     * Интегрирование функции 1/|r - r'| по поверхностной ячейке (треугольник или четырёхугольник)
     * @param point - r
     * @param cell - ячейка (r')
     * @return значение интеграла по этой ячейке
     */
    Types::scalar integrate_1_div_r(const Mesh::point_t &point, const Mesh::IndexedCell &cell);
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_SINGULARINTEGRATION_HPP
