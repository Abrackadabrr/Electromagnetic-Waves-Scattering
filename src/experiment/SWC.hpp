//
// Created by evgen on 09.01.2025.
//

#ifndef SWC_HPP
#define SWC_HPP

#include "types/Types.hpp"
#include "math/fields/SurfaceVectorField.hpp"

namespace EMW::EngineeringCoefficients {

/**
 * Расчитывает коэффициент стоячей волны для поля
 * @param field собственно поле
 * @return
 */
template<typename field1, typename field2>
Containers::vector<Types::scalar> SWC(const field1& direct_f, const field2& inverse_f, const Containers::vector<Mesh::point_t>& points) {
    Containers::vector<Types::scalar> swc;
    for (auto point : points) {
        const auto direct_v = direct_f(point).norm();
        const auto inverse_v = inverse_f(point).norm();
        const Types::scalar value = (direct_v + inverse_v) / (direct_v - inverse_v);
        swc.emplace_back(value);
    }
    return swc;
}

}


#endif //SWC_HPP
