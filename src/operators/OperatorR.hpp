//
// Created by evgen on 25.11.2024.
//

#ifndef OPERATORR_HPP
#define OPERATORR_HPP

#include "Functions.hpp"
#include "math/MathConstants.hpp"
#include "math/fields/SurfaceVectorField.hpp"
#include "math/integration/Quadrature.hpp"
#include "mesh/MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::OperatorR {

// в этом namespace будут функции, позволяющие наиболее просто
// вычислить матрицу дискретизованного оператора R
namespace detail::forMatrix {
/**
 * Нахождение общего интеграла для четырёх блоков в матрице
 */
template <typename Quadrature, typename cell_t>
Types::Vector3c commonIntegralPart(const cell_t &cell_j, const Mesh::point_t &x_i, Types::complex_d k) {
    // замена переменных в интергале позволяет интегрироваться по квадрату [0, 1]^2
    // phi -- это есть функция, переписанная в новых переменных и помноженная на якобиан отображения
    // отображение в нашем случае билинейное
    const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::Vector3c {
        const Types::Vector3d y = cell_j.parametrization(p, q);
        const Types::scalar mul = cell_j.multiplier(p, q);
        return Helmholtz::V(k, x_i, y) * mul;
    };
    return DefiniteIntegrals::integrate<Quadrature>(phi, {0, 0}, {1., 1.});
}

}

} // namespace EMW::Operators

#endif // OPERATORR_HPP
