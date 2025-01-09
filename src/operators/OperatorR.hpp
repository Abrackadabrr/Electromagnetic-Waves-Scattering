//
// Created by evgen on 25.11.2024.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_OPERATORR_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_OPERATORR_HPP

#include "Functions.hpp"
#include "math/integration/Quadrature.hpp"
#include "mesh/MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::OperatorR {

// в этом namespace будут функции, позволяющие наиболее просто
// вычислить матрицу дискретизованного оператора R
namespace detail {
namespace forMatrix {
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

/** Численный расчет оператора R в точке point */
template <typename Quadrature>
Types::Vector3c R(const Mesh::point_t &point, const Types::complex_d k, const Math::SurfaceVectorField &field) {
  const auto &cells = field.getManifold().getCells();
  const auto &f = field.getField();
  Types::Vector3c result = Types::Vector3c::Zero();

  for (int i = 0; i != cells.size(); i++) {
    result += forMatrix::commonIntegralPart<Quadrature>(cells[i], point, k).cross(f[i]);
  }

  return result;
}

}



} // namespace EMW::Operators

#endif // ELECTROMAGNETIC_WAVES_SCATTERING_
