//
// Created by evgen on 22.10.2025.
//

#ifndef TRUANGLEMIDDLERECTANGLES_HPP
#define TRUANGLEMIDDLERECTANGLES_HPP

#include "Types.hpp"
#include "math/integrators/numerical/base_routines/QuadratureUtils.hpp"


namespace numedie::math::integration::base_routines::NewtonCotess {
template <unsigned int dimention_, unsigned int division_size>
struct RectangularPoints {
  static_assert("We have no quadrature of such dimension and order");
};

// --- 1D --- //

template <>
struct RectangularPoints<1, 1> : quadrature_utils::BasePoints<1, 1> {
  static constexpr Containers::array<quadrature_utils::Node<dim>, order> nodes =
      {quadrature_utils::Node<dim>{{1. / 2}, 1}};
};

template <>
struct RectangularPoints<1, 2> : quadrature_utils::BasePoints<1, 2> {
  static constexpr Containers::array<quadrature_utils::Node<dim>, order> nodes =
      {quadrature_utils::Node<dim>{{0.25}, 1. / 2.},
       quadrature_utils::Node<dim>{{0.75}, 1. / 2.}};
};

template <>
struct RectangularPoints<1, 3> : quadrature_utils::BasePoints<1, 3> {
  static constexpr Containers::array<quadrature_utils::Node<dim>, order> nodes =
      {quadrature_utils::Node<dim>{{1. / 6.}, 1. / 3.},
       quadrature_utils::Node<dim>{{3. / 6.}, 1. / 3.},
       quadrature_utils::Node<dim>{{5. / 6.}, 1. / 3.}};
};


// --- 2D (flat triangle) --- //

template <>
struct RectangularPoints<2, 1> : quadrature_utils::BasePoints<2, 1> {
  static constexpr Containers::array<quadrature_utils::Node<dim>, order> nodes =
      {quadrature_utils::Node<dim>{{1. / 3, 1./ 3}, 1.}};
};

template <>
struct RectangularPoints<2, 4> : quadrature_utils::BasePoints<2, 4> {
  static constexpr Containers::array<quadrature_utils::Node<dim>, order> nodes =
      {quadrature_utils::Node<dim>{{1. / 6, 1./ 6}, 1. / 4.},
      quadrature_utils::Node<dim>{{1. / 3, 1./ 3}, 1. / 4.},
      quadrature_utils::Node<dim>{{1. / 6, 2./ 3}, 1. / 4.},
      quadrature_utils::Node<dim>{{2. / 3, 1./ 6}, 1. / 4.}};
};


}
#endif //TRUANGLEMIDDLERECTANGLES_HPP
