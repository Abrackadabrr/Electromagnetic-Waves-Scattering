//
// Created by evgen on 22.10.2025.
//

#ifndef SPECIALGAUSSLEGENERPOINTS_HPP
#define SPECIALGAUSSLEGENERPOINTS_HPP

#include "EMW/Types.hpp"

namespace numedie::math::integration::base_routines::GaussLegender {



// --- 2D --- //


// --- 3D --- //

template <> struct GaussianPoints<3, 1> : quadrature_utils::BasePoints<3, 1> {
  using Base = quadrature_utils::BasePoints<3, 1>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 1>
      nodes = {quadrature_utils::Node<Base::dim>{{0.25, 0.25, 0.25}, 1}};
};

template <> struct GaussianPoints<3, 2> : quadrature_utils::BasePoints<3, 2> {
  using Base = quadrature_utils::BasePoints<3, 2>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 4>
      nodes = {
          quadrature_utils::Node<Base::dim>{
              {0.1381966011250105, 0.1381966011250105, 0.1381966011250105},
              0.25},
          quadrature_utils::Node<Base::dim>{
              {0.5854101966249685, 0.1381966011250105, 0.1381966011250105},
              0.25},
          quadrature_utils::Node<Base::dim>{
              {0.1381966011250105, 0.5854101966249685, 0.1381966011250105},
              0.25},
          quadrature_utils::Node<Base::dim>{
              {0.1381966011250105, 0.1381966011250105, 0.5854101966249685},
              0.25},
      };
};

template <> struct GaussianPoints<3, 3> : quadrature_utils::BasePoints<3, 3> {
  using Base = quadrature_utils::BasePoints<3, 3>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 8>
      nodes = {
          quadrature_utils::Node<Base::dim>{
              {0.3286811466653490, 0.3286811466653490, 0.3286811466653490},
              0.1274913115575064},
          quadrature_utils::Node<Base::dim>{
              {0.013956560003953067, 0.3286811466653490, 0.3286811466653490},
              0.1274913115575064},
          quadrature_utils::Node<Base::dim>{
              {0.3286811466653490, 0.013956560003953067, 0.3286811466653490},
              0.1274913115575064},
          quadrature_utils::Node<Base::dim>{
              {0.3286811466653490, 0.3286811466653490, 0.013956560003953067},
              0.1274913115575064},

          quadrature_utils::Node<Base::dim>{
              {0.1119207275092915, 0.1119207275092915, 0.1119207275092915},
              0.1225086884424935},
          quadrature_utils::Node<Base::dim>{
              {0.6642378174721255, 0.1119207275092915, 0.1119207275092915},
              0.1225086884424935},
          quadrature_utils::Node<Base::dim>{
              {0.1119207275092915, 0.6642378174721255, 0.1119207275092915},
              0.1225086884424935},
          quadrature_utils::Node<Base::dim>{
              {0.1119207275092915, 0.1119207275092915, 0.6642378174721255},
              0.1225086884424935}};
};

template <> struct GaussianPoints<3, 5> : quadrature_utils::BasePoints<3, 5> {
  using Base = quadrature_utils::BasePoints<3, 5>;
  static constexpr Types::scalar l1 = 0.3108859192633006;
  static constexpr Types::scalar l2 = 0.0927352503108912;
  static constexpr Types::scalar l3 = 0.04550370412564965;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 14>
      nodes = {
          quadrature_utils::Node<Base::dim>{{l1, l1, l1}, 0.1274913115575064},
          quadrature_utils::Node<Base::dim>{{1 - 3 * l1, l1, l1},
                                            0.1274913115575064},
          quadrature_utils::Node<Base::dim>{{l1, 1 - 3 * l1, l1},
                                            0.1274913115575064},
          quadrature_utils::Node<Base::dim>{{l1, l1, 1 - 3 * l1},
                                            0.1274913115575064},

          quadrature_utils::Node<Base::dim>{{l2, l2, l2}, 0.07349304311636194},
          quadrature_utils::Node<Base::dim>{{1 - 3 * l2, l2, l2},
                                            0.07349304311636194},
          quadrature_utils::Node<Base::dim>{{l2, 1 - 3 * l2, l2},
                                            0.07349304311636194},
          quadrature_utils::Node<Base::dim>{{l2, l2, 1 - 3 * l2},
                                            0.07349304311636194},

          quadrature_utils::Node<Base::dim>{{l3, l3, 1. / 2 - l3},
                                            0.04254602077708146},
          quadrature_utils::Node<Base::dim>{{l3, 1. / 2 - l3, l3},
                                            0.04254602077708146},
          quadrature_utils::Node<Base::dim>{{1. / 2 - l3, l3, l3},
                                            0.04254602077708146},
          quadrature_utils::Node<Base::dim>{{l3, 1. / 2 - l3, 1. / 2 - l3},
                                            0.04254602077708146},
          quadrature_utils::Node<Base::dim>{{1. / 2 - l3, l3, 1. / 2 - l3},
                                            0.04254602077708146},
          quadrature_utils::Node<Base::dim>{{1. / 2 - l3, 1. / 2 - l3, l3},
                                            0.04254602077708146},
      };
};

} // namespace numedie::math::integration::base_routines::GaussLegender

#endif // SPECIALGAUSSLEGENERPOINTS_HPP
