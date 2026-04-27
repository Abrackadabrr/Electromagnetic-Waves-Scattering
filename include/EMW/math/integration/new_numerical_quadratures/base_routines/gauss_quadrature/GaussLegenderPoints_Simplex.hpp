//
// Created by evgen on 22.10.2025.
//

#ifndef SPECIALGAUSSLEGENERPOINTS_HPP
#define SPECIALGAUSSLEGENERPOINTS_HPP

#include "EMW/Types.hpp"

namespace numedie::math::integration::base_routines::GaussLegender {



// --- 2D --- //

template <> struct GaussianPoints<2, 1> : quadrature_utils::BasePoints<2, 1> {
  using Base = quadrature_utils::BasePoints<2, 1>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>,
                                     Base::order>
      nodes = {
          quadrature_utils::Node<Base::dim>{{1. / 3, 1. / 3}, 1.},
      };
};

template <> struct GaussianPoints<2, 2> : quadrature_utils::BasePoints<2, 2> {
  using Base = quadrature_utils::BasePoints<2, 2>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 3>
      nodes = {
          quadrature_utils::Node<Base::dim>{{1. / 6, 1. / 6}, 1. / 3},
          quadrature_utils::Node<Base::dim>{{1. / 6, 2. / 3}, 1. / 3},
          quadrature_utils::Node<Base::dim>{{2. / 3, 1. / 6}, 1. / 3},
      };
};

template <> struct GaussianPoints<2, 3> : quadrature_utils::BasePoints<2, 3> {
  using Base = quadrature_utils::BasePoints<2, 3>;
  static constexpr Types::scalar l1 = 0.2319333685530305;
  static constexpr Types::scalar l2 = 0.6590276223740922;
  static constexpr Types::scalar l3 = 1 - l2 - l1;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 6>
      nodes = {
          quadrature_utils::Node<Base::dim>{{l1, l2}, 1. / 6},
          quadrature_utils::Node<Base::dim>{{l2, l1}, 1. / 6},
          quadrature_utils::Node<Base::dim>{{l1, l3}, 1. / 6},
          quadrature_utils::Node<Base::dim>{{l3, l1}, 1. / 6},
          quadrature_utils::Node<Base::dim>{{l3, l2}, 1. / 6},
          quadrature_utils::Node<Base::dim>{{l2, l3}, 1. / 6},
      };
};

template <> struct GaussianPoints<2, 4> : quadrature_utils::BasePoints<2, 4> {
  using Base = quadrature_utils::BasePoints<2, 4>;
  static constexpr Types::scalar l1 = 0.44594849091596489;
  static constexpr Types::scalar w1 = 0.22338158967801147;
  static constexpr Types::scalar l2 = 0.091576213509770854;
  static constexpr Types::scalar w2 = 0.10995174365532188;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 6>
      nodes = {
          quadrature_utils::Node<Base::dim>{{l1, l1}, w1},
          quadrature_utils::Node<Base::dim>{{l1, 1 - 2 * l1}, w1},
          quadrature_utils::Node<Base::dim>{{1 - 2 * l1, l1}, w1},
          quadrature_utils::Node<Base::dim>{{l2, l2}, w2},
          quadrature_utils::Node<Base::dim>{{l2, 1 - 2 * l2}, w2},
          quadrature_utils::Node<Base::dim>{{1 - 2 * l2, l2}, w2},
      };
};

template <> struct GaussianPoints<2, 5> : quadrature_utils::BasePoints<2, 5> {
  using Base = quadrature_utils::BasePoints<2, 5>;
  static constexpr Types::scalar l1 = (6 + std::sqrt(15)) / 21;
  static constexpr Types::scalar w1 = (155 + std::sqrt(15)) / 1200;
  static constexpr Types::scalar w2 = (155 - std::sqrt(15)) / 1200;
  static constexpr Types::scalar l2 = (6 - std::sqrt(15)) / 21;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 7>
      nodes = {
          quadrature_utils::Node<Base::dim>{{1. / 3, 1. / 3}, 0.225},

          quadrature_utils::Node<Base::dim>{{l1, l1}, w1},
          quadrature_utils::Node<Base::dim>{{l1, 1 - 2 * l1}, w1},
          quadrature_utils::Node<Base::dim>{{1 - 2 * l1, l1}, w1},

          quadrature_utils::Node<Base::dim>{{l2, l2}, w2},
          quadrature_utils::Node<Base::dim>{{l2, 1 - 2 * l2}, w2},
          quadrature_utils::Node<Base::dim>{{1 - 2 * l2, l2}, w2},
      };
};

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
