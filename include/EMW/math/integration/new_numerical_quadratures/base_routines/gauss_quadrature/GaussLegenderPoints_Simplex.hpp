//
// Created by evgen on 22.10.2025.
//

#ifndef SPECIALGAUSSLEGENERPOINTS_HPP
#define SPECIALGAUSSLEGENERPOINTS_HPP

#include "Types.hpp"
#include "math/integrators/numerical/base_routines/QuadratureUtils.hpp"

namespace numedie::math::integration::base_routines::GaussLegender {

template <unsigned int dimension_, unsigned int order_> struct GaussianPoints {
  static_assert("We have no quadrature of such dimension and order");
};

// --- 1D --- //

template <> struct GaussianPoints<1, 1> : quadrature_utils::BasePoints<1, 1> {
  using Base = quadrature_utils::BasePoints<1, 1>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>,
                                     Base::order>
      nodes = {quadrature_utils::Node<Base::dim>{{1. / 2}, 1}};
};

template <> struct GaussianPoints<1, 3> : quadrature_utils::BasePoints<1, 3> {
  using Base = quadrature_utils::BasePoints<1, 3>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>,
                                     Base::order - 1>
      nodes = {quadrature_utils::Node<Base::dim>{
                   {1. / 2 * (1 - 1. / std::sqrt(3.))}, 1. / 2},
               quadrature_utils::Node<Base::dim>{
                   {1. / 2 * (1 - 1. / std::sqrt(3.))}, 1. / 2}};
};

template <> struct GaussianPoints<1, 5> : quadrature_utils::BasePoints<1, 5> {
  using Base = quadrature_utils::BasePoints<1, 5>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>,
                                     Base::order - 2>
      nodes = {
          quadrature_utils::Node<Base::dim>{
              {1. / 2 - (1. / 2) * std::sqrt(3. / 5)}, 5. / 18},
          quadrature_utils::Node<Base::dim>{
              {1. / 2 + (1. / 2) * std::sqrt(3. / 5)}, 5. / 18},
          quadrature_utils::Node<Base::dim>{{1. / 2}, 4. / 9},
      };
};

template <> struct GaussianPoints<1, 7> : quadrature_utils::BasePoints<1, 7> {
  using Base = quadrature_utils::BasePoints<1, 7>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>,
                                     Base::order - 3>
      nodes = {
          quadrature_utils::Node<Base::dim>{
              {1. / 2 -
               (1. / 2) * std::sqrt(3. / 7 - (2. / 7) * std::sqrt(6. / 5))},
              1. / 4 + std::sqrt(30) / 72},
          quadrature_utils::Node<Base::dim>{
              {1. / 2 +
               (1. / 2) * std::sqrt(3. / 7 - (2. / 7) * std::sqrt(6. / 5))},
              1. / 4 + std::sqrt(30) / 72},
          quadrature_utils::Node<Base::dim>{
              {1. / 2 -
               (1. / 2) * std::sqrt(3. / 7 + (2. / 7) * std::sqrt(6. / 5))},
              1. / 4 - std::sqrt(30) / 72},
          quadrature_utils::Node<Base::dim>{
              {1. / 2 +
               (1. / 2) * std::sqrt(3. / 7 + (2. / 7) * std::sqrt(6. / 5))},
              1. / 4 - std::sqrt(30) / 72},
      };
};

template <> struct GaussianPoints<1, 9> : quadrature_utils::BasePoints<1, 9> {
  using Base = quadrature_utils::BasePoints<1, 9>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>,
                                     Base::order - 4>
      nodes = {
          quadrature_utils::Node<Base::dim>{{1. / 2}, 64. / 225},
          quadrature_utils::Node<Base::dim>{
              {1. / 2 - (1. / 6) * std::sqrt(5 - 2. * std::sqrt(10. / 7))},
              (322 + 13 * std::sqrt(70)) / 1800},
          quadrature_utils::Node<Base::dim>{
              {1. / 2 + (1. / 6) * std::sqrt(5 - 2. * std::sqrt(10. / 7))},
              (322 + 13 * std::sqrt(70)) / 1800},
          quadrature_utils::Node<Base::dim>{
              {1. / 2 - (1. / 6) * std::sqrt(5 + 2. * std::sqrt(10. / 7))},
              (322 - 13 * std::sqrt(70)) / 1800},
          quadrature_utils::Node<Base::dim>{
              {1. / 2 + (1. / 6) * std::sqrt(5 + 2. * std::sqrt(10. / 7))},
              (322 - 13 * std::sqrt(70)) / 1800},
      };
};

// --------- порядок 11 (n = 6 точек) ----------
template <> struct GaussianPoints<1, 11> : quadrature_utils::BasePoints<1, 11> {
  using Base = quadrature_utils::BasePoints<1, 11>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 6>
      nodes = {
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) -
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.932469514203152027812301554493994)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.171324492379170345130443202740000)},
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) -
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.661209386466264513661399595019905)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.360761573048138607569833513837716)},
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) -
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.238619186083196908630501721680711)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.467913934572691047389870343989551)},
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) +
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.238619186083196908630501721680711)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.467913934572691047389870343989551)},
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) +
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.661209386466264513661399595019905)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.360761573048138607569833513837716)},
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) +
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.932469514203152027812301554493994)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.171324492379170345130443202740000)},
      };
};

// --------- порядок 13 (n = 7 точек) ----------
template <> struct GaussianPoints<1, 13> : quadrature_utils::BasePoints<1, 13> {
  using Base = quadrature_utils::BasePoints<1, 13>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 7>
      nodes = {
          // t = -0.949107912342758524526189684047851...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) -
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.949107912342758524526189684047851)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.129484966168869693270611432679000)},
          // t = -0.741531185599394439863864773280788...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) -
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.741531185599394439863864773280788)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.279705391489276667901467771423000)},
          // t = -0.405845151377397166906606412076961...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) -
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.405845151377397166906606412076961)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.381830050505118944950369775489000)},
          // t = 0  (центр)
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)},
              // центральный вес на [-1,1] = 512/1225; на [0,1] делим на 2 =>
              // 256/1225
              static_cast<Types::scalar>(256) /
                  static_cast<Types::scalar>(1225)},
          // t = +0.405845151377397166906606412076961...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) +
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.405845151377397166906606412076961)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.381830050505118944950369775489000)},
          // t = +0.741531185599394439863864773280788...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) +
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.741531185599394439863864773280788)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.279705391489276667901467771423000)},
          // t = +0.949107912342758524526189684047851...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) +
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.949107912342758524526189684047851)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.129484966168869693270611432679000)},
      };
};

// --------- порядок 15 (n = 8 точек) ----------
template <> struct GaussianPoints<1, 15> : quadrature_utils::BasePoints<1, 15> {
  using Base = quadrature_utils::BasePoints<1, 15>;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>, 8>
      nodes = {
          // t = -0.960289856497536231683560868569473...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) -
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.960289856497536231683560868569473)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.101228536290376259152531354309000)},
          // t = -0.796666477413626739591553936475830...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) -
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.796666477413626739591553936475830)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.222381034453374470544355994426000)},
          // t = -0.525532409916328985817739049189246...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) -
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.525532409916328985817739049189246)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.313706645877887287337962201987000)},
          // t = -0.183434642495649804939476142360764...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) -
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.183434642495649804939476142360764)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.362683783378361982965150449277000)},
          // t = +0.183434642495649804939476142360764...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) +
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.183434642495649804939476142360764)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.362683783378361982965150449277000)},
          // t = +0.525532409916328985817739049189246...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) +
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.525532409916328985817739049189246)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.313706645877887287337962201987000)},
          // t = +0.796666477413626739591553936475830...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) +
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.796666477413626739591553936475830)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.222381034453374470544355994426000)},
          // t = +0.960289856497536231683560868569473...
          quadrature_utils::Node<Base::dim>{
              {static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2) +
               (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                   static_cast<Types::scalar>(
                       0.960289856497536231683560868569473)},
              (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) *
                  static_cast<Types::scalar>(
                      0.101228536290376259152531354309000)},
      };
};

template <> struct GaussianPoints<1, 17> : quadrature_utils::BasePoints<1, 17> {
  using Base = quadrature_utils::BasePoints<1, 17>;
  static constexpr Types::scalar k2 = 0.968160239507626;
  static constexpr Types::scalar k3 = 0.836031107326636;
  static constexpr Types::scalar k4 = 0.613371432700590;
  static constexpr Types::scalar k5 = 0.324253423403809;
  static constexpr Containers::array<quadrature_utils::Node<Base::dim>,
                                     Base::order - 8>
      nodes = {
          quadrature_utils::Node<Base::dim>{{1. / 2}, 0.165119677500630},

          quadrature_utils::Node<Base::dim>{{1. / 2 - k2 / 2},
                                            0.040637194180787},
          quadrature_utils::Node<Base::dim>{{1. / 2 + k2 / 2},
                                            0.040637194180787},

          quadrature_utils::Node<Base::dim>{{1. / 2 - k3 / 2},
                                            0.090324080347429},
          quadrature_utils::Node<Base::dim>{{1. / 2 + k3 / 2},
                                            0.090324080347429},

          quadrature_utils::Node<Base::dim>{{1. / 2 - k4 / 2},
                                            0.130305348201468},
          quadrature_utils::Node<Base::dim>{{1. / 2 + k4 / 2},
                                            0.130305348201468},

          quadrature_utils::Node<Base::dim>{{1. / 2 - k5 / 2},
                                            0.156173538520002},
          quadrature_utils::Node<Base::dim>{{1. / 2 + k5 / 2},
                                            0.156173538520002},
      };
};

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
