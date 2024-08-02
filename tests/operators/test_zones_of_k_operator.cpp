//
// Created by evgen on 14.02.24.
//


#include "gtest/gtest.h"

#include "operators/Operators.hpp"
#include "operators/Functions.hpp"
#include "types/Types.hpp"
#include "mesh/MeshTypes.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "integration/newton_cotess/Rectangular.hpp"

using namespace EMW;

template<Types::index N, Types::index M>
using GL2 = DefiniteIntegrals::GaussLegendre::Quadrature<N, M>;

template<Types::index N>
using GL = DefiniteIntegrals::GaussLegendre::Quadrature<N>;

template<Types::index N, Types::index M>
using R2 = DefiniteIntegrals::NewtonCotess::Quadrature<N, M>;

template<Types::index N>
using R = DefiniteIntegrals::NewtonCotess::Quadrature<N>;


class K1_TESTS : public testing::Test {
protected:
    Types::complex_d k{10, 0};
};

TEST_F(K1_TESTS, TEST_MIDDLE_ZONE) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::vector<Mesh::point_t> points = {
            Mesh::point_t{-0.05, -0.05, 0},
            Mesh::point_t{0.05, -0.05, 0},
            Mesh::point_t{0.05, 0.05, 0},
            Mesh::point_t{-0.05, 0.05, 0}
    };
    Mesh::IndexedCell cell{{0, 1, 2, 3}, points};
    cell.collPoint_.J_ = Types::Vector3c{Types::complex_d{0, 0}, Types::complex_d{1, 0}, Types::complex_d{0, 0}};
    // объявляем точку, в которой интегрируем
    const Mesh::point_t point{100, 0, 0};

    // точный расчет интеграла по ячейке
    const Types::Vector3c preciseResult = Operators::K0OverSingularCell<GL<8>>(point, cell.collPoint_.J_, cell, k) +
                                          k * k *
                                          Operators::K1OverSingularCellDivided<GL2<8, 8>>(point, cell.collPoint_.J_,
                                                                                          cell, k);

    std::cout << preciseResult << std::endl;

    // расчет с упрощенным ядром
    Types::Vector3c reducedResult =
            Helmholtz::reducedK_kernel(k, point, cell.collPoint_.point_, cell.collPoint_.J_) * cell.area_;

    std::cout << reducedResult << std::endl;

    const Types::scalar error_real = (preciseResult - reducedResult).real().norm() / preciseResult.real().norm();
    const Types::scalar error_imag = (preciseResult - reducedResult).imag().norm() / preciseResult.imag().norm();
    ASSERT_NEAR(error_real, 0, 1e-2);
    ASSERT_NEAR(error_imag, 0, 1e-2);
}
