//
// Created by evgen on 14.02.24.
//

#include "gtest/gtest.h"

#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/integration/newton_cotess/Rectangular.hpp"
#include "mesh/MeshTypes.hpp"
#include "operators/Functions.hpp"
#include "operators/OperatorK.hpp"
#include "types/Types.hpp"

using namespace EMW;

template <Types::index N, Types::index M> using GL2 = DefiniteIntegrals::GaussLegendre::Quadrature<N, M>;

template <Types::index N> using GL = DefiniteIntegrals::GaussLegendre::Quadrature<N>;

template <Types::index N, Types::index M> using R2 = DefiniteIntegrals::NewtonCotess::Quadrature<N, M>;

template <Types::index N> using R = DefiniteIntegrals::NewtonCotess::Quadrature<N>;

class K1_TESTS : public testing::Test {
  protected:
    Types::complex_d k{1, 0};
};

TEST_F(K1_TESTS, SINGULARITY_EXTRACTION) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::vector<Mesh::point_t> points = {Mesh::point_t{-1, -1, 0}, Mesh::point_t{1, -1, 0},
                                                      Mesh::point_t{1, 1, 0}, Mesh::point_t{-1, 1, 0}};
    Mesh::IndexedCell cell{{0, 1, 2, 3}, points};
    // объявляем точку, в которой интегрируем
    const Mesh::point_t point{0, 0, 0};

    const Types::complex_d result_with_extration =
        OperatorK::detail::K1OverSingularCellRnDWithSingularityExtraction<GL2<8, 8>>(point, cell, k);
    std::cout << result_with_extration << std::endl;
}

TEST_F(K1_TESTS, COMPARISON_OF_INTEGRALS) {
    // объявляем ячейку, по которой происходит интегрирование
    const Types::scalar h = 0.01;
    const Containers::vector<Mesh::point_t> points = {Mesh::point_t{-1, -1, 0} * h, Mesh::point_t{1, -1, 0} * h,
                                                      Mesh::point_t{1, 1, 0} * h, Mesh::point_t{-1, 1, 0} * h};
    Mesh::IndexedCell cell{{0, 1, 2, 2}, points};
    // объявляем точку, в которой интегрируем
    const Mesh::point_t point{0, 0, 0};

    const Types::complex_d result_with_extration =
        OperatorK::detail::K1OverSingularCellRnDWithSingularityExtraction<GL2<8, 8>>(point, cell, k);
    const Types::complex_d result_without =
        OperatorK::detail::K1OverSingularCellReducedAndDivided<GL2<8, 8>>(point, cell, k);
    std::cout.precision(17);
    std::cout << result_with_extration << std::endl;
    std::cout << result_without << std::endl;
}

TEST_F(K1_TESTS, TEST_MIDDLE_ZONE) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::vector<Mesh::point_t> points = {Mesh::point_t{-0.05, -0.05, 0}, Mesh::point_t{0.05, -0.05, 0},
                                                      Mesh::point_t{0.05, 0.05, 0}, Mesh::point_t{-0.05, 0.05, 0}};
    Mesh::IndexedCell cell{{0, 1, 2, 3}, points};
    const Types::Vector3c j = Types::Vector3c{Types::complex_d{0, 0}, Types::complex_d{1, 0}, Types::complex_d{0, 0}};
    // объявляем точку, в которой интегрируем
    const Mesh::point_t point{100, 0, 0};

    // точный расчет интеграла по ячейке
    const Types::Vector3c preciseResult = OperatorK::K0OverSingularCell<GL<8>>(point, j, cell, k) +
                                          k * k * OperatorK::K1OverSingularCellDivided<GL2<8, 8>>(point, j, cell, k);

    std::cout << preciseResult << std::endl;

    // расчет с упрощенным ядром
    Types::Vector3c reducedResult =
        Helmholtz::reducedK_kernel(k, point, cell.collPoint_, j) * cell.area_;

    std::cout << reducedResult << std::endl;

    const Types::scalar error_real = (preciseResult - reducedResult).real().norm() / preciseResult.real().norm();
    const Types::scalar error_imag = (preciseResult - reducedResult).imag().norm() / preciseResult.imag().norm();
    ASSERT_NEAR(error_real, 0, 1e-2);
    ASSERT_NEAR(error_imag, 0, 1e-2);
}
