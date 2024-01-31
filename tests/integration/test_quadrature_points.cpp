//
// Created by evgen on 31.01.24.
//

#include "gtest/gtest.h"

#include "Types.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"

using namespace EMW;

TEST(Quadrature, CompilationTest) {
    constexpr auto product = DefiniteIntegrals::Quadrature<6, 8>::nodes;
    for (const auto item : product) {
        const Types::scalar weight = item.weight;
//        std::cout << item.point[1] << ", ";
    }
    ASSERT_TRUE(true);
}
