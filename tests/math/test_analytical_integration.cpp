//
// Created by evgen on 23.08.24.
//
#include "gtest/gtest.h"

#include "math/integration/analytical/SingularIntegration.hpp"

using namespace EMW;
using namespace EMW::Types;

class ANALYTICAL_INTEGRATION_TESTS : public testing::Test {
};

TEST_F(ANALYTICAL_INTEGRATION_TESTS, ONE_DIV_R_OVER_CELL) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::vector<Mesh::point_t> points = {
            Mesh::point_t{-1, -1, 0},
            Mesh::point_t{1, -1, 0},
            Mesh::point_t{1, 1, 0},
            Mesh::point_t{-1, 1, 0}
    };
    Mesh::IndexedCell cell{{0, 1, 2, 3}, points};

    scalar res = Math::AnalyticalIntegration::integrate_1_div_r(Mesh::point_t{0, 0, 0}, cell);
    ASSERT_NEAR(1 * res, 8 * std::asinh(1.), 1e-14);
}
