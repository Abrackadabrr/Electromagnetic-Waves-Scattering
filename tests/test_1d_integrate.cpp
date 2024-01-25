//
// Created by evgen on 25.01.24.
//

#include "gtest/gtest.h"
#include "Types.hpp"

using EMW::Types::scalar;

class TEST_1D_INTEGRATION : public ::testing::Test {
protected:
    scalar prec = 1e-3;
};

TEST_F(TEST_1D_INTEGRATION, ConstantIntegration) {
    ASSERT_NEAR(prec, 0, prec);
}