//
// Created by evgen on 14.06.24.
//

#include "gtest/gtest.h"

#include "types/Types.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "integration/newton_cotess/Rectangular.hpp"
#include "integration/Quadrature.hpp"

using EMW::Types::scalar;
using namespace EMW;
namespace GL = DefiniteIntegrals::GaussLegendre;
namespace NC = DefiniteIntegrals::NewtonCotess;

Types::complex_d constant(scalar x, scalar y) {
    return 1i;
}

scalar f4(scalar x, scalar y) {
    return x * x * y * y * y * y;;
}

scalar f5(scalar x, scalar y) {
    return x * x * y * y * y * y * y;
}

scalar f(scalar x) {
    return (x + 1) * x * x * x * x;
}

class TRIANGULAR_INTEGRATION : public ::testing::Test {
protected:
    scalar prec = 4e-15;
};


TEST(TRIANGULAR_INTEGRATION, CONSTANT) {
    // 1) Make mesh with 1 triangular node
    // 2) Check if integration over it correct
}
