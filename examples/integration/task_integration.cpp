//
// Created by evgen on 31.01.24.
//

#include "types/Types.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "integration/gauss_quadrature/GaussQuadrature.hpp"

using EMW::Types::scalar;
using namespace EMW;

scalar f4(scalar x, scalar y) {
    return x * x * y * y * y * y;;
}

scalar f5(scalar x, scalar y) {
    return x * x * y * y * y * y * y;
}

scalar f(scalar x) {
    return (x + 1) * x * x * x * x;
}

int main() {
    const scalar y = 1;
    auto f_l = [y] (scalar x) -> scalar {return (x + y) * x * x * x * x;};
    const auto result = DefiniteIntegrals::integrate<DefiniteIntegrals::Quadrature<3>>(f_l, {0}, {2});
    assert(std::abs(2 * 2 * 2 * 2 * 2 * 2. / 6 + 2 * 2 * 2 * 2 * 2. / 5 - result) < 1e-14);
}
