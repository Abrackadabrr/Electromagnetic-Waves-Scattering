//
// Created by evgen on 23.08.24.
//

#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/MathConstants.hpp"

#include <cmath>

namespace EMW::Math::Integration::Analytical {

Types::scalar integrate_1_div_r(const Mesh::point_t &r, const Mesh::IndexedCell &cell) {
    return integrate_1_div_r(r, cell.getVertexAsArray(), cell.normal);
}

Types::scalar self_newtonian_energy_over_cube(Types::scalar side) {
    // const Types::scalar sqrt2 = std::sqrt(static_cast<Types::scalar>(2));
    // const Types::scalar sqrt3 = std::sqrt(static_cast<Types::scalar>(3));
    // const auto pi = Math::Constants::PI<Types::scalar>();
    //
    // const Types::scalar C =
    //     std::log((sqrt2 + 1) * (sqrt3 + 2))
    //     - pi / 3
    //     - (2 * sqrt3 - sqrt2 - 1) / 5;
    constexpr Types::scalar C_mul_2 = 1.88231264438966;
    const Types::scalar L5 = (side * side) * (side * side) * side;

    return L5 * C_mul_2;
}

namespace detail {

inline Types::scalar hypot3(Types::scalar x, Types::scalar y, Types::scalar z) {
    return std::hypot(std::hypot(x, y), z);
}

inline Types::scalar safe_log(Types::scalar value) {
    constexpr Types::scalar eps = std::numeric_limits<Types::scalar>::epsilon();

    if (std::abs(value) < eps) {
        return Types::scalar(0);
    }

    return std::log(value);
}

inline Types::scalar atan_ratio(Types::scalar num, Types::scalar den) {
    const Types::scalar s = std::signbit(den) ? Types::scalar(-1) : Types::scalar(+1);

    return std::atan2(num * s, std::abs(den));
}

inline Types::scalar primitive_F(Types::scalar x, Types::scalar y, Types::scalar z) {
    const Types::scalar R = hypot3(x, y, z);

    if (R == Types::scalar(0)) {
        return Types::scalar(0);
    }

    Types::scalar result = Types::scalar(0);

    result += x * y * safe_log(z + R);
    result += y * z * safe_log(x + R);
    result += z * x * safe_log(y + R);

    const Types::scalar xR = x * R;
    const Types::scalar yR = y * R;
    const Types::scalar zR = z * R;

    result -= Types::scalar(0.5) * x * x * atan_ratio(y * z, xR);

    result -= Types::scalar(0.5) * y * y * atan_ratio(z * x, yR);

    result -= Types::scalar(0.5) * z * z * atan_ratio(x * y, zR);

    return result;
}

} // namespace detail

Types::scalar newtonian_potential_of_parallelepiped(const Types::point_t &point, Types::scalar a, Types::scalar b,
                                                    Types::scalar c) {
    const Types::scalar X = point.x();
    const Types::scalar Y = point.y();
    const Types::scalar Z = point.z();

    const Types::scalar xs[2] = {-a - X, a - X};

    const Types::scalar ys[2] = {-b - Y, b - Y};

    const Types::scalar zs[2] = {-c - Z, c - Z};

    Types::scalar sum = Types::scalar(0);
    Types::scalar compensation = Types::scalar(0);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {

                const Types::scalar sign = ((i + j + k) & 1) ? Types::scalar(-1) : Types::scalar(+1);

                const Types::scalar value = sign * detail::primitive_F(xs[i], ys[j], zs[k]);

                //
                // Kahan summation
                //
                const Types::scalar corrected = value - compensation;

                const Types::scalar new_sum = sum + corrected;

                compensation = (new_sum - sum) - corrected;

                sum = new_sum;
            }
        }
    }
    return -sum;
}
}