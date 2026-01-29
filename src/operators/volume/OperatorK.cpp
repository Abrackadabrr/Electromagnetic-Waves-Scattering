//
// Created by evgen on 17.01.2026.
//

#include "OperatorK.hpp"

#include "math/integration/Quadrature.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"

namespace EMW::Operators::Volume {

namespace Gl = DefiniteIntegrals::GaussLegendre;

Types::Matrix3c operator_K_over_cube_mesh::matrix_3_coef(Types::index k, Types::index p) {
    const auto &cube_k = mesh.getCells()[k];
    const auto &cube_p = mesh.getCells()[p];
    const auto h = mesh.h();
    const Types::complex_d k2 = wave_number * wave_number;

    if ((cube_k.center_ - cube_p.center_).norm() < 6 * h) {
        // интегрирование с вынесением особенности
    } else {
        // обычное интегрирование по объёму дважды
    }
}

Types::scalar operator_K_over_cube_mesh::newton_potential_over_cube(Types::index k, Types::point_t r) {
    const auto &left_down_corner_of_cube_k = mesh.leftDownCorner(k);
    const auto dx = mesh.dx();
    const auto dy = mesh.dy();
    const auto dz = mesh.dz();
    // интегрирование старое (квадратура в квадратуре в квадратуре)
    const auto cut_of_cube_integral = [corner = left_down_corner_of_cube_k, dx,
                                       dy](Types::scalar x, Types::scalar y, Types::scalar z, Types::scalar z_dash) {
        // Интегрирование по сечению куба на высоте z_dash
        Containers::vector<Types::point_t> points1 = {
            {corner.x(), corner.y(), corner.z() + z_dash},
            {corner.x() + dx, corner.y(), z_dash},
            {corner.x() + dx, corner.y() + dy, z_dash},
            {corner.x(), corner.y() + dy, z_dash},
        };
        const Mesh::IndexedCell cell1({0, 1, 2, 3}, points1);
        return Math::AnalyticalIntegration::integrate_1_div_r(Types::point_t{x, y, z}, cell1);
    };

    const auto integrand = [x = r.x(), y = r.y(), z = r.z(), cut_of_cube_integral](Types::scalar z_dash) {
        return cut_of_cube_integral(x, y, z, z_dash);
    };
    return DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(integrand, {0}, {dz/2}) +
            DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(integrand, {dz/2}, {dz/2});
    };
}
