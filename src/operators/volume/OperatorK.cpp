//
// Created by evgen on 17.01.2026.
//

#include "OperatorK.hpp"

#include "math/integration/Quadrature.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"

#include "operators/Functions.hpp"

#include <iostream>

namespace EMW::Operators::Volume {
namespace Gl = DefiniteIntegrals::GaussLegendre;

Types::Matrix3c operator_K_over_cube_mesh::matrix_2_coef(Types::index k, Types::index p) const {
    Types::Matrix3c result = Types::Matrix3c::Zero();
    const auto &cube_k = mesh.getCells()[k];
    const auto &cube_p = mesh.getCells()[p];
    const auto h = mesh.h();

    for (Types::index i = 0; i < 3; i++) {
        Containers::array<Mesh::IndexedCell, 2> faces_k;
        faces_k[0] = cube_k.getFace(static_cast<Mesh::VolumeCells::IndexedCube::Axis>(i),
                                    Mesh::VolumeCells::IndexedCube::Direction::Minus, mesh.getNodes());
        faces_k[1] = cube_k.getFace(static_cast<Mesh::VolumeCells::IndexedCube::Axis>(i),
                                    Mesh::VolumeCells::IndexedCube::Direction::Plus, mesh.getNodes());
        for (Types::index j = 0; j < 3; j++) {
            Containers::array<Mesh::IndexedCell, 2> faces_p;
            faces_p[0] = cube_p.getFace(static_cast<Mesh::VolumeCells::IndexedCube::Axis>(i),
                                        Mesh::VolumeCells::IndexedCube::Direction::Minus, mesh.getNodes());
            faces_p[1] = cube_p.getFace(static_cast<Mesh::VolumeCells::IndexedCube::Axis>(i),
                                        Mesh::VolumeCells::IndexedCube::Direction::Plus, mesh.getNodes());

            // и тут нужно взять четыре одинаковых по вайбу интеграла
            for (auto &&face_k : faces_k)
                for (auto &&face_p : faces_p) {
                    // Делаем.
                    if ((cube_p.center_ - cube_p.center_).norm() < 6 * h) {
                        // то интегрируемся с выделением особенности
                        // 1. Интеграл от ньютонова потенциала 2д ячейки по ней же самой
                        const auto analytical_integrand = [&face_k, &face_p](Types::scalar x, Types::scalar y) {
                            const auto point = face_k.parametrization(x, y);
                            return Math::AnalyticalIntegration::integrate_1_div_r(point, face_p) * face_k.multiplier(
                                       x, y);
                        };
                        result(i, j) += Math::Constants::inverse_4PI<Types::scalar>() * DefiniteIntegrals::integrate<
                            DefiniteIntegrals::GaussLegendre::Quadrature<
                                4, 4>>(analytical_integrand, {0, 0}, {1, 1});
                        // 2. Интеграл от ограниченного остатка
                        const auto residual_integrand = [&face_k, &face_p, wn = wave_number](
                            Types::scalar x1, Types::scalar y1, Types::scalar x2, Types::scalar y2) {
                            const auto x = face_p.parametrization(x1, y1);
                            const auto y = face_k.parametrization(x2, y2);
                            return Helmholtz::F_bounded_part(wn, x, y) * face_p.multiplier(x1, y1) * face_k.multiplier(
                                       x2, y2);
                        };

                        result(i, j) += DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<
                            2, 2, 2, 2>>(
                            residual_integrand, {0, 0, 0, 0}, {1, 1, 1, 1});
                    } else {
                        // Иначе интегрируемся без выделения особенности сразу
                        const auto integrand = [&face_k, &face_p, wn = wave_number](
                            Types::scalar x1, Types::scalar y1, Types::scalar x2, Types::scalar y2) {
                            const auto x = face_p.parametrization(x1, y1);
                            const auto y = face_k.parametrization(x2, y2);
                            return Helmholtz::F(wn, x, y) * face_p.multiplier(x1, y1) * face_k.multiplier(x2, y2);
                        };

                        result(i, j) += DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<
                            2, 2, 2, 2>>(
                            integrand, {0, 0, 0, 0}, {1, 1, 1, 1});
                    }
                }

        }
    }
    return result;
}

Types::complex_d operator_K_over_cube_mesh::matrix_3_coef(Types::index k, Types::index p) const {
    const auto &cube_k = mesh.getCells()[k];
    const auto &cube_p = mesh.getCells()[p];
    const auto h = mesh.h();
    const Types::complex_d k2 = wave_number * wave_number;

    // if (k == p)
    // чтобы это зафорсить, нужно обобщить формулу для самоэнергии
    // куба на параллелограмм и тогда все будить чики бамбони

    if ((cube_k.center_ - cube_p.center_).norm() < 6 * h) {
        // интегрирование с выделением особенности
        // Ограниченная часть от функции
        const auto integrand_bounded_part = [wn = wave_number](Types::scalar x1, Types::scalar y1, Types::scalar z1,
                                                               Types::scalar x2, Types::scalar y2,
                                                               Types::scalar z2) {
            return Helmholtz::F_bounded_part(wn, {x1, y1, z1}, {x2, y2, z2});
        };
        //  Ньютонов потенциал куба в правильном формате
        const auto potential_of_cube = [this, k](Types::scalar x, Types::scalar y, Types::scalar z) {
            return newton_potential_of_cube(k, {x, y, z});
        };

        return k2 * (
                   DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<2, 2, 2, 2, 2, 2>>(
                       integrand_bounded_part, {mesh.leftDownCorner(k)[0], mesh.leftDownCorner(k)[1],
                                                mesh.leftDownCorner(k)[2],
                                                mesh.leftDownCorner(p)[0], mesh.leftDownCorner(p)[1],
                                                mesh.leftDownCorner(p)[2]},
                       {mesh.dx(), mesh.dy(), mesh.dz(), mesh.dx(), mesh.dy(), mesh.dz()}) +
                   Math::Constants::inverse_4PI<Types::scalar>() *
                   DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<4, 4, 4>>(
                       potential_of_cube, {mesh.leftDownCorner(p)[0], mesh.leftDownCorner(p)[1],
                                           mesh.leftDownCorner(p)[2]},
                       {mesh.dx(), mesh.dy(), mesh.dz()}));
    }

    // Интегрирование без выделения особенности:
    // просто берем фундаментальное решение уравнения Гельмгольца
    // и интегрируем его по двум кубам.
    const auto integrand = [wn = wave_number](Types::scalar x1, Types::scalar y1, Types::scalar z1,
                                              Types::scalar x2, Types::scalar y2,
                                              Types::scalar z2) {
        return Helmholtz::F(wn, {x1, y1, z1}, {x2, y2, z2});
    };
    return k2 * DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<2, 2, 2, 2, 2, 2>>(
               integrand, {mesh.leftDownCorner(k)[0], mesh.leftDownCorner(k)[1], mesh.leftDownCorner(k)[2],
                           mesh.leftDownCorner(p)[0], mesh.leftDownCorner(p)[1], mesh.leftDownCorner(p)[2]},
               {mesh.dx(), mesh.dy(), mesh.dz(), mesh.dx(), mesh.dy(), mesh.dz()});

}

Types::scalar operator_K_over_cube_mesh::newton_potential_of_cube(Types::index k, Types::point_t r) const {
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
    return DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(integrand, {0}, {dz / 2}) +
           DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(
               integrand, {dz / 2}, {dz / 2});
}

Types::MatrixXc operator_K_over_cube_mesh::get_galerkin_matrix() const {
    const Types::index n_cubes = mesh.getCells().size();
    Types::MatrixXc result = Types::MatrixXc::Zero(3 * n_cubes, 3 * n_cubes);
    for (auto k = 0u; k < n_cubes; ++k) {
        const auto& cube_k = mesh.getCells()[k];
        for (auto p = 0u; p < n_cubes; ++p) {
            // считаем поверхностную часть
            const auto volume_res = matrix_3_coef(k, p);
            result.block(3 * k, 3 * p, 3, 3) = matrix_2_coef(k, p);
            // и подправляем общую матрицу
            result(3 * k, 3 * p) += volume_res;
            result(3 * k + 1, 3 * p + 1) += volume_res;
            result(3 * k + 2, 3 * p + 2) += volume_res;
            std::cout << k << ' ' << p << '\n';
        }
    }
    return result;
}

void operator_K_over_cube_mesh::get_galerkin_matrix_inplace(Types::MatrixXc * p_mat) const {
    const Types::index n_cubes = mesh.getCells().size();
    *p_mat = Types::MatrixXc::Zero(3 * n_cubes, 3 * n_cubes);
    for (auto k = 0u; k < n_cubes; ++k) {
        const auto& cube_k = mesh.getCells()[k];
        for (auto p = 0u; p < n_cubes; ++p) {
            // считаем поверхностную часть
            const auto volume_res = matrix_3_coef(k, p);
            (*p_mat).block(3 * k, 3 * p, 3, 3) = matrix_2_coef(k, p);
            // и подправляем общую матрицу
            (*p_mat)(3 * k, 3 * p) += volume_res;
            (*p_mat)(3 * k + 1, 3 * p + 1) += volume_res;
            (*p_mat)(3 * k + 2, 3 * p + 2) += volume_res;
        }
    }
}
}
