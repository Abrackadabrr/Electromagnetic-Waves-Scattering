//
// Created by evgen on 17.01.2026.
//

#include "OperatorK.hpp"

#include "math/integration/Quadrature.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"

#include "operators/Functions.hpp"

#include <bitset>
#include <iostream>
#include <bits/fs_fwd.h>

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
            faces_p[0] = cube_p.getFace(static_cast<Mesh::VolumeCells::IndexedCube::Axis>(j),
                                        Mesh::VolumeCells::IndexedCube::Direction::Minus, mesh.getNodes());
            faces_p[1] = cube_p.getFace(static_cast<Mesh::VolumeCells::IndexedCube::Axis>(j),
                                        Mesh::VolumeCells::IndexedCube::Direction::Plus, mesh.getNodes());

            // и тут нужно взять четыре одинаковых по вайбу интеграла
            for (size_t face_k_idx = 0; face_k_idx < 2; face_k_idx++)
                for (size_t face_p_idx = 0; face_p_idx < 2; face_p_idx++) {
                    // Делаем.
                    Types::scalar multiplier = face_k_idx == face_p_idx ? 1 : -1;
                    auto face_k = faces_k[face_k_idx];
                    auto face_p = faces_p[face_p_idx];

                    if ((cube_k.center_ - cube_p.center_).norm() < 6 * h) {
                        // то интегрируемся с выделением особенности
                        // 1. Интеграл от ньютонова потенциала 2д ячейки
                        const auto analytical_integrand = [&face_k, &face_p](Types::scalar x, Types::scalar y) {
                            const auto point = face_k.parametrization(x, y);
                            const auto integrand_value = Math::AnalyticalIntegration::integrate_1_div_r(point, face_p);
                            return integrand_value * face_k.multiplier(
                                       x, y);
                        };
                        const auto singular_part =
                            Math::Constants::inverse_4PI<Types::scalar>() * DefiniteIntegrals::integrate<
                                DefiniteIntegrals::GaussLegendre::Quadrature<
                                    3, 3>>(analytical_integrand, {0, 0}, {1, 1});
                        result(i, j) += multiplier * singular_part;

                        // 2. Интеграл от ограниченного остатка
                        const auto residual_integrand = [&face_k, &face_p, wn = wave_number](
                            Types::scalar x1, Types::scalar y1, Types::scalar x2, Types::scalar y2) {
                            const auto x = face_p.parametrization(x1, y1);
                            const auto y = face_k.parametrization(x2, y2);
                            return Helmholtz::F_bounded_part(wn, x, y) * face_p.multiplier(x1, y1) * face_k.multiplier(
                                       x2, y2);
                        };
                        const auto regular_part = DefiniteIntegrals::integrate<
                            DefiniteIntegrals::GaussLegendre::Quadrature<
                                3, 3, 3, 3>>(
                            residual_integrand, {0, 0, 0, 0}, {1, 1, 1, 1});
                        result(i, j) += multiplier * regular_part;

                    } else {
                        // Иначе интегрируемся без выделения особенности сразу
                        const auto integrand = [&face_k, &face_p, wn = wave_number](
                            Types::scalar x1, Types::scalar y1, Types::scalar x2, Types::scalar y2) {
                            const auto x = face_p.parametrization(x1, y1);
                            const auto y = face_k.parametrization(x2, y2);
                            return Helmholtz::F(wn, x, y) * face_p.multiplier(x1, y1) * face_k.multiplier(x2, y2);
                        };

                        result(i, j) += multiplier * DefiniteIntegrals::integrate<
                            DefiniteIntegrals::GaussLegendre::Quadrature<
                                3, 3, 3, 3>>(
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
        const auto potential_of_cube_k = [this, k](Types::scalar x, Types::scalar y, Types::scalar z) {
            return newton_potential_of_cube(k, {x, y, z});
        };

        const auto regular_part =
            DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<
                2, 2, 2, 2, 2, 2>>(
                integrand_bounded_part, {mesh.leftDownCorner(k)[0], mesh.leftDownCorner(k)[1],
                                         mesh.leftDownCorner(k)[2],
                                         mesh.leftDownCorner(p)[0], mesh.leftDownCorner(p)[1],
                                         mesh.leftDownCorner(p)[2]},
                {mesh.dx(), mesh.dy(), mesh.dz(), mesh.dx(), mesh.dy(), mesh.dz()});
        const auto singular_part =
            Math::Constants::inverse_4PI<Types::scalar>() *
            DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<
                3, 3, 3>>(
                potential_of_cube_k, {mesh.leftDownCorner(p)[0], mesh.leftDownCorner(p)[1],
                                      mesh.leftDownCorner(p)[2]},
                {mesh.dx(), mesh.dy(), mesh.dz()});
        return k2 * (regular_part + singular_part);
    }

    // Интегрирование без выделения особенности:
    // просто берем фундаментальное решение уравнения Гельмгольца
    // и интегрируем его по двум кубам.
    const auto integrand = [wn = wave_number](Types::scalar x1, Types::scalar y1, Types::scalar z1,
                                              Types::scalar x2, Types::scalar y2,
                                              Types::scalar z2) {
        return Helmholtz::F(wn, {x1, y1, z1}, {x2, y2, z2});
    };
    return k2 * DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<3, 3, 3, 3, 3, 3>>(
               integrand, {mesh.leftDownCorner(k)[0], mesh.leftDownCorner(k)[1], mesh.leftDownCorner(k)[2],
                           mesh.leftDownCorner(p)[0], mesh.leftDownCorner(p)[1], mesh.leftDownCorner(p)[2]},
               {mesh.dx(), mesh.dy(), mesh.dz(), mesh.dx(), mesh.dy(), mesh.dz()});

}

Types::scalar operator_K_over_cube_mesh::newton_potential_of_cube(Types::index k, Types::point_t r) const {
    const auto &left_down_corner_of_cube_k = mesh.leftDownCorner(k);
    const auto dx = mesh.dx();
    const auto dy = mesh.dy();
    const auto dz = mesh.dz();
    // интегрирование старое (квадратура поверх 1/r по квадрату)
    const auto cut_of_cube_integral = [corner = left_down_corner_of_cube_k, dx,
            dy](Types::scalar x, Types::scalar y, Types::scalar z, Types::scalar z_dash) {
        // Интегрирование по сечению куба на высоте z_dash
        Containers::vector<Types::point_t> points1 = {
            {corner.x(), corner.y(), corner.z() + z_dash},
            {corner.x() + dx, corner.y(), corner.z() + z_dash},
            {corner.x() + dx, corner.y() + dy, corner.z() + z_dash},
            {corner.x(), corner.y() + dy, corner.z() + z_dash},
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

Types::Matrix3c operator_K_over_cube_mesh::galerkin_block_for_cubes(size_t k, size_t p) const {
    const Types::index n_cubes = mesh.getCells().size();
    // считаем поверхностную часть
    const auto volume_res = matrix_3_coef(k, p);
    Types::Matrix3c result = -matrix_2_coef(k, p);
    // и подправляем общую матрицу
    result(0, 0) += volume_res;
    result(1, 1) += volume_res;
    result(2, 2) += volume_res;

    return result;
}

Types::MatrixXc operator_K_over_cube_mesh::get_galerkin_matrix() const {
    const Types::index n_cubes = mesh.getCells().size();
    Types::MatrixXc result = Types::MatrixXc::Zero(3 * n_cubes, 3 * n_cubes);
#pragma omp parallel for simd schedule(static) num_threads(14)
    for (auto p = 0u; p < n_cubes; ++p) {
        for (auto k = 0u; k < n_cubes; ++k) {
            // считаем поверхностную часть
            const auto volume_res = matrix_3_coef(k, p);
            const auto surface_res = matrix_2_coef(k, p);
            result.block(3 * k, 3 * p, 3, 3) = -surface_res;
            // и подправляем общую матрицу
            result(3 * k, 3 * p) += volume_res;
            result(3 * k + 1, 3 * p + 1) += volume_res;
            result(3 * k + 2, 3 * p + 2) += volume_res;
        }
    }
    return result;
}

void operator_K_over_cube_mesh::get_galerkin_matrix_inplace(Types::MatrixXc *p_mat) const {
    const Types::index n_cubes = mesh.getCells().size();
    *p_mat = Types::MatrixXc::Zero(3 * n_cubes, 3 * n_cubes);
    for (auto k = 0u; k < n_cubes; ++k) {
        const auto &cube_k = mesh.getCells()[k];
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

Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d>
operator_K_over_cube_mesh::compute_galerkin_matrix(Types::scalar l1_basis_function_norm) const {
    const size_t first_layer_toeplitz = mesh.nx() - 1;
    const size_t second_layer_toeplitz = mesh.ny() - 1;
    const size_t third_layer_toeplitz = mesh.nz() - 1;
    decltype(auto) result =
        Math::LinAgl::Matrix::ZeroTripleToeplitzBlock<Types::complex_d>(first_layer_toeplitz, second_layer_toeplitz,
                                                                        third_layer_toeplitz, 3);
    const Types::scalar l2_basis_fn_norm_sqr = l1_basis_function_norm * l1_basis_function_norm;
    // Циклы для расчета трижды теплицевой матрицы
    // Эта штука так же долго считает, но занимает уже меньше памяти
    for (size_t i3 = 0; i3 < third_layer_toeplitz; ++i3)
        for (size_t j3 = 0; j3 < third_layer_toeplitz; ++j3)
            for (size_t i2 = 0; i2 < second_layer_toeplitz; ++i2)
                for (size_t j2 = 0; j2 < second_layer_toeplitz; ++j2)
                    for (size_t i1 = 0; i1 < first_layer_toeplitz; ++i1)
                        for (size_t j1 = 0; j1 < first_layer_toeplitz; ++j1) {
                            auto& working_block = result.get_block(i3, j3).
                                   get_block(i2, j2).
                                   get_block(i1, j1);
                            // Ускорение заполнения матрицы за счет отсутствия
                            // пересчёта одинаковых блоков
                            // TODO: сделать нормальный расчет, то есть аналитически вывести все формулки
                            if (working_block.norm() == 0) {
                                // Ищем кубы по трёхмерному индексу
                                const auto idx1 = mesh.cube_idx(i1, i2, i3);
                                const auto idx2 = mesh.cube_idx(j1, j2, j3);
                                // Счёт
                                const auto volume_res = matrix_3_coef(idx1, idx2);
                                Types::Matrix3c surface_res = -matrix_2_coef(idx1, idx2);
                                // и подправляем общую матрицу
                                surface_res(0, 0) += volume_res;
                                surface_res(1, 1) += volume_res;
                                surface_res(2, 2) += volume_res;

                                working_block = surface_res / l2_basis_fn_norm_sqr;
                            }
                        }
    return result;
}

#if 0
Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d>
operator_K_over_cube_mesh::compute_galerkin_matrix_wise() const {
    const size_t first_layer_toeplitz = mesh.nx() - 1;
    const size_t second_layer_toeplitz = mesh.ny() - 1;
    const size_t third_layer_toeplitz = mesh.nz() - 1;
    decltype(auto) result =
        Math::LinAgl::Matrix::ZeroTripleToeplitzBlock<Types::complex_d>(first_layer_toeplitz, second_layer_toeplitz,
                                                                        third_layer_toeplitz, 3);
    // Циклы для расчета трижды теплицевой матрицы

    // Функция, которая будет необходима для расчета возьми блоков в тёплицевой структуре
    const auto calc_index_of_cube = [](size_t i1, size_t i2, size_t i3, uint8_t mask) {
        return std::array{ (mask & 1) * i1, (mask & 2) * i2, (mask & 4) * i3};
    };

    for (Types::integer i3 = -third_layer_toeplitz + 1; i3 < third_layer_toeplitz; ++i3)
        for (Types::integer i2 = -second_layer_toeplitz + 1; i2 < second_layer_toeplitz; ++i2)
            for (Types::integer i1 = -first_layer_toeplitz + 1; i1 < first_layer_toeplitz; ++i1) {
                // Находим какой-то куб внутри сетки
                const auto current_cube_idx = mesh.cube_idx(i1, i2, i3);
                // Считаем два варианта: взаимодействие reference_cube_idx на current_cube_idx
                // и наоборот
                const auto ref_to_cur = galerkin_block_for_cubes(reference_cube_idx, current_cube_idx);
                const auto cur_to_ref = galerkin_block_for_cubes(reference_cube_idx, current_cube_idx);
                // И теперь пытаемся понять куда этот расчет запихать в контейнер
                for (char mask = 0; mask < 8; ++mask) {
                    const auto idxs = calc_index_of_cube(i1, i2, i3, mask);
                }
                // И тут есть восемь возможных вариантов для расчета
                // (как раз (2N_1 - 1)(2N_2 - 1)(2N_3 - 1) уникальных блоков в матрице)
                auto &&col_third_level_block =
                    result.get_block(i3, 0);
                get_block(i2, 0).
                    get_block(i1, 0) = surface_res;

                // Заполняем первую строку
                result.get_block(0, i3).
                       get_block(0, i2).
                       get_block(0, i1) = surface_res;
            }

    return result;
}
#endif
}
