//
// Created by evgen on 17.01.2026.
//

#include "operators/volume/OperatorK.hpp"

#include "math/integration/NumericalIntegration.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"

#include "operators/Functions.hpp"

#include "math/matrix/decompositions/Decompositions.hpp"

namespace EMW::Operators::Volume {
namespace Gl = DecartIntegration::GaussLegendre;

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

                    if ((cube_k.center_ - cube_p.center_).norm() < nearnes_tresholds * h) {
                        // то интегрируемся с выделением особенности
                        // 1. Интеграл от ньютонова потенциала 2д ячейки
                        const auto analytical_integrand = [&face_k, &face_p](Types::scalar x, Types::scalar y) {
                            const auto point = face_k.parametrization(x, y);
                            const auto integrand_value = Math::Integration::Analytical::integrate_1_div_r(point, face_p);
                            return integrand_value * face_k.multiplier(
                                       x, y);
                        };
                        const auto singular_part =
                            Math::Constants::inverse_4PI<Types::scalar>() * DecartIntegration::integrate<
                                DecartIntegration::NewtonCotess::Quadrature<
                                    2, 2>>(analytical_integrand, {0, 0}, {1, 1});
                        result(i, j) += multiplier * singular_part;

                        // 2. Интеграл от ограниченного остатка
                        const auto residual_integrand = [&face_k, &face_p, wn = wave_number](
                            Types::scalar x1, Types::scalar y1, Types::scalar x2, Types::scalar y2) {
                            const auto x = face_p.parametrization(x1, y1);
                            const auto y = face_k.parametrization(x2, y2);
                            return Helmholtz::F_bounded_part(wn, x, y) * face_p.multiplier(x1, y1) * face_k.multiplier(
                                       x2, y2);
                        };
                        const auto regular_part = DecartIntegration::integrate<
                            DecartIntegration::GaussLegendre::Quadrature<
                                4, 4, 4, 4>>(
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

                        result(i, j) += multiplier * DecartIntegration::integrate<
                            DecartIntegration::GaussLegendre::Quadrature<
                                4, 4, 4, 4>>(
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
    const auto &k_corner = mesh.leftDownCorner(k);
    const auto &p_corner = mesh.leftDownCorner(p);
    const auto h = mesh.h();
    const Types::complex_d k2 = wave_number * wave_number;

    // if (k == p)
    // чтобы это зафорсить, нужно обобщить формулу для самоэнергии
    // куба на параллелограмм и тогда все будить чики бамбони

    if ((cube_k.center_ - cube_p.center_).norm() < nearnes_tresholds * h) {
        // интегрирование с выделением особенности
        // Ограниченная часть от функции
        const auto integrand_bounded_part = [wn = wave_number](Types::scalar x1, Types::scalar y1, Types::scalar z1,
                                                               Types::scalar x2, Types::scalar y2, Types::scalar z2) {
            return Helmholtz::F_bounded_part(wn, {x1, y1, z1}, {x2, y2, z2});
        };
        //  Ньютонов потенциал куба в правильном формате
        const auto potential_of_cube_k = [this, k](Types::scalar x, Types::scalar y, Types::scalar z) {
            return newton_potential_of_cube(k, {x, y, z});
        };

        const auto regular_part =
            DecartIntegration::integrate<DecartIntegration::GaussLegendre::Quadrature<
                2, 2, 2, 2, 2, 2>>(
                integrand_bounded_part, {k_corner.x(), k_corner.y(), k_corner.z(),
                                         p_corner.x(), p_corner.y(), p_corner.z()},
                {mesh.dx(), mesh.dy(), mesh.dz(), mesh.dx(), mesh.dy(), mesh.dz()});
        const auto singular_part =
            Math::Constants::inverse_4PI<Types::scalar>() *
            DecartIntegration::integrate<DecartIntegration::GaussLegendre::Quadrature<
                4, 4, 4>>(
                potential_of_cube_k, {p_corner.x(), p_corner.y(), p_corner.z()},
                {mesh.dx(), mesh.dy(), mesh.dz()});
        return k2 * (regular_part + singular_part);
    }

    // Интегрирование без выделения особенности:
    // просто берем фундаментальное решение уравнения Гельмгольца
    // и интегрируем его по двум кубам.
    const auto integrand = [wn = wave_number](Types::scalar x1, Types::scalar y1, Types::scalar z1,
                                              Types::scalar x2, Types::scalar y2, Types::scalar z2) {
        return Helmholtz::F(wn, {x1, y1, z1}, {x2, y2, z2});
    };
    return k2 * DecartIntegration::integrate<DecartIntegration::GaussLegendre::Quadrature<2, 2, 2, 2, 2, 2>>(
               integrand, {k_corner.x(), k_corner.y(), k_corner.z(), p_corner.x(), p_corner.y(), p_corner.z()},
               {mesh.dx(), mesh.dy(), mesh.dz(), mesh.dx(), mesh.dy(), mesh.dz()});
}

Types::Matrix3c operator_K_over_cube_mesh::far_zone_interaction(Types::index k, Types::index p) const {
    // Просто интегрируем выражение для поля в дальней зоне
    const auto &k_corner = mesh.leftDownCorner(k);
    const auto &p_corner = mesh.leftDownCorner(p);
    Types::Matrix3c interaction_block = Types::Matrix3c::Zero();
    for (size_t i = 0; i < 3; ++i) {
        Types::Vector3c j = Types::Vector3c::Zero();
        j[i] = {1., 0};
        const auto integrand = [wn = wave_number, j](Types::scalar x1, Types::scalar y1, Types::scalar z1,
                                                     Types::scalar x2, Types::scalar y2, Types::scalar z2) {
            return Helmholtz::far_zone_integral_kernel(wn, Types::point_t{x1, y1, z1} - Types::point_t{x2, y2, z2}, j);
        };
        interaction_block.col(i) = DecartIntegration::integrate<DecartIntegration::GaussLegendre::Quadrature<
            2, 2, 2, 2, 2, 2>>(
            integrand, {k_corner.x(), k_corner.y(), k_corner.z(), p_corner.x(), p_corner.y(), p_corner.z()},
            {mesh.dx(), mesh.dy(), mesh.dz(), mesh.dx(), mesh.dy(), mesh.dz()});
    }
    return Math::Constants::inverse_4PI<Types::scalar>() * interaction_block;
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
        return Math::Integration::Analytical::integrate_1_div_r(Types::point_t{x, y, z}, cell1);
    };

    const auto integrand = [x = r.x(), y = r.y(), z = r.z(), cut_of_cube_integral](Types::scalar z_dash) {
        return cut_of_cube_integral(x, y, z, z_dash);
    };
    return DecartIntegration::integrate<DecartIntegration::GaussLegendre::Quadrature<3>>(integrand, {0}, {dz / 2}) +
           DecartIntegration::integrate<DecartIntegration::GaussLegendre::Quadrature<3>>(
               integrand, {dz / 2}, {dz / 2});
}

Types::Matrix3c operator_K_over_cube_mesh::galerkin_block_for_cubes(size_t k, size_t p) const {
#if 0
    // 1. Если кубы далеко, то считаем через far_zone
    if (mesh.distance(k, p) > 10)
        // Ну например 10. ...
        return far_zone_interaction(k, p);
#endif

    // 2. Иначе считаем через преобразование сингулярного оператора
    const auto volume_res = matrix_3_coef(k, p);
    Types::Matrix3c surface_res = -matrix_2_coef(k, p);
    // и подправляем общую матрицу
    surface_res(0, 0) += volume_res;
    surface_res(1, 1) += volume_res;
    surface_res(2, 2) += volume_res;

    return surface_res;
}

Types::MatrixXc operator_K_over_cube_mesh::compute_galerkin_matrix_dense(Types::scalar basis_function_module) const {
    const Types::index n_cubes = mesh.getCells().size();
    Types::MatrixXc result = Types::MatrixXc::Zero(3 * n_cubes, 3 * n_cubes);
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
    return result * (basis_function_module * basis_function_module);
}

void operator_K_over_cube_mesh::compute_galerkin_matrix_dense_inplace(Types::MatrixXc *p_mat) const {
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
operator_K_over_cube_mesh::compute_galerkin_matrix(Types::scalar basis_function_module) const {
    const size_t first_layer_toeplitz = mesh.nx() - 1;
    const size_t second_layer_toeplitz = mesh.ny() - 1;
    const size_t third_layer_toeplitz = mesh.nz() - 1;
    decltype(auto) result =
        Math::LinAgl::Matrix::ZeroTripleToeplitzBlock<Types::complex_d>(first_layer_toeplitz, second_layer_toeplitz,
                                                                        third_layer_toeplitz, 3);
    const Types::scalar basis_fn_module_sqr = basis_function_module * basis_function_module;
    // Циклы для расчета трижды теплицевой матрицы
    for (size_t i3 = 0; i3 < third_layer_toeplitz; ++i3)
        for (size_t j3 = 0; j3 < third_layer_toeplitz; ++j3) {
            for (size_t i2 = 0; i2 < second_layer_toeplitz; ++i2)
                for (size_t i1 = 0; i1 < first_layer_toeplitz; ++i1)
                    for (size_t j2 = 0; j2 < second_layer_toeplitz; ++j2)
                        for (size_t j1 = 0; j1 < first_layer_toeplitz; ++j1) {
                            auto &&working_block = result.get_block(i3, j3).
                                                          get_block(i2, j2).
                                                          get_block(i1, j1);
                            // Ускорение заполнения матрицы за счет отсутствия
                            // пересчёта одинаковых блоков
                            // TODO: сделать нормальный расчет, то есть аналитически вывести все формулки
                            if (working_block.norm() == 0) {
                                // Ищем кубы по трёхмерному индексу
                                const auto idx1 = mesh.cube_idx(i1, i2, i3);
                                const auto idx2 = mesh.cube_idx(j1, j2, j3);
                                working_block = galerkin_block_for_cubes(idx1, idx2) * basis_fn_module_sqr;
                            }
                        }
        }
    return result;
}

[[nodiscard]] Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d>
operator_K_over_cube_mesh::compute_galerkin_matrix(
    Idx3d start_i, Idx3d start_j, Idx3d sizes, Types::scalar basis_fn_module) const {

    // Делаем нулевую трижды тёплицеву матрицу
    const size_t first_layer_toeplitz = sizes.Nx;
    const size_t second_layer_toeplitz = sizes.Ny;
    const size_t third_layer_toeplitz = sizes.Nz;
    constexpr size_t inner_size = 3;
    decltype(auto) result =
        Math::LinAgl::Matrix::ZeroTripleToeplitzBlock<Types::complex_d>(first_layer_toeplitz, second_layer_toeplitz,
                                                                        third_layer_toeplitz, inner_size);
    const Types::scalar basis_fn_module_sqr = basis_fn_module * basis_fn_module;

    // Циклы для расчета трижды теплицевой матрицы
    for (size_t i3 = 0; i3 < third_layer_toeplitz; ++i3)
        for (size_t j3 = 0; j3 < third_layer_toeplitz; ++j3)
            for (size_t i2 = 0; i2 < second_layer_toeplitz; ++i2)
                for (size_t i1 = 0; i1 < first_layer_toeplitz; ++i1)
                    for (size_t j2 = 0; j2 < second_layer_toeplitz; ++j2)
                        for (size_t j1 = 0; j1 < first_layer_toeplitz; ++j1) {
                            auto &&working_block = result.get_block(i3, j3).
                                                          get_block(i2, j2).
                                                          get_block(i1, j1);
                            // Ускорение заполнения матрицы за счет отсутствия
                            // пересчёта одинаковых блоков
                            // TODO: сделать нормальный расчет, то есть аналитически вывести все формулки
                            if (working_block.norm() == 0) {
                                // Ищем кубы по трёхмерному индексу
                                const auto idx1 = mesh.cube_idx(start_i.Nx + i1, start_i.Ny + i2, start_i.Nz + i3);
                                const auto idx2 = mesh.cube_idx(start_j.Nx + j1, start_j.Ny + j2, start_j.Nz + j3);
                                // Счёт
                                working_block = galerkin_block_for_cubes(idx1, idx2) * basis_fn_module_sqr;
                            }
                        }
    return result;
}

operator_K_over_cube_mesh::matrix_and_permutation<Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d>>
operator_K_over_cube_mesh::compute_galerkin_matrix_custom_blocksize(size_t Nx, size_t Ny, size_t Nz,
                                                                    Types::scalar basis_fn_module) const {
    // Проверка, что делится нацело
    if ((mesh.nx() - 1) % Nx != 0 || (mesh.ny() - 1) % Ny != 0 || (mesh.nz() - 1) % Nz != 0) {
        throw std::invalid_argument("OperatorK::compute_galerkin_matrix_custom_blocksize: "
            "Nx, Ny, Nz is not consistent with mesh size dimentions");
    }
    const Idx3d sizes{Nx, Ny, Nz};
    const size_t first_layer_toeplitz = mesh.nCubesX() / Nx;
    const size_t second_layer_toeplitz = mesh.nCubesY() / Ny;
    const size_t third_layer_toeplitz = mesh.nCubesZ() / Nz;
    const size_t inner_size = 3 * Nx * Ny * Nz;

    decltype(auto) result =
        Math::LinAgl::Matrix::ZeroTripleToeplitzBlock<Types::complex_d>(first_layer_toeplitz, second_layer_toeplitz,
                                                                        third_layer_toeplitz, inner_size);

#pragma omp parallel for num_threads(14) shared(result) firstprivate(Nx, Ny, Nz, basis_fn_module)
    for (size_t j3 = 0; j3 < third_layer_toeplitz; ++j3) {
        // цикл по первой строке в матрице
        size_t i3 = 0;
        auto &&working_block_on_tl = result.get_block(i3, j3);

        for (size_t i2 = 0; i2 < second_layer_toeplitz; ++i2)
            for (size_t i1 = 0; i1 < first_layer_toeplitz; ++i1)
                for (size_t j2 = 0; j2 < second_layer_toeplitz; ++j2)
                    for (size_t j1 = 0; j1 < first_layer_toeplitz; ++j1) {
                        auto &&working_block = working_block_on_tl.
                                               get_block(i2, j2).
                                               get_block(i1, j1);
                        // Ускорение заполнения матрицы за счет отсутствия
                        // пересчёта одинаковых блоков
                        // TODO: сделать нормальный расчет, то есть аналитически вывести все формулки
                        // TODO: тогда тут будет 3 цикла вместо 6
                        if (working_block.norm() == 0) {
                            // Расчет триджы-тёплицевой матрицы для соответствующих коллекций кубов
                            // (в плотном формате) и запись в соответствующий блок большой матрицы
                            const Idx3d start_i = {i1 * sizes.Nx, i2 * sizes.Ny, i3 * sizes.Nz};
                            const Idx3d start_j = {j1 * sizes.Nx, j2 * sizes.Ny, j3 * sizes.Nz};
                            working_block = compute_galerkin_matrix(start_i, start_j, sizes, basis_fn_module)
                                .to_dense();
                            // Из самого забавного: тут получается 12 вложенных циклов for.
                            // Что-то мне не очень это нравится.
                        }
                    }
    }

#pragma omp parallel for num_threads(14) shared(result) firstprivate(Nx, Ny, Nz, basis_fn_module)
    for (size_t i3 = 1; i3 < third_layer_toeplitz; ++i3) {
        // цикл по первому столбцу в матрице
        size_t j3 = 0;
        auto &&working_block_on_tl = result.get_block(i3, j3);

        for (size_t i2 = 0; i2 < second_layer_toeplitz; ++i2)
            for (size_t i1 = 0; i1 < first_layer_toeplitz; ++i1)
                for (size_t j2 = 0; j2 < second_layer_toeplitz; ++j2)
                    for (size_t j1 = 0; j1 < first_layer_toeplitz; ++j1) {
                        auto &&working_block = working_block_on_tl.
                                               get_block(i2, j2).
                                               get_block(i1, j1);
                        // Ускорение заполнения матрицы за счет отсутствия
                        // пересчёта одинаковых блоков
                        // TODO: сделать нормальный расчет, то есть аналитически вывести все формулки
                        // TODO: тогда тут будет 3 цикла вместо 6
                        if (working_block.norm() == 0) {
                            // Расчет триджы-тёплицевой матрицы для соответствующих коллекций кубов
                            // (в плотном формате) и запись в соответствующий блок большой матрицы
                            const Idx3d start_i = {i1 * sizes.Nx, i2 * sizes.Ny, i3 * sizes.Nz};
                            const Idx3d start_j = {j1 * sizes.Nx, j2 * sizes.Ny, j3 * sizes.Nz};
                            working_block = compute_galerkin_matrix(start_i, start_j, sizes, basis_fn_module)
                                .to_dense();
                            // Из самого забавного: тут получается 12 вложенных циклов for.
                            // Что-то мне не очень это нравится.
                        }
                    }
    }
    return {result, mesh.getPermutation(Nx, Ny, Nz)};
}

operator_K_over_cube_mesh::matrix_and_permutation<Math::LinAgl::Matrix::TripleToeplitzFactoredBlock<Types::complex_d>>
operator_K_over_cube_mesh::compute_galerkin_matrix_custom_blocksize_compressed(size_t Nx, size_t Ny, size_t Nz,
                                                                               Types::scalar basis_fn_module,
                                                                               Types::scalar epsilon) const {
    // Проверка, что делится нацело
    if ((mesh.nx() - 1) % Nx != 0 || (mesh.ny() - 1) % Ny != 0 || (mesh.nz() - 1) % Nz != 0) {
        throw std::invalid_argument("OperatorK::compute_galerkin_matrix_custom_blocksize: "
            "Nx, Ny, Nz is not consistent with mesh size dimentions");
    }
    const Idx3d sizes{Nx, Ny, Nz};
    const size_t first_layer_toeplitz = mesh.nCubesX() / Nx;
    const size_t second_layer_toeplitz = mesh.nCubesY() / Ny;
    const size_t third_layer_toeplitz = mesh.nCubesZ() / Nz;
    const size_t inner_size = 3 * Nx * Ny * Nz;

    decltype(auto) result =
        Math::LinAgl::Matrix::ZeroTripleToeplitzFactoredBlock<Types::complex_d>(
            first_layer_toeplitz, second_layer_toeplitz,
            third_layer_toeplitz, inner_size);

    Types::scalar norm_of_self_interation_block = 1;

#pragma omp parallel for num_threads(14) shared(result) firstprivate(Nx, Ny, Nz, basis_fn_module)
    for (size_t j3 = 0; j3 < third_layer_toeplitz; ++j3) {
        // цикл по первой строке в матрице
        size_t i3 = 0;
        auto &&working_block_on_tl = result.get_block(i3, j3);

        for (size_t i2 = 0; i2 < second_layer_toeplitz; ++i2)
            for (size_t i1 = 0; i1 < first_layer_toeplitz; ++i1)
                for (size_t j2 = 0; j2 < second_layer_toeplitz; ++j2)
                    for (size_t j1 = 0; j1 < first_layer_toeplitz; ++j1) {
                        auto &&working_block = working_block_on_tl.
                                               get_block(i2, j2).
                                               get_block(i1, j1);
                        // Ускорение заполнения матрицы за счет отсутствия
                        // пересчёта одинаковых блоков
                        // TODO: сделать нормальный расчет, то есть аналитически вывести все формулки
                        // TODO: тогда тут будет 3 цикла вместо 5
                        if (working_block.factor_number() == 0) {
                            // Расчет триджы-тёплицевой матрицы для соответствующих коллекций кубов
                            // (в плотном формате) и запись в соответствующий блок большой матрицы
                            const Idx3d start_i = {i1 * sizes.Nx, i2 * sizes.Ny, i3 * sizes.Nz};
                            const Idx3d start_j = {j1 * sizes.Nx, j2 * sizes.Ny, j3 * sizes.Nz};
                            auto dense_block = compute_galerkin_matrix(start_i, start_j, sizes, basis_fn_module)
                                .to_dense();
                            if (start_i == start_j) {
                                working_block = Math::LinAgl::Matrix::DynamicFactoredMatrix<decltype(dense_block)>{
                                    {std::move(dense_block)}};
                                norm_of_self_interation_block = working_block.get<0>().norm();
                            } else {
                                // Для начала подкрутим точность относительно диагонального
                                const auto local_epsilon =
                                    norm_of_self_interation_block / dense_block.norm() * epsilon;
                                // Теперь делаем всё для креста
                                const auto row_fun = [&dense_block](Types::index m)-> Types::VectorXc {
                                    return dense_block.row(m);
                                };
                                const auto col_fun = [&dense_block](Types::index m)-> Types::VectorXc {
                                    return dense_block.col(m);
                                };
                                working_block = Math::LinAgl::Decompositions::ComplexACA::svd_postcompression(
                                    Math::LinAgl::Decompositions::ComplexACA::compute(
                                        row_fun, col_fun, dense_block.rows(), dense_block.cols(), local_epsilon),
                                    local_epsilon);

                                // Для пущей важности можно посмотреть на ранги
                                // std::cout << "rank = " << working_block.get<0>().cols() << '\n';
                            }
                        }

                    }
    }

#pragma omp parallel for num_threads(14) shared(result) firstprivate(Nx, Ny, Nz, basis_fn_module)
    for (size_t i3 = 1; i3 < third_layer_toeplitz; ++i3) {
        // цикл по первой строке в матрице
        size_t j3 = 0;
        auto &&working_block_on_tl = result.get_block(i3, j3);

        for (size_t i2 = 0; i2 < second_layer_toeplitz; ++i2)
            for (size_t i1 = 0; i1 < first_layer_toeplitz; ++i1)
                for (size_t j2 = 0; j2 < second_layer_toeplitz; ++j2)
                    for (size_t j1 = 0; j1 < first_layer_toeplitz; ++j1) {
                        auto &&working_block = working_block_on_tl.
                                               get_block(i2, j2).
                                               get_block(i1, j1);
                        // Ускорение заполнения матрицы за счет отсутствия
                        // пересчёта одинаковых блоков
                        // TODO: сделать нормальный расчет, то есть аналитически вывести все формулки
                        // TODO: тогда тут будет 3 цикла вместо 5
                        if (working_block.factor_number() == 0) {
                            // Расчет триджы-тёплицевой матрицы для соответствующих коллекций кубов
                            // (в плотном формате) и запись в соответствующий блок большой матрицы
                            const Idx3d start_i = {i1 * sizes.Nx, i2 * sizes.Ny, i3 * sizes.Nz};
                            const Idx3d start_j = {j1 * sizes.Nx, j2 * sizes.Ny, j3 * sizes.Nz};
                            auto dense_block = compute_galerkin_matrix(start_i, start_j, sizes, basis_fn_module)
                                .to_dense();
                            if (start_i == start_j) {
                                working_block = Math::LinAgl::Matrix::DynamicFactoredMatrix<decltype(dense_block)>{
                                    {std::move(dense_block)}};
                                norm_of_self_interation_block = working_block.get<0>().norm();
                            } else {
                                // Для начала подкрутим точность относительно диагонального
                                const auto local_epsilon =
                                    norm_of_self_interation_block / dense_block.norm() * epsilon;
                                // Теперь делаем всё для креста
                                const auto row_fun = [&dense_block](Types::index m)-> Types::VectorXc {
                                    return dense_block.row(m);
                                };
                                const auto col_fun = [&dense_block](Types::index m)-> Types::VectorXc {
                                    return dense_block.col(m);
                                };
                                working_block = Math::LinAgl::Decompositions::ComplexACA::svd_postcompression(
                                    Math::LinAgl::Decompositions::ComplexACA::compute(
                                        row_fun, col_fun, dense_block.rows(), dense_block.cols(), local_epsilon),
                                    local_epsilon);

                                // Для пущей важности можно посмотреть на ранги
                                // std::cout << "rank = " << working_block.get<0>().cols() << '\n';
                            }
                        }
                    }

    }
    return {result, mesh.getPermutation(Nx, Ny, Nz)};
}


operator_K_over_cube_mesh::matrix_and_permutation<Math::LinAgl::Matrix::TripleToeplitzFactoredBlock<Types::complex_d>>
operator_K_over_cube_mesh::compute_galerkin_matrix_custom_blocksize_compressed(size_t Nx, size_t Ny, size_t Nz,
                                                                               Types::scalar basis_fn_module,
                                                                               Types::scalar epsilon,
                                                                               Math::LinAgl::Matrix::TripleToeplitzBlock
                                                                               <Types::complex_d> &dense_mat) const {
    // Проверка, что делится нацело
    if ((mesh.nx() - 1) % Nx != 0 || (mesh.ny() - 1) % Ny != 0 || (mesh.nz() - 1) % Nz != 0) {
        throw std::invalid_argument("OperatorK::compute_galerkin_matrix_custom_blocksize: "
            "Nx, Ny, Nz is not consistent with mesh size dimentions");
    }
    const Idx3d sizes{Nx, Ny, Nz};
    const size_t first_layer_toeplitz = mesh.nCubesX() / Nx;
    const size_t second_layer_toeplitz = mesh.nCubesY() / Ny;
    const size_t third_layer_toeplitz = mesh.nCubesZ() / Nz;
    const size_t inner_size = 3 * Nx * Ny * Nz;

    decltype(auto) result =
        Math::LinAgl::Matrix::ZeroTripleToeplitzFactoredBlock<Types::complex_d>(
            first_layer_toeplitz, second_layer_toeplitz,
            third_layer_toeplitz, inner_size);

    dense_mat = Math::LinAgl::Matrix::ZeroTripleToeplitzBlock<Types::complex_d>(
        first_layer_toeplitz, second_layer_toeplitz,
        third_layer_toeplitz, inner_size);

    Types::scalar norm_of_self_interation_block = 1;

#pragma omp parallel for num_threads(14) shared(result, dense_mat) firstprivate(Nx, Ny, Nz, basis_fn_module)
    for (size_t j3 = 0; j3 < third_layer_toeplitz; ++j3) {
        // цикл по первой строке в матрице
        size_t i3 = 0;
        auto &&working_block_on_tl = result.get_block(i3, j3);

        for (size_t i2 = 0; i2 < second_layer_toeplitz; ++i2)
            for (size_t i1 = 0; i1 < first_layer_toeplitz; ++i1)
                for (size_t j2 = 0; j2 < second_layer_toeplitz; ++j2)
                    for (size_t j1 = 0; j1 < first_layer_toeplitz; ++j1) {
                        auto &&working_block = working_block_on_tl.get_block(i2, j2).
                                                                   get_block(i1, j1);
                        auto &&dense_mat_block = dense_mat.get_block(i3, j3).
                                                           get_block(i2, j2).
                                                           get_block(i1, j1);
                        // Ускорение заполнения матрицы за счет отсутствия
                        // пересчёта одинаковых блоков
                        // TODO: сделать нормальный расчет, то есть аналитически вывести все формулки
                        // TODO: тогда тут будет 3 цикла вместо 6
                        if (working_block.factor_number() == 0) {
                            // Расчет триджы-тёплицевой матрицы для соответствующих коллекций кубов
                            // (в плотном формате) и запись в соответствующий блок большой матрицы
                            const Idx3d start_i = {i1 * sizes.Nx, i2 * sizes.Ny, i3 * sizes.Nz};
                            const Idx3d start_j = {j1 * sizes.Nx, j2 * sizes.Ny, j3 * sizes.Nz};
                            dense_mat_block = compute_galerkin_matrix(start_i, start_j, sizes, basis_fn_module)
                                .to_dense();
                            if (start_i == start_j) {
                                working_block = Math::LinAgl::Matrix::DynamicFactoredMatrix<Types::MatrixXc>{
                                    {dense_mat_block}};
                                norm_of_self_interation_block = working_block.get<0>().norm();
                            } else {
                                // Для начала подкрутим точность относительно диагонального
                                const auto local_epsilon = epsilon;
                                // Теперь делаем всё для креста
                                const auto row_fun = [&dense_mat_block](Types::index m)-> Types::VectorXc {
                                    return dense_mat_block.row(m);
                                };
                                const auto col_fun = [&dense_mat_block](Types::index m)-> Types::VectorXc {
                                    return dense_mat_block.col(m);
                                };
                                working_block = Math::LinAgl::Decompositions::ComplexACA::svd_postcompression(
                                    Math::LinAgl::Decompositions::ComplexACA::compute(
                                        row_fun, col_fun, dense_mat_block.rows(), dense_mat_block.cols(),
                                        local_epsilon),
                                    local_epsilon);

                                // Для пущей важности можно посмотреть на ранги
                                // std::cout << "rank = " << working_block.get<0>().cols() << '\n';
                            }
                        }
                    }
    }

#pragma omp parallel for num_threads(14) shared(result) firstprivate(Nx, Ny, Nz, basis_fn_module)
    for (size_t i3 = 1; i3 < third_layer_toeplitz; ++i3) {
        // цикл по первой строке в матрице
        size_t j3 = 0;
        auto &&working_block_on_tl = result.get_block(i3, j3);

        for (size_t i2 = 0; i2 < second_layer_toeplitz; ++i2)
            for (size_t i1 = 0; i1 < first_layer_toeplitz; ++i1)
                for (size_t j2 = 0; j2 < second_layer_toeplitz; ++j2)
                    for (size_t j1 = 0; j1 < first_layer_toeplitz; ++j1) {
                        auto &&working_block = working_block_on_tl.get_block(i2, j2).
                                                                   get_block(i1, j1);
                        auto &&dense_mat_block = dense_mat.get_block(i3, j3).
                                                           get_block(i2, j2).
                                                           get_block(i1, j1);
                        // Ускорение заполнения матрицы за счет отсутствия
                        // пересчёта одинаковых блоков
                        // TODO: сделать нормальный расчет, то есть аналитически вывести все формулки
                        // TODO: тогда тут будет 3 цикла вместо 6
                        if (working_block.factor_number() == 0) {
                            // Расчет триджы-тёплицевой матрицы для соответствующих коллекций кубов
                            // (в плотном формате) и запись в соответствующий блок большой матрицы
                            const Idx3d start_i = {i1 * sizes.Nx, i2 * sizes.Ny, i3 * sizes.Nz};
                            const Idx3d start_j = {j1 * sizes.Nx, j2 * sizes.Ny, j3 * sizes.Nz};
                            dense_mat_block = compute_galerkin_matrix(start_i, start_j, sizes, basis_fn_module)
                                .to_dense();
                            if (start_i == start_j) {
                                working_block = Math::LinAgl::Matrix::DynamicFactoredMatrix<Types::MatrixXc>{
                                    {dense_mat_block}};
                                norm_of_self_interation_block = working_block.get<0>().norm();
                            } else {
                                // Для начала подкрутим точность относительно диагонального
                                const auto local_epsilon =
                                    norm_of_self_interation_block / dense_mat_block.norm() * epsilon;
                                // Теперь делаем всё для креста
                                const auto row_fun = [&dense_mat_block](Types::index m)-> Types::VectorXc {
                                    return dense_mat_block.row(m);
                                };
                                const auto col_fun = [&dense_mat_block](Types::index m)-> Types::VectorXc {
                                    return dense_mat_block.col(m);
                                };
                                working_block = Math::LinAgl::Decompositions::ComplexACA::svd_postcompression(
                                    Math::LinAgl::Decompositions::ComplexACA::compute(
                                        row_fun, col_fun, dense_mat_block.rows(), dense_mat_block.cols(),
                                        local_epsilon),
                                    local_epsilon);

                                // Для пущей важности можно посмотреть на ранги
                                // std::cout << "rank = " << working_block.get<0>().cols() << '\n';
                            }
                        }
                    }
    }
    return {result, mesh.getPermutation(Nx, Ny, Nz)};
}


}
