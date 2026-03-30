//
// Created by evgen on 23.03.2026.
//

#include <gtest/gtest.h>

#include "types/Types.hpp"

#include "math/matrix/Matrix.hpp"

#include <mesh/volume_mesh/CubeMeshWithData.hpp>
#include <operators/volume/OperatorK.hpp>

#include <random>

using namespace EMW;

std::pair<Types::index, Types::index> full_index_from_block_index(Types::index i1, Types::index j1,
                                                                  Types::index i2, Types::index j2,
                                                                  Types::index i3, Types::index j3,
                                                                  Types::index fl, Types::index sl, Types::index tl) {
    Types::index internal_size = 1;
    Types::index fl_stride = internal_size;
    Types::index sl_stride = fl * fl_stride;
    Types::index tl_stride = sl * sl_stride;

    Types::index row = i3 * tl_stride + i2 * sl_stride + i1;
    Types::index col = j3 * tl_stride + j2 * sl_stride + j1;

    return {row, col};
}

Types::scalar scalar_toeplitz_matrix_element(Types::index i1, Types::index j1,
                                             Types::index i2, Types::index j2, Types::index i3, Types::index j3) {
    return (Types::scalar(i1) - Types::scalar(j1) + 3) *
           (Types::scalar(i2) - Types::scalar(j2) + 3) *
           (Types::scalar(i3) - Types::scalar(j3) + 3);
}

decltype(auto) scalar_triple_toeplitz_matrix(Types::index fl, Types::index sl, Types::index tl) {
    Types::index total_size = fl * sl * tl;

    auto full_matrix = Math::LinAgl::Matrix::ZeroTripleToeplitzBlock<Types::scalar>(fl, sl, tl, 1);
    for (size_t i3 = 0; i3 < tl; i3++)
        for (size_t j3 = 0; j3 < tl; j3++)
            for (size_t i2 = 0; i2 < sl; i2++)
                for (size_t j2 = 0; j2 < sl; j2++)
                    for (size_t i1 = 0; i1 < fl; i1++)
                        for (size_t j1 = 0; j1 < fl; j1++) {
                            full_matrix.get_block(i3, j3).get_block(i2, j2).get_block(i1, j1)(0, 0) =
                                scalar_toeplitz_matrix_element(i1, j1, i2, j2, i3, j3);
                        }
    return full_matrix;
}

Types::MatrixXd full_scalar_triple_toeplitz_matrix(Types::index fl, Types::index sl, Types::index tl) {
    Types::index total_size = fl * sl * tl;
    Types::index fl_stride = 1;
    Types::index sl_stride = fl * fl_stride;
    Types::index tl_stride = sl * sl_stride;

    Types::MatrixXd full_matrix = Types::MatrixXd::Zero(total_size, total_size);
    for (size_t i = 0; i < total_size; i++)
        for (size_t j = 0; j < total_size; j++) {
            auto i3 = i / tl_stride;
            auto j3 = j / tl_stride;
            auto i2 = (i % tl_stride) / sl_stride;
            auto j2 = (j % tl_stride) / sl_stride;
            auto i1 = ((i % tl_stride) % sl_stride) / fl_stride;
            auto j1 = ((j % tl_stride) % sl_stride) / fl_stride;
            const auto idx = full_index_from_block_index(i1, j1, i2, j2, i3, j3, fl, sl, tl);
            if (!(idx.first == i && idx.second == j)) {
                std::cout << i << " " << j << " " << idx.first << " " << idx.second << std::endl;
            }
            full_matrix(i, j) = scalar_toeplitz_matrix_element(i1, j1, i2, j2, i3, j3);
        }
    return full_matrix;
}


TEST(TOEPLITZ_MATRIX_TESTS, TRIPLE_TOEPLITZ_SCALAR) {
    // Описываем структуру матрицы
    const Types::index first_layer_rows = 30;
    const Types::index first_layer_cols = first_layer_rows;

    const Types::index second_layer_rows = 12;
    const Types::index second_layer_cols = second_layer_rows;

    const Types::index third_layer_rows = 3;
    const Types::index third_layer_cols = third_layer_rows;

    // Строим обычную матрицу, без специального хранения элементов
    const Types::MatrixXd full_matrix =
        full_scalar_triple_toeplitz_matrix(first_layer_rows, second_layer_rows, third_layer_rows);

    // Собираем дважды тёплицеву матрицу
    const auto trip_toep_matrix = scalar_triple_toeplitz_matrix(first_layer_rows, second_layer_rows, third_layer_rows);

    // Проверка совпадения матриц
    ASSERT_NEAR((full_matrix - trip_toep_matrix.to_dense()).norm(), 0, 1e-14);
}

constexpr Types::scalar SPHERE_RADUIS = 0.5;

Types::scalar permittivity_distribution(const Types::point_t &x) {
    return x.norm() < SPHERE_RADUIS ? 2.56 : 1;
}

TEST(TOEPLITZ_MATRIX_TESTS, VIE_DIFFERENT_ENUMRATION_TEST) {
    constexpr Types::scalar cube_length = 2.1 * SPHERE_RADUIS;
    constexpr Types::index Nx = 13;
    constexpr Types::index Ny = 13;
    constexpr Types::index Nz = 13;
    constexpr Types::scalar mesh_one_axis_size = cube_length / (Nx - 1);
    constexpr Types::scalar basis_fn_norm = 1. / (mesh_one_axis_size * std::sqrt(mesh_one_axis_size));
    Mesh::VolumeMesh::CubeMeshWithData mesh{Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2},
                                            (Nx - 1) * mesh_one_axis_size,
                                            (Ny - 1) * mesh_one_axis_size, (Nz - 1) * mesh_one_axis_size, Nx, Ny, Nz};

    Operators::Volume::operator_K_over_cube_mesh operator_K{{1., 0.}, mesh};
    auto toeplitz = operator_K.compute_galerkin_matrix(basis_fn_norm);

    auto [mat, perm] = operator_K.compute_galerkin_matrix_custom_blocksize(
        4, 4, 4, basis_fn_norm);

    const Types::VectorXc x = Types::VectorXc::Random(toeplitz.cols());

    // Плотная умножалка
    const Types::VectorXc x_permuted = perm * x;
    const Types::VectorXc mat_result = perm.transpose() * mat.matvec(x_permuted);
    const auto toeplitz_result = toeplitz.matvec(x);

    ASSERT_NEAR((mat_result - toeplitz_result).norm(), 0.0, 1e-6 * toeplitz_result.norm());
}

TEST(TOEPLITZ_MATRIX_TESTS, VIE_SKELTONIZATION_TEST) {
    constexpr Types::scalar cube_length = 2.1 * SPHERE_RADUIS;
    constexpr Types::index Nx = 13;
    constexpr Types::index Ny = 13;
    constexpr Types::index Nz = 13;
    constexpr Types::scalar mesh_one_axis_size = cube_length / (Nx - 1);
    constexpr Types::scalar basis_fn_norm = 1. / (mesh_one_axis_size * std::sqrt(mesh_one_axis_size));
    Mesh::VolumeMesh::CubeMeshWithData mesh{Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2},
                                            (Nx - 1) * mesh_one_axis_size,
                                            (Ny - 1) * mesh_one_axis_size, (Nz - 1) * mesh_one_axis_size, Nx, Ny, Nz};

    Operators::Volume::operator_K_over_cube_mesh operator_K{{1., 0.}, mesh};
    size_t nx, ny, nz;
    nx = 4;
    ny = 4;
    nz = 4;
    auto [mat_dense, perm_dense] = operator_K.compute_galerkin_matrix_custom_blocksize(
        nx, ny, nz, basis_fn_norm);
    auto [mat, perm] = operator_K.compute_galerkin_matrix_custom_blocksize_compressed(
        nx, ny, nz, basis_fn_norm, 1e-3);

    std::cout << "Approximation error = " << (mat.to_dense() - mat_dense.to_dense()).norm() / mat_dense.to_dense().norm() << std::endl;
}
