//
// Created by evgen on 17.01.2026.
//

#include <gtest/gtest.h>

#include "operators/volume/OperatorK.hpp"
#include "operators/volume/ProjectorOnMesh.hpp"

#include "mesh/volume_mesh/CubeMesh.hpp"

#include "experiment/PhysicalCondition.hpp"

#include <Utils.hpp>


using namespace EMW;

class VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS : public ::testing::Test {
protected:
    Types::scalar cube_length{};
    Types::complex_d k{};

    void SetUp() override {
        cube_length = 0.1;
        k = {2.0 * Math::Constants::PI<Types::scalar>() / cube_length, 0};
    }
};

TEST_F(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, SingularIntegrationTest) {
    constexpr size_t N = 40;
    cube_length = 1;
    const Types::complex_d low_freq_k = {6 * Math::Constants::PI<Types::scalar>(), 0.};
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2},
                                    cube_length, N};
    Operators::Volume::operator_K_over_cube_mesh op_K{low_freq_k, mesh};
    // Сравнение расчета через сингулярное интегрирование и через интегрирование в дальней зоне
    const Types::index idx1 = 0, idx2 = mesh.getCells().size() - 1;
    std::cout << "distance = " << mesh.distance(idx1, idx2) << std::endl;
    const auto sing_res = op_K.galerkin_block_for_cubes(idx1, idx2);
    const auto far_res = op_K.far_zone_interaction(idx1, idx2);

    const Types::scalar rel_err = std::abs((sing_res - far_res).norm() / far_res.norm());
    std::cout << rel_err << std::endl;
}

TEST_F(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, SimpleTripleBlockToeplitzTest) {
    constexpr size_t Nx2 = 2;
    constexpr size_t Nx3 = Nx2 + 1;
    constexpr size_t internal_n = 3;
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{0, 0, 0}, cube_length, Nx3};
    Operators::Volume::operator_K_over_cube_mesh op_K{k, mesh};
    const auto galerkin_matrix = op_K.compute_galerkin_matrix_dense(1);
    // Банальное соотвествие размерам
    {
        ASSERT_EQ(galerkin_matrix.rows(), Nx2 * Nx2 * Nx2 * internal_n);
        ASSERT_EQ(galerkin_matrix.cols(), Nx2 * Nx2 * Nx2 * internal_n);
        std::cout << "Rows: " << galerkin_matrix.rows() << std::endl;
        std::cout << "Cols: " << galerkin_matrix.cols() << std::endl;
    }
    // Проверка локального расчета для очевидных пары кубов (0 -> 1, 2 -> 3, 4 -> 5, 6 -> 7, и отдельно 0 -> 4, 3 -> 7)
    const auto res12_1 = op_K.matrix_3_coef(0, 1);
    const auto res23_1 = op_K.matrix_3_coef(2, 3);
    const auto res45_1 = op_K.matrix_3_coef(4, 5);
    const auto res67_1 = op_K.matrix_3_coef(6, 7);
    ASSERT_NEAR(std::abs(res12_1 - res23_1), 0, 1e-14 * std::abs(res12_1));
    ASSERT_NEAR(std::abs(res12_1 - res45_1), 0, 1e-14 * std::abs(res12_1));
    ASSERT_NEAR(std::abs(res12_1 - res67_1), 0, 1e-14 * std::abs(res12_1));
    const auto res04_1 = op_K.matrix_3_coef(0, 4);
    const auto res37_1 = op_K.matrix_3_coef(3, 7);
    ASSERT_NEAR(std::abs(res04_1 - res37_1), 0, 1e-14 * std::abs(res37_1));

    // Для кубической сетки 2х2х2 блочная тёплицевость проверяется только на диагональных блоках
    // Верхний уровень тёплицевости
    {
        const auto block1 = galerkin_matrix.block(0, 0, Nx2 * Nx2 * internal_n, Nx2 * Nx2 * internal_n);
        const auto block2 = galerkin_matrix.block(Nx2 * Nx2 * internal_n, Nx2 * Nx2 * internal_n,
                                                  Nx2 * Nx2 * internal_n, Nx2 * Nx2 * internal_n);
        const Types::scalar absErr = (block1 - block2).norm();
        ASSERT_NEAR(absErr, 0, 1e-15 * (block1).norm());
    }
    // Средний уровень тёплицевости
    {
        const auto block1 = galerkin_matrix.block(0, 0, Nx2 * internal_n, Nx2 * internal_n);
        const auto block2 = galerkin_matrix.block(Nx2 * internal_n, Nx2 * internal_n, Nx2 * internal_n,
                                                  Nx2 * internal_n);
        const auto block3 = galerkin_matrix.block(2 * Nx2 * internal_n, 2 * Nx2 * internal_n, Nx2 * internal_n,
                                                  Nx2 * internal_n);
        const auto block4 = galerkin_matrix.block(3 * Nx2 * internal_n, 3 * Nx2 * internal_n, Nx2 * internal_n,
                                                  Nx2 * internal_n);
        ASSERT_NEAR((block1 - block2).norm(), 0, 2e-15 * block1.norm());
        ASSERT_NEAR((block1 - block3).norm(), 0, 2e-15 * block1.norm());
        ASSERT_NEAR((block1 - block4).norm(), 0, 2e-15 * block1.norm());
    }
    // Внутренний уровень тёплицевости
    {
        const auto block_ref = galerkin_matrix.block(0, 0, internal_n, internal_n);
        for (size_t idx = 0; idx < Nx2 * Nx2 * Nx2; idx++) {
            const auto block_to_compare = galerkin_matrix.block(idx * internal_n, idx * internal_n, internal_n,
                                                                internal_n);
            ASSERT_NEAR((block_ref - block_to_compare).norm(), 0, 1e-14 * block_ref.norm());
        }
    }
}

TEST_F(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, EQUALITY_OF_MATRIX_ELEMENTS) {
    constexpr Types::scalar total_mesh_size = 2;
    constexpr Types::index Nx = 40;
    constexpr Types::index Ncubes = Nx - 1;
    constexpr Types::scalar cube_size = total_mesh_size / (Nx - 1);
    constexpr Types::scalar basis_fn_module = 1. / sqrt(cube_size * cube_size * cube_size);
    constexpr Types::scalar rel_tol = 5e-13;

    // Берем кубическую сетку на кубе
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{0, 0, 0}, total_mesh_size, Nx};
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};

    // Цикл по парам кубов, для которых элементы в матрице должны быть одинаковы
    for (size_t h_idx = 0; h_idx < Ncubes; h_idx++) {
        for (size_t col_idx = 1; col_idx < Ncubes; col_idx++) {
            auto cube_1 = mesh.cube_idx(0, 0, h_idx);
            auto cube_2 = mesh.cube_idx(0, col_idx, h_idx);
            const auto reference_block = operator_K.galerkin_block_for_cubes(cube_1, cube_2);
            // Теперь спускаемся вниз и считаем блоки для кубов, которые в точной арифметике имеют reference_block
            for (size_t row_idx = 0; row_idx < Ncubes; row_idx++) {
                cube_1 = mesh.cube_idx(row_idx, 0, h_idx);
                cube_2 = mesh.cube_idx(row_idx, col_idx, h_idx);
                const auto current_block = operator_K.galerkin_block_for_cubes(cube_1, cube_2);
                ASSERT_NEAR((current_block - reference_block).norm(), 0,
                            rel_tol * reference_block.norm()) << col_idx << ' ' << row_idx;
            }
        }
    }

    // Ещё одни пары кубов с одинаковыми блоками
    {
        auto cube_1 = mesh.cube_idx(0, 0, 0);
        auto cube_2 = mesh.cube_idx(0, Ncubes - 1, Ncubes - 1);
        const auto reference_block = operator_K.galerkin_block_for_cubes(cube_1, cube_2);
        // Двигаем кубы по оси X
        for (size_t idx = 0; idx < Ncubes; idx++) {
            cube_1 = mesh.cube_idx(idx, 0, 0);
            cube_2 = mesh.cube_idx(idx, Ncubes - 1, Ncubes - 1);
            const auto current_block = operator_K.galerkin_block_for_cubes(cube_1, cube_2);
            ASSERT_NEAR((current_block - reference_block).norm(), 0,
                        rel_tol * reference_block.norm()) << idx;
        }
    }

    // И ещё одни пары кубов с одинаковыми блоками
    {
        for (size_t jdx = 0; jdx < Ncubes; jdx++) {
            auto cube_1 = mesh.cube_idx(0, 0, 0);
            auto cube_2 = mesh.cube_idx(jdx, Ncubes - 1, Ncubes - 1);
            const auto reference_block = operator_K.galerkin_block_for_cubes(cube_1, cube_2);
            // Двигаем кубы по оси X
            for (size_t idx = 0; idx < Ncubes - jdx; idx++) {
                cube_1 = mesh.cube_idx(idx, 0, 0);
                cube_2 = mesh.cube_idx(idx + jdx, Ncubes - 1, Ncubes - 1);
                const auto current_block = operator_K.galerkin_block_for_cubes(cube_1, cube_2);
                ASSERT_NEAR((current_block - reference_block).norm(), 0,
                            rel_tol * reference_block.norm()) << idx;
            }
        }
    }
}

TEST_F(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, TOEPLITZ_DENSE_COMPARISON) {
    constexpr Types::scalar total_mesh_size = 0.1;
    constexpr Types::index Nx = 5;
    constexpr Types::scalar cube_size = total_mesh_size / (Nx - 1);
    constexpr Types::scalar basis_fn_module = 1. / sqrt(cube_size * cube_size * cube_size);
    constexpr Types::scalar rel_tol = 1e-14;

    // берем кубическую сетку на кубе
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{0, 0, 0}, cube_length, Nx};
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};
    // собираем полную матрицу
    auto full_mat = operator_K.compute_galerkin_matrix_dense(basis_fn_module);
    std::cout << "Matrix (" << full_mat.rows() << ", " << full_mat.cols() << ')' << std::endl;

    // Собираем тёплицеву матрицу
    auto toep_mat = operator_K.compute_galerkin_matrix(basis_fn_module);
    std::cout << "Toep Matrix (" << toep_mat.rows() << ", " << toep_mat.cols() << ')' << std::endl;
    std::cout << full_mat.norm() << std::endl;

    // Сравниваем типов
    ASSERT_NEAR((full_mat - toep_mat.to_dense()).norm(), 0, rel_tol * full_mat.norm());
}

TEST_F(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, TOEPLITZ_TOEPLITZ_COMPARISON) {
    constexpr Types::scalar total_mesh_size = 0.1;
    constexpr Types::index Nx = 17;
    constexpr Types::scalar cube_size = total_mesh_size / (Nx - 1);
    constexpr Types::scalar basis_fn_module = 1. / sqrt(cube_size * cube_size * cube_size);
    constexpr Types::scalar rel_tol = 2e-14;

    // берем кубическую сетку на кубе
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{0, 0, 0}, cube_length, Nx};
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};

    const auto toeplitz_structure = [](
        const Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d> &mat) -> std::vector<size_t> {
        return {mat.rows_in_toeplitrz(),
                mat.get_block(0, 0).rows_in_toeplitrz(),
                mat.get_block(0, 0).get_block(0, 0).rows_in_toeplitrz(),
                static_cast<size_t>(mat.get_block(0, 0).get_block(0, 0).get_block(0, 0).rows())};
    };

    // Собираем тёплицеву матрицу
    auto [toep_mat, _] = operator_K.compute_galerkin_matrix_custom_blocksize(1, 1, 1, basis_fn_module);
    auto toep_mat_str = toeplitz_structure(toep_mat);
    std::cout << "Toep Matrix (" << toep_mat.rows() << ", " << toep_mat.cols() << ')' << std::endl;
    std::cout << "3: " << toep_mat_str[0] << "\n2: " << toep_mat_str[1] << "\n1: " << toep_mat_str[2] <<
        "\ninner: " << toep_mat_str[3] << std::endl;

    // Собираем тёплицеву матрицу, но с блоками побольше
    auto [big_toep_mat, perm] = operator_K.compute_galerkin_matrix_custom_blocksize(4, 2, 2, basis_fn_module);
    toep_mat_str = toeplitz_structure(big_toep_mat);
    std::cout << "Big Matrix (" << big_toep_mat.rows() << ", " << big_toep_mat.cols() << ')' << std::endl;
    std::cout << "3: " << toep_mat_str[0] << "\n2: " << toep_mat_str[1] << "\n1: " << toep_mat_str[2] <<
        "\ninner: " << toep_mat_str[3] << std::endl;

    // Сравниваем обе матрицы c учетом перестановки строк и столбцов у первой
    ASSERT_NEAR((perm.transpose() * big_toep_mat.to_dense() * perm - toep_mat.to_dense()).norm(), 0,
                rel_tol * big_toep_mat.to_dense().norm());
}

TEST_F(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, TOEPLITZ_COMPARISON_WITH_ACA_COMPRESSED) {
    Eigen::setNbThreads(1);
    openblas_set_num_threads(1);
    constexpr Types::scalar total_mesh_size = 0.1;
    constexpr Types::index Nx = 17;
    constexpr Types::scalar cube_size = total_mesh_size / (Nx - 1);
    constexpr Types::scalar basis_fn_module = 1. / sqrt(cube_size * cube_size * cube_size);
    constexpr Types::scalar rel_tol = 2e-14;

    // берем кубическую сетку на кубе
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{0, 0, 0}, cube_length, Nx};
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};
    std::cout << "Computing toeplitz matrix ..." << std::endl;
    auto toep_mat = operator_K.compute_galerkin_matrix_custom_blocksize(4, 4, 4, basis_fn_module);
    std::cout << "End. Computing compressed matrix ..." << std::endl;
    auto compressed_mat = operator_K.
        compute_galerkin_matrix_custom_blocksize_compressed(4, 4, 4, basis_fn_module, 1e-5);
    std::cout << "End. Computing error" << std::endl;
    const auto full_toep = toep_mat.matrix.to_dense();
    const auto full_compressed = compressed_mat.matrix.to_dense();

    std::cout << (full_toep - full_compressed).norm() / full_toep.norm() << std::endl;
}

TEST_F(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, TOEPLITZ_MATVEC_COMPARISON_WITH_ACA_COMPRESSED) {
    Eigen::setNbThreads(1);
    openblas_set_num_threads(1);
    constexpr Types::scalar total_mesh_size = 0.1;
    constexpr Types::index Nx = 17;
    constexpr Types::scalar cube_size = total_mesh_size / (Nx - 1);
    constexpr Types::scalar basis_fn_module = 1. / sqrt(cube_size * cube_size * cube_size);
    constexpr Types::scalar rel_tol = 2e-14;

    // берем кубическую сетку на кубе
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{0, 0, 0}, cube_length, Nx};
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};
    std::cout << "Computing toeplitz matrix ..." << std::endl;
    auto toep_mat = operator_K.compute_galerkin_matrix_custom_blocksize(4, 4, 4, basis_fn_module);
    std::cout << "End. Computing compressed matrix ..." << std::endl;
    auto compressed_mat = operator_K.
        compute_galerkin_matrix_custom_blocksize_compressed(4, 4, 4, basis_fn_module, 1e-5);
    std::cout << Utils::get_memory_usage(compressed_mat.matrix) << std::endl;
    std::cout << compressed_mat.matrix.cols() << " " <<
        Utils::get_elements_for_parametrization(compressed_mat.matrix) / compressed_mat.matrix.cols() << std::endl;
    std::cout << "End. Computing error" << std::endl;

    Types::VectorXc x = Types::VectorXc::Random(toep_mat.matrix.cols());

    const auto full_toep = toep_mat.matrix * x;;
    const auto full_compressed = compressed_mat.matrix * x;

    std::cout << (full_toep - full_compressed).norm() / full_toep.norm() << std::endl;

    std::cout << compressed_mat.matrix.get_block(0, 2).get_block(0, 2).get_block(0, 2).to_dense().norm() << std::endl;
}
