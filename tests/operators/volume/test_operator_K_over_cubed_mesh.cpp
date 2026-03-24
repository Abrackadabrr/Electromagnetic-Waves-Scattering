//
// Created by evgen on 17.01.2026.
//

#include <gtest/gtest.h>

#include "operators/volume/OperatorK.hpp"
#include "operators/volume/ProjectorOnMesh.hpp"

#include "mesh/volume_mesh/CubeMesh.hpp"

#include "experiment/PhysicalCondition.hpp"


using namespace EMW;

class VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS : public ::testing::Test {
protected:
    Types::scalar cube_length{};
    size_t Nx5{};
    Types::complex_d k{};

    void SetUp() override {
        cube_length = 0.1;
        Nx5 = 5;
        k = {2.0 * Math::Constants::PI<Types::scalar>() / cube_length, 0};
    }
};

TEST_F(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, SimpleTripleBlockToeplitzTest) {
    constexpr size_t Nx2 = 2;
    constexpr size_t Nx3 = Nx2 + 1;
    constexpr size_t internal_n = 3;
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{0, 0, 0}, cube_length, Nx3};
    Operators::Volume::operator_K_over_cube_mesh op_K{k, mesh};
    const auto galerkin_matrix = op_K.get_galerkin_matrix();
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
        const Types::scalar relErr = (block1 - block2).norm() / (block1).norm();
        ASSERT_NEAR(relErr, 0, 1e-15);
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

TEST_F(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, TEST_SOLVING) {
    constexpr Types::scalar cube_length = 0.1;
    constexpr Types::index Nx = 8;
    // берем кубическую сетку на кубе
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{0, 0, 0}, cube_length, Nx};
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};
    // собираем полную матрицу
    auto full_mat = operator_K.get_galerkin_matrix();
    std::cout << "Matrix (" << full_mat.rows() << ", " << full_mat.cols() << ')' << std::endl;

    // Собираем тёплицеву матрицу
    auto toep_mat = operator_K.compute_galerkin_matrix();
    std::cout << "Toep Matrix (" << toep_mat.rows() << ", " << toep_mat.cols() << ')' << std::endl;
    std::cout << full_mat.norm() << std::endl;
    // Сравниваем типов
    ASSERT_NEAR((full_mat - toep_mat.to_dense()).norm(), 0, 1e-14 * full_mat.norm());
}
