//
// Created by evgen on 17.01.2026.
//

#include <gtest/gtest.h>

#include "operators/volume/OperatorK.hpp"
#include "operators/volume/ProjectorOnMesh.hpp"

#include "mesh/volume_mesh/CubeMesh.hpp"

#include "experiment/PhysicalCondition.hpp"


using namespace EMW;

TEST(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, COMPILATION_TEST) {
    constexpr Types::scalar cube_length = 1.0;
    constexpr Types::complex_d k{2.0 * Math::Constants::PI<Types::scalar>() / cube_length, 0};
    constexpr Types::index Nx = 10;
    // берем кубическую сетку на кубе
    std::cout << "Nx = " << Nx << std::endl;
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{0, 0, 0}, cube_length, Nx};
    // строим матрицу оператора методом галеркина
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};
    std::cout << "Nx = " << Nx << std::endl;
    const auto mat = operator_K.get_galerkin_matrix();
    std::cout << mat.rows() << ' ' << mat.cols() << std::endl;
}

TEST(VOLUME_OPERATOR_OVER_CUBE_MESH_TESTS, TEST_SOLVING) {
    constexpr Types::scalar cube_length = 1.0;
    constexpr Types::complex_d k{2.0 * Math::Constants::PI<Types::scalar>() / cube_length, 0};
    constexpr Types::index Nx = 5;
    // берем кубическую сетку на кубе
    Mesh::VolumeMesh::CubeMesh mesh{Types::point_t{0, 0, 0}, cube_length, Nx};
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};
    // собираем матрицу
    auto mat = operator_K.get_galerkin_matrix();
    std::cout << "Matrix (" << mat.rows() << ", " << mat.cols() << ')' << std::endl;
    mat -= Types::MatrixXc::Identity(mat.rows(), mat.cols());
    // собираем правую часть
    Physics::planeWaveCase incident_field{Types::Vector3d{1, 0, 0}, k, Types::Vector3d{0, 1, 0}};
    Operators::Volume::ProjectorOnMesh proj{mesh};
    const auto rhs = proj([incident_field](Types::point_t p) { return incident_field.value(p); });
    std::cout << "RHS size: " << rhs.rows() << std::endl;
}
