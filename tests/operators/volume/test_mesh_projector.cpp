//
// Created by evgen on 27.01.2026.
//
#include <gtest/gtest.h>

#include "operators/volume/ProjectorOnMesh.hpp"

#include "mesh/volume_mesh/CubeMeshWithData.hpp"

using namespace EMW;

TEST(PROJECTOR_ON_MESH, CONSTANT_PROJECTION) {
    // Параметры куба
    constexpr Types::scalar cube_length = 0.025;
    constexpr Types::index Nx = 3;
    constexpr Types::index Ny = static_cast<size_t>(1 / cube_length) + 1;
    constexpr Types::index Nz = Ny;

    // Сетка
    Mesh::VolumeMesh::CubeMeshWithData mesh{Types::point_t{0, 0, 0}, (Nx - 1) * cube_length,
                                            (Ny - 1) * cube_length, (Nz - 1) * cube_length, Nx, Ny, Nz};

    // Проектор
    Operators::Volume::ProjectorOnMesh proj{mesh};
    auto rhs = proj([](Types::point_t p) {
            return Types::Vector3c{std::complex{0., 0.}, std::complex{0., 0.}, std::complex{1., 0.}};
        });

    // Проверка правильности расчета
    const auto cube_measure = mesh.dx() * mesh.dy() * mesh.dz();
    for (Types::index i = 0; i < mesh.getCells().size(); i++) {
        ASSERT_NEAR(rhs[3 * i].real(), 0, 1e-15);
        ASSERT_NEAR(rhs[3 * i + 1].real(), 0, 1e-15);
        ASSERT_NEAR(rhs[3 * i + 2].real(), cube_measure, 1e-15);
    }
}
