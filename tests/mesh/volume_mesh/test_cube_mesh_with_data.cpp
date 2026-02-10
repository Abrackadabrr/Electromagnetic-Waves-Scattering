//
// Created by evgen on 10.02.2026.
//

#include "mesh/volume_mesh/CubeMeshWithData.hpp"

#include "visualisation/VTKFunctions.hpp"

#include <gtest/gtest.h>

using namespace EMW;

TEST(CUBE_MESH_WITH_DATA, DUMP_TEST) {
    Types::point_t corner{0, 0, 0};
    const Types::scalar xs = 1;
    Mesh::VolumeMesh::CubeMeshWithData mesh{corner, xs, 10};

    // добавляем данные на сетку с данными
    const auto epsilon = [](Types::point_t x) { return Types::complex_d{x.squaredNorm(), 0};};
    mesh.invokeScalarData("eps", epsilon);

    VTK::volume_mesh_withdata_snapshot(
        mesh, "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/tests/mesh/volume_mesh/");
}

TEST(CUBE_MESH_WITH_DATA, CHECK_CORRECT_ORDER) {
    Types::point_t corner{0, 0, 0};
    const Types::scalar xs = 1;
    Mesh::VolumeMesh::CubeMeshWithData mesh{corner, xs, 10};

    const Types::index idx = mesh.cube_idx(5, 6, 7);

    // добавляем данные на сетку с данными
    const auto epsilon = [](Types::point_t x) { return Types::complex_d{x.squaredNorm(), 0};};
    mesh.invokeScalarData("eps", epsilon);

    const auto epsvec = mesh.getScalarDataAsVector("eps");

    for (Types::index i = 0; i < mesh.getCells().size(); i++) {
        ASSERT_NEAR(epsvec[i].real(), epsilon(mesh.getCells()[i].center_).real(), 1e-310);
    }
}
