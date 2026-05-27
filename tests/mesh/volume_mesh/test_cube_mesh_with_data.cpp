//
// Created by evgen on 10.02.2026.
//

#include "mesh/volume_mesh/CubeMeshWithData.hpp"

#include "visualisation/include/VTKFunctions.hpp"

#include <gtest/gtest.h>

using namespace EMW;

void AssertCubeMeshWithDataEqual(const Mesh::VolumeMesh::CubeMeshWithData& lhs,
                                 const Mesh::VolumeMesh::CubeMeshWithData& rhs, const double tol = 1e-12) {
    const auto& lhs_nodes = lhs.getNodes();
    const auto& rhs_nodes = rhs.getNodes();
    ASSERT_EQ(lhs_nodes.size(), rhs_nodes.size()) << "Nodes count mismatch";
    for (Types::index i = 0; i < lhs_nodes.size(); ++i) {
        for (int comp = 0; comp < 3; ++comp) {
            ASSERT_NEAR(lhs_nodes[i][comp], rhs_nodes[i][comp], tol) << "Node mismatch at idx=" << i;
        }
    }

    const auto& lhs_cells = lhs.getCells();
    const auto& rhs_cells = rhs.getCells();
    ASSERT_EQ(lhs_cells.size(), rhs_cells.size()) << "Cells count mismatch";
    for (Types::index i = 0; i < lhs_cells.size(); ++i) {
        for (int k = 0; k < 8; ++k) {
            ASSERT_EQ(lhs_cells[i].nodes_[k], rhs_cells[i].nodes_[k]) << "Cell connectivity mismatch at idx=" << i;
        }
    }

    const auto& lhs_scalar_data = lhs.getScalarData();
    const auto& rhs_scalar_data = rhs.getScalarData();
    ASSERT_EQ(lhs_scalar_data.size(), rhs_scalar_data.size()) << "Scalar fields count mismatch";
    for (const auto& [name, lhs_data] : lhs_scalar_data) {
        ASSERT_TRUE(rhs_scalar_data.contains(name)) << "Missing scalar field: " << name;
        const auto& rhs_data = rhs_scalar_data.find(name)->second;
        ASSERT_EQ(lhs_data.size(), rhs_data.size()) << "Scalar field size mismatch for field: " << name;
        for (Types::index i = 0; i < lhs_data.size(); ++i) {
            ASSERT_NEAR(lhs_data[i].real(), rhs_data[i].real(), tol) << "Scalar real mismatch in field: " << name;
            ASSERT_NEAR(lhs_data[i].imag(), rhs_data[i].imag(), tol) << "Scalar imag mismatch in field: " << name;
        }
    }

    const auto& lhs_vector_data = lhs.getVectorData();
    const auto& rhs_vector_data = rhs.getVectorData();
    ASSERT_EQ(lhs_vector_data.size(), rhs_vector_data.size()) << "Vector fields count mismatch";
    for (const auto& [name, lhs_data] : lhs_vector_data) {
        ASSERT_TRUE(rhs_vector_data.contains(name)) << "Missing vector field: " << name;
        const auto& rhs_data = rhs_vector_data.find(name)->second;
        ASSERT_EQ(lhs_data.size(), rhs_data.size()) << "Vector field size mismatch for field: " << name;
        for (Types::index i = 0; i < lhs_data.size(); ++i) {
            for (int comp = 0; comp < 3; ++comp) {
                ASSERT_NEAR(lhs_data[i][comp].real(), rhs_data[i][comp].real(), tol)
                    << "Vector real mismatch in field: " << name;
                ASSERT_NEAR(lhs_data[i][comp].imag(), rhs_data[i][comp].imag(), tol)
                    << "Vector imag mismatch in field: " << name;
            }
        }
    }
}

TEST(CUBE_MESH_WITH_DATA, DUMP_TEST) {
    Types::point_t corner{0, 0, 0};
    const Types::scalar xs = 1;
    Mesh::VolumeMesh::CubeMeshWithData mesh{corner, xs, 10};

    // добавляем данные на сетку с данными
    const auto epsilon = [](Types::point_t x) { return Types::complex_d{x.squaredNorm(), 0}; };
    mesh.invokeScalarData("eps", epsilon);

    VTK::volume_mesh_withdata_snapshot(
        mesh, "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/tests/mesh/volume_mesh/");
}

TEST(CUBE_MESH_WITH_DATA, READ_TEST) {

    const std::string path =
        "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/tests/mesh/volume_mesh/";
    const std::string mesh_name = "dumped_mesh";
    Types::point_t corner{0, 0, 0};
    const Types::scalar xs = 1;
    Mesh::VolumeMesh::CubeMeshWithData mesh{corner, xs, 100};
    mesh.setName(mesh_name);

    // добавляем данные на сетку с данными
    const auto epsilon = [](Types::point_t x) { return Types::complex_d{x.squaredNorm(), 0}; };
    mesh.invokeScalarData("eps", epsilon);

    VTK::volume_mesh_withdata_snapshot(mesh, path);

    auto read_mesh = VTK::volume_mesh_withdata_from_vtu(path + mesh_name + ".vtu");

    // Mesh equivelance
    AssertCubeMeshWithDataEqual(read_mesh, mesh);

    read_mesh.setName("redumped_mesh");
    VTK::volume_mesh_withdata_snapshot(read_mesh, path);
}

TEST(CUBE_MESH_WITH_DATA, CHECK_CORRECT_ORDER) {
    Types::point_t corner{0, 0, 0};
    const Types::scalar xs = 1;
    Mesh::VolumeMesh::CubeMeshWithData mesh{corner, xs, 10};

    const Types::index idx = mesh.cube_idx(5, 6, 7);

    // добавляем данные на сетку с данными
    const auto epsilon = [](Types::point_t x) { return Types::complex_d{x.squaredNorm(), 0}; };
    mesh.invokeScalarData("eps", epsilon);

    const auto epsvec = mesh.getScalarDataAsVector("eps");

    for (Types::index i = 0; i < mesh.getCells().size(); i++) {
        ASSERT_NEAR(epsvec[i].real(), epsilon(mesh.getCells()[i].center_).real(), 1e-310);
    }
}
