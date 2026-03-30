//
// Created by evgen on 17.01.2026.
//

#include "VTKFunctions.hpp"

#include <gtest/gtest.h>

#include "types/Types.hpp"

#include "mesh/volume_mesh/CubeMesh.hpp"
#include "mesh/volume_mesh/VolumeCells.hpp"

using namespace EMW;

TEST(VOLUME_MESH_TESTS, GET_FACE_ON_CUBE) {
    const Containers::vector<Types::point_t> full_points{
        Types::point_t{0, 0, 0}, Types::point_t{1, 0, 0}, Types::point_t{0, 1, 0}, Types::point_t{1, 1, 0},
        Types::point_t{0, 0, 1}, Types::point_t{1, 0, 1}, Types::point_t{0, 1, 1}, Types::point_t{1, 1, 1},
    };

    const Containers::array<Types::index, 8> idx = {0, 1, 2, 3, 4, 5, 6, 7};

    Mesh::VolumeCells::IndexedCube test_cell{full_points, idx};

    const auto xface = test_cell.getFace(Mesh::VolumeCells::IndexedCube::Axis::Z,
                                         Mesh::VolumeCells::IndexedCube::Direction::Plus, full_points);

    std::cout << xface.normal << std::endl;
    const auto vertecies = xface.getVertex();
    std::cout << vertecies.d << std::endl;
}

TEST(VOLUME_MESH_TESTS, SIMPLE_CREATION_AND_SNAPSHOT) {
    Types::point_t corner{0, 0, 0};
    const Types::scalar xs = 1;
    const Mesh::VolumeMesh::CubeMeshWithData mesh{corner, xs, 3};

    VTK::volume_mesh_withdata_snapshot(
        mesh, "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/tests/mesh/volume_mesh/");
}

TEST(VOLUME_MESH_TESTS, PERMUTATION_CALCULATION) {
    Types::point_t corner{0, 0, 0};
    const Types::scalar xs = 1;
    constexpr size_t N = 3;
    const Mesh::VolumeMesh::CubeMeshWithData mesh{corner, xs, N};

    const auto p1 = mesh.getPermutationForCubes(1, 2, 2);
    std::cout << p1.indices() << '\n' << std::endl;
    Eigen::VectorXi vec = Eigen::VectorXi::Ones(p1.indices().size());
    for (size_t i = 0; i < p1.indices().size(); i++) {
        vec[i] = i;
    }
    std::cout << vec <<"\n"<<std::endl;
    std::cout << p1 * vec << std::endl;
}
