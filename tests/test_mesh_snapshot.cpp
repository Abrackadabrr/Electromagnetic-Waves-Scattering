//
// Created by evgen on 30.01.24.
//

#include "gtest/gtest.h"
#include "Types.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/MeshTypes.hpp"
#include "visualisation/VTKFunctions.hpp"

using namespace EMW;

TEST(MESH, VTK_SNAPSHOT_TEST) {
    Containers::vector<Mesh::Point> points{
            Mesh::Point{0, 0, 0},
            Mesh::Point{1, 0, 0},
            Mesh::Point{2, 0, 0},
            Mesh::Point{0, 1, 0},
            Mesh::Point{1, 1, 0},
            Mesh::Point{2, 1, 0},
            Mesh::Point{0, 2, 0},
            Mesh::Point{1, 2, 0},
            Mesh::Point{2, 2, 0}
    };

    Containers::vector<Mesh::IndexedCell::nodes_t> cells{
            Mesh::IndexedCell::nodes_t{0, 1, 4, 3},
            Mesh::IndexedCell::nodes_t{1, 2, 5, 4},
            Mesh::IndexedCell::nodes_t{3, 4, 7, 6},
            Mesh::IndexedCell::nodes_t{4, 5, 8, 7}
    };

    const Mesh::SurfaceMesh mesh{points, cells};
    VTK::test_snapshot(0, mesh, "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/");
}
