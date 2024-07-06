//
// Created by evgen on 20.02.24.
//

#include <vector>
#include <utility>
#include <tuple>
#include <ranges>
#include <fstream>

#include "gtest/gtest.h"
#include "types/Types.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/MeshTypes.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "mesh/Parser.hpp"

using namespace EMW;

TEST(MESH, VTK_SNAPSHOT_TEST) {
    const std::string nodes = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/2002_nodes.csv";
    const std::string cells = "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/examples/2050_cells.csv";

    Mesh::SurfaceMesh mesh = EMW::Parser::parseMesh(nodes, cells, 2002, 2050);
    mesh.setName("mesh_test_parsing");
    VTK::surface_snapshot(0, mesh, "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/");
}
