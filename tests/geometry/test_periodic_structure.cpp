//
// Created by evgen on 18.01.2025.
//

#include "VTKFunctions.hpp"

#include "gtest/gtest.h"
#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/Utils.hpp"

#include "mesh/Parser.hpp"

#include "geometry/PeriodicStructure.hpp"
#include "geometry/ShiftedPeriodicStructure.hpp"

using namespace EMW;

int main() {
    // считываем сетку на антенне
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/8000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/2000_cells.csv";
    constexpr EMW::Types::index nNodes = 8000;
    constexpr EMW::Types::index nCells = 2000;

    // собираем сетки
    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};

    constexpr Types::index N1 = 3;
    constexpr Types::index N2 = 3;
    constexpr Types::index N1_x_N2 = N1 * N2;

    const Types::scalar a = 0.07;  // расстояние между центрами сеток "на диагонали 1"
    const Types::scalar b = 0.14;  // расстояние между центрами сеток "на диагонали 2"
    const Types::scalar step = std::sqrt(a * a + b * b / 4);

    const Types::scalar alpha = std::atan2(b, 2 * a);
    const Types::scalar y_coord_of_vector = (2 * a / b);
    const Types::Vector3d dir1 = Types::Vector3d{1, y_coord_of_vector, 0}.normalized();
    const Types::Vector3d dir2 = Types::Vector3d{1, -y_coord_of_vector, 0}.normalized();

    const Geometry::ShiftedStructure<N1, N2> geometry{dir1, dir2, step, step, mesh_base};

    const auto new_g = geometry.expand_without_saving_nice_origin();

    VTK::geometry_snapshot(new_g, "/home/evgen/Education/MasterDegree/thesis/"
                               "Electromagnetic-Waves-Scattering/tests/geometry/shifted_lattice.vtu");

}
