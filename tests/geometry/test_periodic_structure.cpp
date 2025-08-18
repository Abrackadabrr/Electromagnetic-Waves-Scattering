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

#include <geometry/HomogeneousStructure.hpp>

using namespace EMW;

int main() {
    // считываем сетку на антенне
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/8000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/2000_cells.csv";
    // собираем сетки
    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};

    constexpr Types::index N1 = 4;  // количество строк в решетке
    constexpr Types::index N2 = 4;  // количество столбцов в решетке
    constexpr Types::index N1_x_N2 = N1 * N2;
    // out waveguides
    Containers::set<Types::index> out_waveguides{2, 3, 7, 8, 12, 13};

    const Types::scalar a_hat = 0.07;  // расстояние между центрами сеток "на диагонали 1"
    const Types::scalar b_hat = 0.14;  // расстояние между центрами сеток "на диагонали 2"
    const Types::scalar step = std::sqrt(a_hat * a_hat + b_hat * b_hat / 4);

    const Types::scalar alpha = std::atan2(b_hat, 2 * a_hat);
    const Types::scalar y_coord_of_vector = (2 * a_hat / b_hat);
    const Types::Vector3d dir1 = Types::Vector3d{1, -1, 0}.normalized();
    const Types::Vector3d dir2 = Types::Vector3d{1, 1, 0}.normalized();

    Geometry::PeriodicStructure<N1, N2> geometry{dir1, dir2, step, step, mesh_base};

    Containers::vector<Types::Vector3d> origins;

    for (int i = 0; i < N1_x_N2; i++)
        if (!out_waveguides.contains(i))
            origins.push_back(geometry.get_origin(i));

    const Geometry::HomogeneousStructure new_geom{origins, mesh_base};

    VTK::geometry_snapshot(new_geom, "/home/evgen/Education/MasterDegree/thesis/"
                               "Electromagnetic-Waves-Scattering/tests/geometry/shifted_hex_lattice.vtu");
}
