//
// Created by evgen on 13.10.2024.
//

#include "mesh/Utils.hpp"
#include "Parser.hpp"

namespace EMW::Mesh::Utils {

Mesh::SurfaceMesh loadTriangularMesh(int nNodes, int nCells, std::string type) {
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "plate/triangulated/" +
                                  type + "/nodes/" + std::to_string(nNodes) + "_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "plate/triangulated/" +
                                  type + "/cells/" + std::to_string(nCells) + "_cells.csv";

    auto preMesh = EMW::Parser::parseMesh(nodesFile, cellsFile);
    auto surfaceMesh = EMW::Mesh::SurfaceMesh(
        preMesh.first, preMesh.second,
        [&](const Containers::array<Types::index, 4> &p, const Containers::vector<point_t> &fp) -> point_t {
            return (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(3)) * (fp[p[0]] + fp[p[1]] + fp[p[2]]);
        });
    surfaceMesh.setName(type + "_triangular_mesh_" + std::to_string(nCells));
    return surfaceMesh;
}

Mesh::SurfaceMesh move_by_vector(const Mesh::SurfaceMesh &mesh, const Types::Vector3d &v) {
    // собираем новые узлы сетки
    const auto new_nodes_view = mesh.getNodes() | std::views::transform([&](auto &node) {
        return Types::Vector3d{node + v};
    });
    const Containers::vector<point_t> nodes{new_nodes_view.begin(), new_nodes_view.end()};
    // собираем старое их объединение в ячейки
    const auto new_cells_view = mesh.getCells() | std::views::transform([](auto &cell) {return cell.points_;});
    const Containers::vector<IndexedCell::nodes_t> cells{new_cells_view.begin(), new_cells_view.end()};
    return {nodes, cells};
}

} // namespace EMW::Mesh::Utils
