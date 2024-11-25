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

    auto preMesh = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
    auto surfaceMesh = EMW::Mesh::SurfaceMesh(
        preMesh.first, preMesh.second,
        [&](const Containers::array<Types::index, 4> &p, const Containers::vector<point_t> &fp) -> point_t {
            return (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(3)) * (fp[p[0]] + fp[p[1]] + fp[p[2]]);
        });
    surfaceMesh.setName(type + "_triangular_mesh_" + std::to_string(nCells));
    return surfaceMesh;
}
} // namespace EMW::Mesh::Utils