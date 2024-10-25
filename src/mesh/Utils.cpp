//
// Created by evgen on 13.10.2024.
//

#include "mesh/Utils.hpp"
#include "Parser.hpp"

namespace EMW::Mesh::Utils {

Mesh::SurfaceMesh loadTriangularMesh(int nNodes, int nCells, std::string type) {
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "plate/triangulated/" + type + "/nodes/" +
                                  std::to_string(nNodes) + "_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "plate/triangulated/" + type + "/cells/" +
                                  std::to_string(nCells) + "_cells.csv";

    auto surfaceMesh = Mesh::SurfaceMesh{EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells)};
    surfaceMesh.setName(type + "_triangular_mesh_" + std::to_string(nCells));
    return surfaceMesh;
}
} // namespace EMW::Mesh::Utils