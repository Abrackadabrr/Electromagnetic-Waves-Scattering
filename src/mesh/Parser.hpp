//
// Created by evgen on 25.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_OBJPARSER_H
#define ELECTROMAGNETIC_WAVES_SCATTERING_OBJPARSER_H

#include "../types/Types.hpp"
#include "MeshTypes.hpp"
#include "SurfaceMesh.hpp"
#include <string>

namespace EMW::Parser {
std::pair<Containers::vector<Mesh::point_t>, Containers::vector<Containers::array<Types::index, 4>>>
parseMesh(const std::string &csvNodes, const std::string &csvCells, int nNodes, int nCells);

}

#endif // ELECTROMAGNETIC_WAVES_SCATTERING_OBJPARSER_H
