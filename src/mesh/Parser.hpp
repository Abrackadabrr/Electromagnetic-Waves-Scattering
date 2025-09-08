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

struct ParserOutput {
    Containers::vector<Mesh::point_t> nodes;
    Containers::vector<Containers::array<Types::index, 4>> cells;
    Containers::vector<std::string> tags;

    void shrink_to_fit() {
        nodes.shrink_to_fit();
        cells.shrink_to_fit();
        tags.shrink_to_fit();
    }
};

ParserOutput parseMesh(const std::string &csvNodes, const std::string &csvCells);
ParserOutput parse_mesh_without_tag(const std::string &csvNodes, const std::string &csvCells);

} // namespace EMW::Parser

#endif // ELECTROMAGNETIC_WAVES_SCATTERING_OBJPARSER_H
