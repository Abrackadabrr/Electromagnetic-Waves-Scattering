//
// Created by evgen on 14.10.2024.
//

#include "Parser.hpp"
#include "third_party/csv/csv.h"
namespace EMW::Parser {

ParserOutput parseMesh(const std::string &csvNodes, const std::string &csvCells) {
    ParserOutput output;

    io::CSVReader<3> nodesF(csvNodes);
    io::CSVReader<5> cellsF(csvCells);

    nodesF.read_header(io::ignore_extra_column, "x", "y", "z");
    cellsF.read_header(io::ignore_extra_column, "f", "s", "t", "fou", "group");
    Types::scalar x, y, z;
    Types::index f, s, t, fou;
    std::string tag;
    while (nodesF.read_row(x, y, z)) {
        output.nodes.emplace_back(x, y, z);
    }

    while (cellsF.read_row(f, s, t, fou, tag)) {
        output.cells.push_back({f, s, t, fou});
        output.tags.push_back(tag);
    }
    output.shrink_to_fit();
    return output;
}

ParserOutput parse_mesh_without_tag(const std::string &csvNodes, const std::string &csvCells) {
    ParserOutput output;

    io::CSVReader<3> nodesF(csvNodes);
    io::CSVReader<4> cellsF(csvCells);

    nodesF.read_header(io::ignore_extra_column, "x", "y", "z");
    cellsF.read_header(io::ignore_extra_column, "f", "s", "t", "fou");
    Types::scalar x, y, z;
    Types::index f, s, t, fou;
    while (nodesF.read_row(x, y, z)) {
        output.nodes.emplace_back(x, y, z);
    }

    while (cellsF.read_row(f, s, t, fou)) {
        output.cells.push_back({f, s, t, fou});
    }
    output.shrink_to_fit();
    return output;
}

} // namespace EMW::Parser
