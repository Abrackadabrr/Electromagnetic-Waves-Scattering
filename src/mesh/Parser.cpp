//
// Created by evgen on 14.10.2024.
//

#include "Parser.hpp"
#include "third_party/csv/csv.h"
namespace EMW::Parser {
std::pair<Containers::vector<Mesh::point_t>, Containers::vector<Containers::array<Types::index, 4>>>
parseMesh(const std::string &csvNodes, const std::string &csvCells, int nNodes, int nCells) {
    Containers::vector<Mesh::point_t> nodes;
    nodes.reserve(nNodes);
    Containers::vector<Containers::array<Types::index, 4>> cells;
    cells.reserve(nCells);

    io::CSVReader<3> nodesF(csvNodes);
    io::CSVReader<4> cellsF(csvCells);

    nodesF.read_header(io::ignore_extra_column, "x", "y", "z");
    cellsF.read_header(io::ignore_extra_column, "f", "s", "t", "fou");
    Types::scalar x, y, z;
    Types::index f, s, t, fou;

    while (nodesF.read_row(x, y, z)) {
        nodes.emplace_back(x, y, z);
    }

    while (cellsF.read_row(f, s, t, fou)) {
        cells.push_back({f, s, t, fou});
    }
    nodes.shrink_to_fit();
    cells.shrink_to_fit();
    return {nodes, cells};
}
} // namespace EMW::Parser