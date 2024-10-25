//
// Created by evgen on 14.10.2024.
//

#include "Parser.hpp"

EMW::Mesh::SurfaceMesh EMW::Parser::parseMesh(const std::string & csvNodes, const std::string & csvCells, int nNodes, int nCells) {
    Containers::vector<Mesh::point_t> nodes;
    nodes.reserve(nNodes);
    Containers::vector<Containers::array<Types::index, 4>> cells;
    cells.reserve(nCells);

    Containers::vector<Mesh::Node::field_t> testField; testField.reserve(nCells);

    io::CSVReader<3> nodesF(csvNodes);
    io::CSVReader<4> cellsF(csvCells);

    nodesF.read_header(io::ignore_extra_column, "x", "y", "z");
    cellsF.read_header(io::ignore_extra_column, "f", "s", "t", "fou");
    Types::scalar x, y, z;
    Types::index f, s, t, fou;

    while(nodesF.read_row(x, y, z)){
        nodes.emplace_back(x, y, z);
    }

    int counter = 0;
    while(cellsF.read_row(f, s, t, fou)){
        cells.push_back({f-1, s-1, t-1, fou-1});
        testField.emplace_back(Types::complex_d{static_cast<double>(counter), 0}, Types::complex_d{0, 0}, Types::complex_d{0, 0});
        counter++;
    }

    return Mesh::SurfaceMesh{nodes, cells};
}