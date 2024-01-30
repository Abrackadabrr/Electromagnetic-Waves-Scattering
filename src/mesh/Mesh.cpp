//
// Created by evgen on 30.01.24.
//

#include "Mesh.hpp"

namespace EMW::Mesh {
    SurfaceMesh::SurfaceMesh(Containers::vector<Point> nodes, Containers::vector<IndexedCell> cells,
                             Containers::vector<Node> collocationPoints) : nodes_(nodes), cells_(cells) {

    };
}