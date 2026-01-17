//
// Created by evgen on 17.01.2026.
//

#ifndef MESHBASE_HPP
#define MESHBASE_HPP

#include "types/Types.hpp"

namespace EMW::Mesh {

template<typename cell_t>
class MeshBase {

    struct mesh_info_t {
        Types::index nodes_size;
        Types::index cells_size;
        std::string name;
    };

    Containers::vector<Types::point_t> nodes_;
    Containers::vector<cell_t> cells_;
    std::string name = "default_mesh_name";
};

}
#endif //MESHBASE_HPP
