//
// Created by evgen on 05.02.2025.
//

#ifndef CONSEPTS_HPP
#define CONSEPTS_HPP

#include <type_traits>

#include "types/Types.hpp"

namespace EMW::Consepts {
template<typename Topology>
concept manifold_like = requires(Topology topology) {
    {topology.getCells()} -> Containers::vector<typename Topology::CellType>;
};
}

#endif //CONSEPTS_HPP
