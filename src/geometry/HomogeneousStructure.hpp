//
// Created by evgen on 18.01.2025.
//

#ifndef HOMOGENEOUSSTRUCTURE_HPP
#define HOMOGENEOUSSTRUCTURE_HPP

#include "Scene.hpp"

#include "types/Types.hpp"

namespace EMW::Geometry {
class HomogeneousStructure : Scene {
    Containers::vector<Types::Vector3d> origins;
};
} // namespace EMW::Geometry

#endif // HOMOGENEOUSSTRUCTURE_HPP
