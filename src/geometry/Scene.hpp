//
// Created by evgen on 18.01.2025.
//

#ifndef SCENE_HPP
#define SCENE_HPP

#include "types/Types.hpp"

#include "mesh/SurfaceMesh.hpp"

#include <memory>
#include <utility>

namespace EMW::Geometry {
/**
 * Класс, описывающий геометрию сцены, которую нужно будет рассчитывать
 * Потребность в таком классе заключается в моделировании периодических структур
 * Scene -- это базовый класс. Далее от него будет отнаследовано что-то более кнокретное
 */
class Scene {
    using array_of_surfaces = Containers::vector<std::unique_ptr<Mesh::SurfaceMesh>>;
    array_of_surfaces surfaces;
public:
    Scene() = default;
    explicit Scene(array_of_surfaces  array) : surfaces(std::move(array)) {}
};
}

#endif //SCENE_HPP
