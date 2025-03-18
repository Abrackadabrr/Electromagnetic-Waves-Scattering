//
// Created by evgen on 25.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_SURFACEMESH_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_SURFACEMESH_HPP

#include "MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::Mesh {
/**
 * Класс расчетной поверхностной сетки, дискретный аналог поверхности
 */
class SurfaceMesh {

    struct mesh_info_t {
       Types::index nodes_size;
       Types::index cells_size;
       std::string name;
    };

    Containers::vector<point_t> nodes_;
    Containers::vector<IndexedCell> cells_;
    std::string name = "default_mesh_name";

    SurfaceMesh(const Containers::vector<point_t> &nodes, const Containers::vector<IndexedCell> &cells)
        : nodes_(nodes), cells_(cells) {}

  public:
    /**
     * Дефолт-конструирование пустой поверхности
     * Когда все пустое, то поверхность просто описывает пустое множество, так что состояние консистентно
     */
    SurfaceMesh() = default;

    /**
     * Конструирование сетки по узлам и индексным ячейкам
     */
    SurfaceMesh(Containers::vector<point_t> nodes, Containers::vector<IndexedCell::nodes_t> cells);

    /**
     * Конструирование поверхностной сетки со специальными точками коллокации
     */
    SurfaceMesh(Containers::vector<point_t> nodes, Containers::vector<Containers::array<Types::index, 4>> cells,
                const std::function<point_t(const Containers::array<Types::index, 4> &,
                                            const Containers::vector<point_t> &)> &getPoint);

    /**
     * Задать пользовательские локальные базисы
     */
    void customLocalBasis(const std::function<std::array<Types::Vector3d, 3>(const Mesh::IndexedCell &)> &func);

    /**
     * Задать пользовательские точки коллокации
     */
    void customCollocationPpoints(const std::function<Types::Vector3d(const Mesh::IndexedCell &)> &func);

    /* Далее идут селекторы */
    [[nodiscard]] constexpr const Containers::vector<IndexedCell> &getCells() const { return cells_; }

    [[nodiscard]] constexpr const Containers::vector<point_t> &getNodes() const { return nodes_; }

    [[nodiscard]] std::string getName() const { return name; }

    void setName(const std::string &n) {name = n;}

    /* Немного странные костыльные методы */

    /**
     * Вернуть часть сетки, которая помечана специальным тэгом
     */
    SurfaceMesh getSubmesh(IndexedCell::Tag tag) const;

    mesh_info_t getSubmeshInfo(IndexedCell::Tag tag) const;

    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_SURFACEMESH_HPP
