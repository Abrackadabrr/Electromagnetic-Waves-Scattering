//
// Created by evgen on 17.01.2026.
//

#ifndef VOLUMECELLS_HPP
#define VOLUMECELLS_HPP

#include "../MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::Mesh::VolumeCells {

namespace Utils {
enum Axis { X = 0, Y = 1, Z = 2 };
enum Direction { Minus = 0, Plus = 1 };

Mesh::IndexedCell getFace(Axis ax, Direction dir, const Containers::vector<point_t> &fp);

Mesh::IndexedCell newGetFace(Axis ax, Direction dir, const Containers::vector<point_t> &fp);

Mesh::IndexedCell getXface(Direction dir, const Containers::vector<point_t> &fp);

Mesh::IndexedCell getYface(Direction dir, const Containers::vector<point_t> &fp);

Mesh::IndexedCell getZface(Direction dir, const Containers::vector<point_t> &fp);

} // namespace Utils

/**
 * Этот класс должен быть приватным классом кубической сетки, конечно же, потому что его можно спокойно создать в
 * неконсистентном состоянии. Но пока что это просто класс какой-то тут. Нужен для быстро посмотрения расчетов.
 */
struct IndexedCube {
    enum Axis { X = 0, Y = 1, Z = 2 };
    enum Direction { Minus = 0, Plus = 1 };

    using nodes_t = Containers::array<Types::index, 8>;
    using vertexes_t = Containers::vector<Types::point_t>;
    using full_points_t = Containers::vector<point_t>;

    point_t center_;
    Types::scalar volume_{};
    nodes_t nodes_;
    vertexes_t vertexes_;

    IndexedCube(const full_points_t &full_points, const nodes_t &full_indices);

    Mesh::IndexedCell getFace(Axis ax, Direction dir, const full_points_t &fp) const;

    Mesh::IndexedCell newGetFace(Axis ax, Direction dir, const full_points_t &fp) const;

    Mesh::IndexedCell getXface(Direction dir, const full_points_t &fp) const;

    Mesh::IndexedCell getYface(Direction dir, const full_points_t &fp) const;

    Mesh::IndexedCell getZface(Direction dir, const full_points_t &fp) const;
};

// Убрали много: посмотрим, как это повлияет на перфоманс
struct VertexParallelepiped {
    enum Axis { X = 0, Y = 1, Z = 2 };
    enum Direction { Minus = 0, Plus = 1 };

    using vertexes_t = Containers::vector<Types::point_t>;
    using full_points_t = Containers::vector<point_t>;
    using nodes_t = Containers::array<Types::index, 8>;

    // Количество хранимой памяти: 3 double* = 48 байт = 0.75 кэш-линий
    vertexes_t vertexes_;

    VertexParallelepiped(const full_points_t &full_points, const nodes_t &full_indices);

    [[nodiscard]] Mesh::IndexedCell getFace(Axis ax, Direction dir, const full_points_t &fp) const;
};
// сильно на перфоманс не повлияло

// Минимально, что нужно хранить -- это нижний угол куба, остальное одинаковое для всех кубов
struct MinimalParallelepiped {
    point_t min_corner;
};
// Такие "кубы" влезут по 4 шутки на 3 кэш-линии


} // namespace EMW::Mesh::VolumeCells

#endif // VOLUMECELLS_HPP
