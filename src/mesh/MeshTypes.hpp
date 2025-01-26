//
// Created by evgen on 29.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP

#include <utility>

#include "types/Types.hpp"

namespace EMW::Mesh {
using point_t = EMW::Types::Vector3d;

struct CellStructure {
    point_t A; // тут перерасход памяти на хранение лишней точки, будет хорошо это исправить
    Types::Vector3d ort1; // B - A
    Types::Vector3d ort2; // D - A
    Types::Vector3d diff; // A + C - B - D

    bool operator==(const CellStructure&) const = default;
};

struct IntegrationParameters {
    // Parameters for surface integral
    Types::Vector3d a; // ort1 x ort2
    Types::Vector3d b; // ort1 x diff
    Types::Vector3d c; // diff x ort2

    // Parameters for curve integral
    Containers::array<Types::Vector3d, 4> mul;
};

struct Cell {
    point_t a, b, c, d;
};

/**
 * Тип аппроксимированной ячейки сетки по четырём индексам
 * Содержит в себе вершины четырехугольника, точку коллокации, площадь
 * Хранение точек происходит в соотвествии с локальным полем нормалей поверхности ("правило буравчика")
 */
struct IndexedCell {

    enum class Tag {
        NO_TAG = 0,
        SIGMA = 1,
        WAVEGUIDE_CROSS_SECTION = 2,
        MAGNETIC_CURRENTS = 3,
    };

    using nodes_t = Containers::array<Types::index, 4>;
    nodes_t points_{};
    Types::scalar area_{};
    point_t collPoint_;
    Types::Vector3d normal;
    Containers::array<Types::Vector3d, 2> tau;
    CellStructure cellStructure;
    IntegrationParameters integrationParameters;

    Tag tag = Tag::SIGMA;

    IndexedCell() = default;

    IndexedCell(const nodes_t &points, const Containers::vector<point_t> &fullPoints,
                const std::function<point_t(const nodes_t &,
                                            const Containers::vector<point_t> &)> &getPoint);

    IndexedCell(const Containers::array<Types::index, 4> &points, const Containers::vector<point_t> &fullPoints)
        : IndexedCell(
              points, fullPoints,
              [&](const Containers::array<Types::index, 4> &p, const Containers::vector<point_t> &fp) -> point_t {
                  return (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(4)) *
                         (fp.at(p[0]) + fp.at(p[1]) + fp.at(p[2]) + fp.at(p[3]));
              }) {}

    [[nodiscard]] point_t parametrization(Types::scalar p, Types::scalar q) const noexcept {
        return cellStructure.A + p * cellStructure.ort1 + q * cellStructure.ort2 + p * q * cellStructure.diff;
    }

    [[nodiscard]] Types::scalar multiplier(Types::scalar p, Types::scalar q) const noexcept {
        return (integrationParameters.a + p * integrationParameters.b + q * integrationParameters.c).norm();
    }

    [[nodiscard]] Cell getVertex() const;

    [[nodiscard]]Containers::array<Mesh::point_t, 4> getVertexAsArray() const;
};
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP
