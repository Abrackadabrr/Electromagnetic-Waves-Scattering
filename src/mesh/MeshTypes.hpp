//
// Created by evgen on 29.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP

#include <utility>

#include "types/Types.hpp"

namespace EMW::Mesh {
using point_t = Types::Vector3d;

/**
 * Тип расчетного узла
 */
struct [[deprecated]] Node {
    using field_t = Types::Vector3c;
    point_t point_;
    field_t E_;
    field_t H_;
    field_t J_;

  public:
    Node() = default;

    explicit Node(point_t point)
        : point_(std::move(point)),
          E_(Types::Vector3c{Types::complex_d{0, 0}, Types::complex_d{0, 0}, Types::complex_d{0, 0}}),
          H_(Types::Vector3c{Types::complex_d{0, 0}, Types::complex_d{0, 0}, Types::complex_d{0, 0}}),
          J_(Types::Vector3c{Types::complex_d{0, 0}, Types::complex_d{0, 0}, Types::complex_d{0, 0}}){};

    Node(Types::scalar x, Types::scalar y, Types::scalar z, field_t E, field_t H, field_t J)
        : point_(x, y, z), E_(std::move(E)), H_(std::move(H)), J_(std::move(J)) {}

    void SetE(const field_t &E) { E_ = E; }

    void SetH(const field_t &H) { H_ = H; }

    void SetJ(const field_t &J) { H_ = J; }
};

struct CellStructure {
    point_t A; // тут перерасход памяти на хранение лишней точки, будет хорошо это исправить
    Types::Vector3d ort1; // B - A
    Types::Vector3d ort2; // D - A
    Types::Vector3d diff; // A + C - B - D
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
    using nodes_t = Containers::array<Types::index, 4>;
    nodes_t points_{};
    Types::scalar area_{};
    point_t collPoint_;
    Types::Vector3d normal;
    Containers::array<Types::Vector3d, 2> tau;
    CellStructure cellStructure;
    IntegrationParameters integrationParameters;

    IndexedCell() = default;

    IndexedCell(const Containers::array<Types::index, 4> &points, const Containers::vector<point_t> &fullPoints,
                const std::function<point_t(const Containers::array<Types::index, 4> &,
                                            const Containers::vector<point_t> &)> &getPoint);

    IndexedCell(const Containers::array<Types::index, 4> &points, const Containers::vector<point_t> &fullPoints)
        : IndexedCell(
              points, fullPoints,
              [&](const Containers::array<Types::index, 4> &p, const Containers::vector<point_t> &fp) -> point_t {
                  return (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(4)) *
                         (fp[p[0]] + fp[p[1]] + fp[p[2]] + fp[p[3]]);
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
