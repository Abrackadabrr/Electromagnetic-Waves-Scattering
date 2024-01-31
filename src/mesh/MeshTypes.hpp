//
// Created by evgen on 29.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP

#include "types/Types.hpp"

namespace EMW::Mesh {
    using Point = Types::Vector3d;

    /**
     * Тип точки коллокации
     */
    struct Node {
        using F_t = Types::Vector3d;
        Point point_;
        F_t E_;
        F_t H_;
        F_t J_;

    public:
        Node() : point_(0, 0, 0), E_(0.0, 0.0, 0.0), H_(0.0, 0.0, 0.0), J_(0.0, 0.0, 0.0) {}

        Node(Types::scalar x, Types::scalar y, Types::scalar z, F_t E, F_t H, F_t J)
                : point_(x, y, z), E_(std::move(E)), H_(std::move(H)), J_(std::move(J)) {}

        void SetE(const F_t &E) { E_ = E; }

        void SetH(const F_t &H) { H_ = H; }

        void SetJ(const F_t &J) { H_ = J; }
    };

    /**
     * Тип аппроксимированной ячейки сетки по четырём индексам
     * Содержит в себе вершины четырехугольника, точку коллокации, площадь
     * Хранение точек происходит в соотвествии с локальным полем нормалей поверхности ("правило буравчика")
     */
    struct IndexedCell {
        using nodes_t = Containers::array<Types::index, 4>;
        nodes_t points_;
        Types::scalar area_;
        Node collPoint_;
        Types::Vector3d normal;
        Types::Vector3d tau1;
        Types::Vector3d tau2;

        IndexedCell() = default;

        IndexedCell(Containers::array<Types::index, 4> points,
                    Types::scalar area,
                    Node collPoint) : points_(points), area_(area), collPoint_(std::move(collPoint)) {};

        IndexedCell(Containers::array<Types::index, 4> points, const Containers::vector<Point> &fullPoints);
    };

    /**
     * Тип аппроксимированной ячейки сетки по четырём точкам
     * Содержит в себе вершины четырехугольника, точку коллокации, площадь
     * Хранение точек происходит в соотвествии с локальным полем нормалей поверхности ("правило буравчика")
     */
    struct Cell {
        Containers::array<Point, 4> points_;
        Types::scalar area_;
        Node collPoint_;
        Types::Vector3d normal;
        Types::Vector3d tau1;
        Types::Vector3d tau2;

        Cell() = default;

        explicit Cell(Containers::array<Point, 4> points);
    };

}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP
