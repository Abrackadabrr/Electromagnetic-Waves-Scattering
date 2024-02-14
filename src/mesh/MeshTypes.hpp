//
// Created by evgen on 29.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP

#include "types/Types.hpp"

namespace EMW::Mesh {
    using Point = Types::Vector3d;

    /**
     * Тип расчетного узла
     */
    struct Node {
        using F_t = Types::Vector3c;
        Point point_;
        F_t E_;
        F_t H_;
        F_t J_;

    public:
        Node() = default;

        explicit Node(const Point &point) : point_(point),
                                            E_(Types::Vector3c{Types::complex_d{0, 0}, Types::complex_d{0, 0},
                                                               Types::complex_d{0, 0}}),
                                            H_(Types::Vector3c{Types::complex_d{0, 0}, Types::complex_d{0, 0},
                                                               Types::complex_d{0, 0}}),
                                            J_(Types::Vector3c{Types::complex_d{0, 0}, Types::complex_d{0, 0},
                                                               Types::complex_d{0, 0}}) {};

        Node(Types::scalar x, Types::scalar y, Types::scalar z, F_t E, F_t H, F_t J)
                : point_(x, y, z), E_(std::move(E)), H_(std::move(H)), J_(std::move(J)) {}

        void SetE(const F_t &E) { E_ = E; }

        void SetH(const F_t &H) { H_ = H; }

        void SetJ(const F_t &J) { H_ = J; }
    };

    struct CellStructure {
        Point A;  // тут перерасход памяти на хранение лишней точки, будет хорошо это исправить
        Types::Vector3d ort1;  // B - A
        Types::Vector3d ort2;  // D - A
        Types::Vector3d diff;  // A + C - B - D
    };

    struct IntegrationParameters {
        // Parameters for surface integral
        Types::Vector3d a;  // ort1 x ort2
        Types::Vector3d b;  // ort1 x diff
        Types::Vector3d c;  // diff x ort2

        // Parameters for curve integral
        Containers::array<Types::Vector3d, 4> mul;
    };

    /**
     * Тип аппроксимированной ячейки сетки по четырём индексам
     * Содержит в себе вершины четырехугольника, точку коллокации, площадь
     * Хранение точек происходит в соотвествии с локальным полем нормалей поверхности ("правило буравчика")
     *
     * В функцию для интегрирования заходит по константной ссылке (оспариваемое решение)
     */
    struct IndexedCell {
        using nodes_t = Containers::array<Types::index, 4>;
        nodes_t points_;
        Types::scalar area_;
        Node collPoint_;
        Types::Vector3d normal;
        Containers::array<Types::Vector3d, 2> tau;
        CellStructure cellStructure;
        IntegrationParameters integrationParameters;

        IndexedCell() = default;

        IndexedCell(Containers::array<Types::index, 4> points, const Containers::vector<Point> &fullPoints);

        [[nodiscard]] Point parametrization(Types::scalar p, Types::scalar q) const noexcept {
            return cellStructure.A + p * cellStructure.ort1 + q * cellStructure.ort2 + p * q * cellStructure.diff;
        }

        [[nodiscard]] Types::scalar multiplier(Types::scalar p, Types::scalar q) const noexcept {
            return (integrationParameters.a + p * integrationParameters.b + q * integrationParameters.c).norm();
        }
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
        Containers::array<Types::Vector3d, 2> tau;
        CellStructure cellStructure;
        IntegrationParameters integrationParameters;

        Cell() = default;

        explicit Cell(Containers::array<Point, 4> points);

        [[nodiscard]] Point parametrization(Types::scalar p, Types::scalar q) const noexcept {
            return cellStructure.A + p * cellStructure.ort1 + q * cellStructure.ort2 + p * q * cellStructure.diff;
        }

        [[nodiscard]] Types::scalar multiplier(Types::scalar p, Types::scalar q) const noexcept {
            return (integrationParameters.a + p * integrationParameters.b + q * integrationParameters.c).norm();
        }
    };

}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_MESHTYPES_HPP
