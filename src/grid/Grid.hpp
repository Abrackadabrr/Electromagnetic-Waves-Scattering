//
// Created by evgen on 25.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_GRID_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_GRID_HPP

#include <cmath>
#include <utility>
#include "Types.hpp"
#include "GridTypes.hpp"

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

namespace EMW::Grid {
    /**
     * Тип точки (узел сетки, точка коллокации)
     */
    struct Node {
        using F_t = Types::Vector3d;
        Point point_;
        F_t E_;
        F_t H_;

    public:
        Node() : point_(0, 0, 0), E_(0.0, 0.0, 0.0), H_(0.0, 0.0, 0.0) {}

        Node(double x, double y, double z, F_t E, F_t H)
                : point_(x, y, z), E_(std::move(E)), H_(std::move(H)) {}

        void SetE(const F_t &E) { E_ = E; }
        void SetH(const F_t &H) { H_ = H; }
    };

    /**
     * Тип ячейки сетки по индексу
     */
    struct IndexedCell {
        /**
         * Массив индексов точек сетки в
         */
        Containers::array<Types::index, 4> points;
    };

    /**
     * Тип ячейки сетки по двум
     */
    struct Cell {
        /**
         * Массив индексов точек сетки, которые являтся
         */
        Containers::array<Point, 4> points;
    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_GRID_HPP
