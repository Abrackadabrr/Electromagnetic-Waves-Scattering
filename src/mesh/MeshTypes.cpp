//
// Created by evgen on 30.01.24.
//

#include "MeshTypes.hpp"

namespace EMW::Mesh {
    Cell::Cell(Containers::array<Point, 4> points) : points_(
            std::move(points)), collPoint_() {
        collPoint_.point_ = (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(4)) *
                            (points_[0] + points_[1] + points_[2] + points_[3]);
        const Types::Vector3d ac = (points_[2] - points_[0]);
        const Types::Vector3d bd = (points_[3] - points_[1]);
        const Types::Vector3d normalVector = ac.cross(bd);
        const Types::scalar n = normalVector.norm();
        area_ = (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) * n;

        // Задаем локальный базис на ПГП
        normal = normalVector / n;
        tau1 = ac / ac.norm();
        tau2 = normal.cross(tau1);
    };

    IndexedCell::IndexedCell(Containers::array<Types::index, 4> points, const Containers::vector<Point> &fullPoints)
            : points_(points), collPoint_() {
        collPoint_.point_ = (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(4)) *
                            (fullPoints[points_[0]] + fullPoints[points_[1]] + fullPoints[points_[2]] +
                             fullPoints[points_[3]]);

        cellStructure.A = fullPoints[points_[0]];
        cellStructure.ort1 = fullPoints[points_[1]] - fullPoints[points_[0]];
        cellStructure.ort1 = fullPoints[points_[3]] - fullPoints[points_[0]];
        cellStructure.diff =
                fullPoints[points_[0]] + fullPoints[points_[2]] - fullPoints[points_[1]] - fullPoints[points_[3]];

        integrationParameters.a = cellStructure.ort1.cross(cellStructure.ort2);
        integrationParameters.b = cellStructure.ort1.cross(cellStructure.diff);
        integrationParameters.c = cellStructure.diff.cross(cellStructure.ort2);



        const Types::Vector3d ac = (fullPoints[points_[2]] - fullPoints[points_[0]]);
        const Types::Vector3d bd = (fullPoints[points_[3]] - fullPoints[points_[1]]);
        const Types::Vector3d normalVector = ac.cross(bd);
        const Types::scalar n = normalVector.norm();
        area_ = (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) * n;

        // Задаем локальный базис на ПГП
        normal = normalVector / n;
        tau1 = ac / ac.norm();
        tau2 = normal.cross(tau1);
    };
}
