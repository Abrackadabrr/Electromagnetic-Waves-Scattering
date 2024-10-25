//
// Created by evgen on 30.01.24.
//

#include "MeshTypes.hpp"

namespace EMW::Mesh {
    IndexedCell::IndexedCell(Containers::array<Types::index, 4> points, const Containers::vector<point_t> &fullPoints)
            : points_(points) {
        // собрали точку коллокации из имеющихся базисных векторов
        collPoint_.point_ = (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(4)) *
                            (fullPoints[points_[0]] + fullPoints[points_[1]] + fullPoints[points_[2]] +
                             fullPoints[points_[3]]);

        cellStructure.A = fullPoints[points_[0]];
        cellStructure.ort1 = fullPoints[points_[1]] - fullPoints[points_[0]];
        cellStructure.ort2 = fullPoints[points_[3]] - fullPoints[points_[0]];
        cellStructure.diff =
                fullPoints[points_[0]] + fullPoints[points_[2]] - fullPoints[points_[1]] - fullPoints[points_[3]];

        integrationParameters.a = cellStructure.ort1.cross(cellStructure.ort2);  // == normal * area
        integrationParameters.b = cellStructure.ort1.cross(cellStructure.diff);
        integrationParameters.c = cellStructure.diff.cross(cellStructure.ort2);

        const Types::Vector3d ac = (fullPoints[points_[2]] - fullPoints[points_[0]]);
        const Types::Vector3d bd = (fullPoints[points_[3]] - fullPoints[points_[1]]);
        const Types::Vector3d normalVector = ac.cross(bd);
        const Types::scalar n = normalVector.norm();
        area_ = (static_cast<Types::scalar>(1) / static_cast<Types::scalar>(2)) * n;

        // Задаем локальный базис на ПГП
        normal = normalVector / n;
        tau[0] = ac.normalized();
        tau[1] = normal.cross(tau[0]);

        // parameters for integral over A->B
        integrationParameters.mul[0] = cellStructure.ort1.norm() * (cellStructure.ort1.cross(
                normal)).normalized();   // == omega_i * -(ort_2).normalized()
        // parameters for integral over B->C
        integrationParameters.mul[1] = (fullPoints[points_[2]] - fullPoints[points_[1]]).norm() *
                                       ((fullPoints[points_[2]] - fullPoints[points_[1]]).cross(normal)).normalized();
        // parameters for integral over C->D
        integrationParameters.mul[2] = (fullPoints[points_[3]] - fullPoints[points_[2]]).norm() *
                                       ((fullPoints[points_[3]] - fullPoints[points_[2]]).cross(normal)).normalized();
        // parameters for integral over D->A
        integrationParameters.mul[3] = cellStructure.ort2.norm() * (normal.cross(
                cellStructure.ort2)).normalized();  // == omega_i * (ort_1).normalized()

        // костыль для плоской геометрии в OYZ
//        if (points[0] % 2 != 0) {
//            tau[0] = {0, 1, 0};
//            tau[1] = {0, 0, 1};
//        } else {
//            tau[0] = Types::Vector3d{0, 0, 1}.normalized();
//            tau[1] = Types::Vector3d{0, -1, 0}.normalized();
//        }
    };


    Cell IndexedCell::getVertex() const{
        Cell result;
        result.a = cellStructure.A;
        result.b = result.a + cellStructure.ort1;
        result.d = result.a + cellStructure.ort2;
        result.c = cellStructure.diff - result.a + result.b + result.d;
        return result;
    }

    Containers::array<Mesh::point_t, 4> IndexedCell::getVertexAsArray() const{
        Cell result;
        result.a = cellStructure.A;
        result.b = result.a + cellStructure.ort1;
        result.d = result.a + cellStructure.ort2;
        result.c = cellStructure.diff - result.a + result.b + result.d;
        return {result.a, result.b, result.c, result.d};
    }
}
