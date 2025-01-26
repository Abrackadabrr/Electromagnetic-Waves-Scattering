//
// Created by evgen on 19.08.24.
//
#if 1
#include "Algorithms.hpp"

namespace EMW::Mesh::Algorithm {
Types::scalar sign(Mesh::point_t p1, Mesh::point_t p2, Mesh::point_t p3) {
    return (p1.x() - p3.x()) * (p2.y() - p3.y()) - (p2.x() - p3.x()) * (p1.y() - p3.y());
}

bool PointInTriangle(const Mesh::point_t &pt, const Cell &cell) {
    Types::scalar d1, d2, d3, d4;
    bool has_neg, has_pos;

    d1 = sign(pt, cell.a, cell.b);
    d2 = sign(pt, cell.b, cell.c);
    d3 = sign(pt, cell.c, cell.d);
    d4 = sign(pt, cell.d, cell.a);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0) || (d4 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0) || (d4 > 0);

    return !(has_neg && has_pos);
}
}
#endif
