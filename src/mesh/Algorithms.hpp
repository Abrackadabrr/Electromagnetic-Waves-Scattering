//
// Created by evgen on 19.08.24.
//

#ifndef PRECONDITIONING_ALGORITHMS_HPP
#define PRECONDITIONING_ALGORITHMS_HPP

#include "MeshTypes.hpp"
#include "types/Types.hpp"

namespace EMW::Mesh::Algorithm {
#if 0
    template <typename I>
    constexpr I logical_nor(I lhs, I rhs) noexcept { return not (lhs or rhs); }

    /**
     * Ориентированная площадь треугольника на трёх точках p1 p2 p3
     * Фича: если больше 0, то p1 лежит справа от p1 - p2
     * @return (p1 - p3) x (p2 - p3)
     */
    Types::scalar sign(Mesh::point_t p1, Mesh::point_t p2, Mesh::point_t p3);

    /**
     * Проверка утверждения "точка находится в ячейке"
     * @param pt точка
     * @param cell ячейка (c == d)
     */
    bool PointInTriangle(const Mesh::point_t &pt, const Cell &cell);
#endif
}

#endif //PRECONDITIONING_ALGORITHMS_HPP
