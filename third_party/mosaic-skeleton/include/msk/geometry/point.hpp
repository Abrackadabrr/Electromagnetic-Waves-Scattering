/*
 * Copyright 2023-2026 Bulat Valiakhmetov
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


/**
 * @file
 *
 * @brief Geometric primitives.
 */


#ifndef MOSAIC_SKELETON_GEOMETRY_POINT_HPP
#define MOSAIC_SKELETON_GEOMETRY_POINT_HPP


#include <cstdlib>
#include <cmath>
#include <complex>
#include <cstdint>
#include <ostream>
#include <array>


namespace MSk {

namespace geom {


/**
 * @brief Point (vector) in real Eucledean space.
 *
 * @tparam RealType Type of Point's coordinates.
 * @tparam dim Point (space) dimension.
 *
 * @note Complex type Point is fully supported as vector in a vector space,
 * except geometrical operations such as hyperplane side determinition.
 */
template<typename RealType, uint64_t dim>
class Point {
private:
    std::array<RealType, dim> coords;

public:
    using value_type = RealType;

public:
    /**
     * @brief Construct point with only first non-zero coordinate (vector of the given length along first axis).
     *
     * @param[in] val Value of the first coordinate.
     */
    Point(RealType val = 0) {
        coords[0] = val;
        for (uint64_t i = 1; i < dim; ++i) {
            coords[i] = 0;
        }
    }

    /**
     * @brief Initialize Point's first coordinates, others are filled with zeroes.
     *
     * @param[in] values Values of several first coordinates.
     */
    Point(std::initializer_list<RealType> values) {
        auto it = values.begin();
        uint64_t i{0};
        for (; i < dim && it != values.end(); ++i, ++it) {
            coords[i] = *it;
        }
        for (; i < dim; ++i) {
            coords[i] = 0;
        }
    }

    /**
     * @brief Fill Point's coordinates from a sequence.
     *
     * @tparam Iterator Sequence iterator type.
     *
     * @param[in] it_begin Beginning iterator.
     * @param[in] it_end Ending iterator.
     */
    template<class Iterator>
    explicit Point(Iterator it_begin, Iterator it_end) {
        uint64_t i = 0;
        for (Iterator it = it_begin; it != it_end && i < dim; ++it, ++i) {
            coords[i] = *(it);
        }
        for (; i < dim; ++i) {
            coords[i] = 0;
        }
    }

    /**
     * @brief Holder norm of the vector.
     *
     * @param[in] p Holder parameter, if zero then 'infinity' norm is returned.
     * @return Norm.
     */
    RealType norm(uint64_t p = 2) const {
        RealType nrm{0};
        if (p == 0) {
            for (uint64_t i = 0; i < dim; ++i) {
                if (std::abs(coords[i]) > nrm) {
                    nrm = std::abs(coords[i]);
                }
            }
        } else if (p == 1) {
            for (uint64_t i = 0; i < dim; ++i) {
                nrm += std::abs(coords[i]);
            }
        } else if (p == 2) {
            for (uint64_t i = 0; i < dim; ++i) {
                nrm += coords[i] * coords[i];
            }
            nrm = std::sqrt(nrm);
        } else {
            for (uint64_t i = 0; i < dim; ++i) {
                nrm += std::pow(std::abs(coords[i]), p);
            }
            nrm = std::pow(nrm, RealType{1} / p);
        }
        return nrm;
    }

    Point(const Point &other) = default;

    Point(Point &&other) noexcept = default;

    ~Point() = default;

    Point &operator=(const Point &other) = default;

    Point &operator=(Point &&other) noexcept = default;

    Point &operator+=(const Point &other) {
        for (uint64_t i = 0; i < dim; ++i) {
            coords[i] += other.coords[i];
        }
        return *this;
    }

    Point &operator-=(const Point &other) {
        for (uint64_t i = 0; i < dim; ++i) {
            coords[i] -= other.coords[i];
        }
        return *this;
    }

    Point &operator*=(RealType mult) {
        for (uint64_t i = 0; i < dim; ++i) {
            coords[i] *= mult;
        }
        return *this;
    }

    friend Point operator+(Point lhs, const Point &rhs) {
        lhs += rhs;
        return lhs;
    }

    friend Point operator-(Point lhs, const Point &rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend Point operator*(RealType p, const Point &rhs) {
        Point tmp(rhs);
        return tmp *= p;
    }

    friend Point operator*(const Point &lhs, RealType p) {
        Point tmp(lhs);
        return tmp *= p;
    }

    RealType &operator[](std::size_t idx) {
        return coords[idx];
    }

    const RealType &operator[](std::size_t idx) const {
        return coords[idx];
    }

    /**
     * @brief Euclidean dot (inner) product of two vectors.
     *
     * @param[in] p1 First vector.
     * @param[in] p2 Second vector.
     * @return Dot protuct
     */
    friend RealType dot_product(const Point &p1, const Point &p2) {
        RealType dot{0};
        for (uint64_t i = 0; i < dim; ++i) {
            dot += p1.coords[i] * p2.coords[i];
        }
        return dot;
    }

    /**
     * @brief Determine on which side of the given hyperplane the Point lies.
     *
     * @param[in] dir Directional vector of the hyperplane.
     * @param[in] start One point of the hyperplane.
     * @return Boolean side.
     */
    bool hplane_side(const Point &normal, const Point &start) const {
        Point p = *this - start;
        RealType proj = std::real(dot_product(p, normal));
        return proj >= 0;
    }

    friend std::ostream &operator<<(std::ostream &os, const Point &obj) {
        os << "(";
        for (uint64_t i = 0; i < dim - 1; ++i) {
            os << obj.coords[i] << ", ";
        }
        os << obj.coords[dim - 1] << ")";
        return os;
    }

    /**
     * @brief Dimension of th Point.
     *
     * @return Dimension.
     */
    uint64_t size() const {
        return dim;
    }

};


/**
 * @brief Generate uniform grid of points on a segment.
 *
 * @tparam RealType Points coordinates type.
 * @tparam dim Points dimension.
 *
 * @param[out] points Resulting points array.
 * @param[in] n_pts Number of points on the segment (boundaries included).
 * @param[in] start First boundary point.
 * @param[in] finish Second boundary point.
 */
template<typename RealType, uint64_t dim>
void generate_line(Point<RealType, dim> *points, uint64_t n_pts,
                   const Point<RealType, dim> &start = {0}, const Point<RealType, dim> &finish = {1}) {
    uint64_t n_cuts{n_pts - 1};
    for (uint64_t i = 0; i <= n_cuts; ++i) {
        points[i] = (RealType(n_cuts - i) / n_cuts) * start + (RealType(i) / n_cuts) * finish;
    }
}

/**
 * @brief Generate uniform grid on a two dimensional circle.
 *
 * @tparam RealType Points coordinates type.
 *
 * @param[out] points Resulting points array.
 * @param[in] n_pts Number of points on the circle.
 * @param[in] radius Circle radius.
 * @param[in] center Circle center.
 */
template<typename RealType>
void generate_circle(Point<RealType, 2> *points, uint64_t n_pts,
                     const RealType radius, const Point<RealType, 2> &center = {0}) {
    for (uint64_t i = 0; i < n_pts; ++i) {
        double phi = (2.0 * M_PI * i) / n_pts;
        points[i][0] = center[0] + radius * std::cos(phi);
        points[i][1] = center[1] + radius * std::sin(phi);
    }
}

/**
 * @brief Generate uniform grid of points in a multidimensional cube.
 *
 * @tparam RealType Points coordinates type.
 * @tparam dim Points dimension.
 *
 * @param[out] points Resulting points array.
 * @param[in] n_side Number of points along each cube's side (boundaries included).
 * @param[in] cube_dim Dimension of a cube (1 - line, 2 - square, 3 - cube, ...)
 * @param[in] corner1 First corner of the cube.
 * @param[in] corner2 Second corner of the cube.
 */
template<typename RealType, uint64_t dim>
void generate_cube(Point<RealType, dim> *points, uint64_t n_side, uint64_t cube_dim,
                   const Point<RealType, dim> &corner1 = {0}, const Point<RealType, dim> &corner2 = {1}) {
    if (cube_dim == 0) {
        return;
    }
    if (cube_dim == 1) {
        generate_line(points, n_side, corner1, corner2);
    } else {
        uint64_t n_pts{n_side};
        for (uint64_t i = 1; i < cube_dim - 1; ++i) {
            n_pts *= n_side;
        }
        RealType last_dim = corner2[cube_dim - 1] - corner1[cube_dim - 1];
        uint64_t n_cuts{n_side - 1};
        for (uint64_t i = 0; i <= n_cuts; ++i) {
            auto new_corner1 = corner1;
            auto new_corner2 = corner2;
            new_corner1[cube_dim - 1] += last_dim * ((RealType{1.0} * i) / n_cuts);
            new_corner2[cube_dim - 1] = new_corner1[cube_dim - 1];
            generate_cube(points + i * n_pts, n_side, cube_dim - 1, new_corner1, new_corner2);
        }
    }
}


extern template class Point<float, 2>;

extern template class Point<double, 2>;


extern template class Point<float, 3>;

extern template class Point<double, 3>;


extern template void generate_line<float, 2>(Point<float, 2> *points, uint64_t n_pts,
        const Point<float, 2> &start, const Point<float, 2> &finish);

extern template void generate_line<double, 2>(Point<double, 2> *points, uint64_t n_pts,
        const Point<double, 2> &start, const Point<double, 2> &finish);


extern template void generate_line<float, 3>(Point<float, 3> *points, uint64_t n_pts,
        const Point<float, 3> &start, const Point<float, 3> &finish);

extern template void generate_line<double, 3>(Point<double, 3> *points, uint64_t n_pts,
        const Point<double, 3> &start, const Point<double, 3> &finish);


extern template void generate_circle<float>(Point<float, 2> *points, uint64_t n_pts,
        const float radius, const Point<float, 2> &center);

extern template void generate_circle<double>(Point<double, 2> *points, uint64_t n_pts,
        const double radius, const Point<double, 2> &center);


extern template void generate_cube<float, 2>(Point<float, 2> *points, uint64_t n_side, uint64_t cube_dim,
                   const Point<float, 2> &corner1, const Point<float, 2> &corner2);

extern template void generate_cube<double, 2>(Point<double, 2> *points, uint64_t n_side, uint64_t cube_dim,
                   const Point<double, 2> &corner1, const Point<double, 2> &corner2);


extern template void generate_cube<float, 3>(Point<float, 3> *points, uint64_t n_side, uint64_t cube_dim,
                   const Point<float, 3> &corner1, const Point<float, 3> &corner2);

extern template void generate_cube<double, 3>(Point<double, 3> *points, uint64_t n_side, uint64_t cube_dim,
                   const Point<double, 3> &corner1, const Point<double, 3> &corner2);


}
}


#endif // MOSAIC_SKELETON_GEOMETRY_POINT_HPP

