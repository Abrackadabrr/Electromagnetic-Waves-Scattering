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
//
// Created by b.valiakhmetov on 2/28/24.
//


/**
 * @file
 *
 * @brief Geometric metrics.
 */


#ifndef MOSAIC_SKELETON_GEOMETRY_METRICS_HPP
#define MOSAIC_SKELETON_GEOMETRY_METRICS_HPP


#include <cstdint>
#include <algorithm>
#include <functional>
#include <limits>
#include <utility>


namespace MSk {

namespace geom {


template<typename Iter, typename Point>
Point centroid(std::function<Point(Iter)> &point_getter, Iter begin, Iter end) {
    Point center{0};
    uint64_t len{0};
    for (Iter it = begin; it != end; ++it) {
        center += point_getter(it);
        ++len;
    }
    using DataType = decltype(center.norm());
    center *= DataType{1.0} / len;
    return center;
}


template<typename Iter, typename Point>
std::pair<Point, Point> bound_box(std::function<Point(Iter)> &point_getter, Iter begin, Iter end) {
    const uint64_t dim = Point().size();
    Point corn_min{point_getter(begin)}, corn_max{point_getter(begin)};
    for (Iter it = begin; it != end; ++it) {
        Point point = point_getter(it);
        for (uint64_t i = 0; i < dim; ++i) {
            corn_min[i] = std::min(corn_min[i], point[i]);
            corn_max[i] = std::max(corn_max[i], point[i]);
        }
    }
    return {corn_min, corn_max};
}


template<typename Iter, typename Point>
Point bound_sphere(std::function<Point(Iter)> &point_getter, Iter begin, Iter end, const Point &center, uint64_t p_norm = 2) {
    // Calculate farthest from center
    Point radial{0};
    decltype(radial.norm(p_norm)) max_rad = 0;
    for (Iter it = begin; it != end; ++it) {
        Point new_radial = point_getter(it) - center;
        auto rad = new_radial.norm(p_norm);
        if (rad > max_rad) {
            max_rad = rad;
            radial = new_radial;
        }
    }
    return radial;
}


template<typename Iter, typename Point>
Point iter_elong(std::function<Point(Iter)> &point_getter, Iter begin, Iter end, const Point &start_point, uint64_t p_norm = 2,
                 uint64_t max_iter = 4, decltype(Point().norm()) tol = std::numeric_limits<decltype(Point().norm())>::epsilon()) {
    Point p_far{start_point}, p_from{start_point};
    auto max_dist = (p_far - p_from).norm(p_norm);
    bool iter_update{true};
    for (uint64_t i = 0; iter_update && i < max_iter; ++i) {
        iter_update = false;
        // Find farthest point from p_from
        for (Iter it = begin; it != end; ++it) {
            Point p_far_new = point_getter(it);
            auto dist = (p_far_new - p_from).norm(p_norm);
            if (dist > max_dist + tol) {
                max_dist = dist;
                p_far = p_far_new;
                iter_update = true;
            }
        }
        p_from = p_far; // Set farthest point as the new starting point
    }
    p_far -= p_from;
    if (p_far.norm(p_norm) < tol) {
        return Point{1};
    }
    using DataType = decltype(p_far.norm());
    p_far *= DataType{1.0} / p_far.norm(p_norm);
    return p_far;
}


}

}



#endif //MOSAIC_SKELETON_GEOMETRY_METRICS_HPP
