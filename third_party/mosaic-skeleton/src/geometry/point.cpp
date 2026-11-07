/*
 * Copyright 2024-2026 Bulat Valiakhmetov
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
// Created by b.valiakhmetov on 4/10/24.
//


#include "msk/geometry/point.hpp"


template class MSk::geom::Point<float, 2>;

template class MSk::geom::Point<double, 2>;


template class MSk::geom::Point<float, 3>;

template class MSk::geom::Point<double, 3>;


template void MSk::geom::generate_line<float, 2>(Point<float, 2> *points, uint64_t n_pts,
        const Point<float, 2> &start, const Point<float, 2> &finish);

template void MSk::geom::generate_line<double, 2>(Point<double, 2> *points, uint64_t n_pts,
        const Point<double, 2> &start, const Point<double, 2> &finish);


template void MSk::geom::generate_line<float, 3>(Point<float, 3> *points, uint64_t n_pts,
        const Point<float, 3> &start, const Point<float, 3> &finish);

template void MSk::geom::generate_line<double, 3>(Point<double, 3> *points, uint64_t n_pts,
        const Point<double, 3> &start, const Point<double, 3> &finish);


template void MSk::geom::generate_circle<float>(Point<float, 2> *points, uint64_t n_pts,
        const float radius, const Point<float, 2> &center);

template void MSk::geom::generate_circle<double>(Point<double, 2> *points, uint64_t n_pts,
        const double radius, const Point<double, 2> &center);


template void MSk::geom::generate_cube<float, 2>(Point<float, 2> *points, uint64_t n_side, uint64_t cube_dim,
                   const Point<float, 2> &corner1, const Point<float, 2> &corner2);

template void MSk::geom::generate_cube<double, 2>(Point<double, 2> *points, uint64_t n_side, uint64_t cube_dim,
                   const Point<double, 2> &corner1, const Point<double, 2> &corner2);


template void MSk::geom::generate_cube<float, 3>(Point<float, 3> *points, uint64_t n_side, uint64_t cube_dim,
                   const Point<float, 3> &corner1, const Point<float, 3> &corner2);

template void MSk::geom::generate_cube<double, 3>(Point<double, 3> *points, uint64_t n_side, uint64_t cube_dim,
                   const Point<double, 3> &corner1, const Point<double, 3> &corner2);
