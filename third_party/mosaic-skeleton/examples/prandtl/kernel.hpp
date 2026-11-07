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
// Created by b.valiakhmetov on 5/21/24.
//


/**
 * @file
 *
 * @brief Prandtl kernel.
 */


#ifndef MOSAIC_SKELETON_EXAMPLES_PRANDTL_KERNEL_HPP
#define MOSAIC_SKELETON_EXAMPLES_PRANDTL_KERNEL_HPP

#define _USE_MATH_DEFINES


#include <cstdint>
#include <cmath>

#include "msk/matrix/matrix.hpp"
#include "msk/geometry/point.hpp"


/**
 * @brief Magnetic dipole field detection model matrix.
 *
 * @tparam RealType Element type.
 */
template<class RealType>
class Prandtl: public MSk::matrix::Matrix<RealType> {
private:
    bool _chebyshev{true};
    uint64_t _size_x, _size_y;
    RealType _length{1.0};

public:
    using obj_type = MSk::geom::Point<double, 2>; ///< Spatial point type

private:
    RealType node(RealType level) const {
        return _length * level;;
    }

    RealType c_node(RealType level) const {
        return _length * (1.0 - std::cos(M_PI * level)) / 2.0;
    }

    RealType part(const MSk::geom::Point<RealType, 2> &vec) const {
        return vec.norm() / (vec[0] * vec[1]);
    }

public:
    Prandtl(uint64_t size_x, uint64_t size_y, const std::string &type = "uniform") :
            MSk::matrix::Matrix<RealType>(size_x * size_y, size_x * size_y),
            _size_x{size_x}, _size_y{size_y},
            _chebyshev{type == "chebyshev" || type == "cheb" ? true : false} {}

    ~Prandtl() override = default;

    RealType element(uint64_t i, uint64_t j) const override {
        uint64_t src_idx{j / _size_y}, src_idy{j % _size_y};
        obj_type colloc = row_point(i);
        RealType xA, xB, yA, yB;
        if (_chebyshev) {
            xA = c_node(1.0 * src_idx / _size_x);
            xB = c_node(1.0 * (src_idx + 1) / _size_x);
            yA = c_node(1.0 * src_idy / _size_y);
            yB = c_node(1.0 * (src_idy + 1) / _size_y);
        } else {
            xA = node(1.0 * src_idx / _size_x);
            xB = node(1.0 * (src_idx + 1) / _size_x);
            yA = node(1.0 * src_idy / _size_y);
            yB = node(1.0 * (src_idy + 1) / _size_y);
        }
        return part(obj_type{xA, yA} - colloc) + part(obj_type{xB, yB} - colloc) - part(obj_type{xA, yB} - colloc) - part(obj_type{xB, yA} - colloc);
    }

    /**
     * @brief Receiver point corresponding to a row index.
     *
     * @param[in] i Row index.
     * @return geom::Point.
     */
    obj_type row_point(uint64_t i) const {
        uint64_t rcv_idx{i / _size_y}, rcv_idy{i % _size_y};
        RealType xA, xB, yA, yB;
        if (_chebyshev) {
            xA = c_node(1.0 * rcv_idx / _size_x);
            xB = c_node(1.0 * (rcv_idx + 1) / _size_x);
            yA = c_node(1.0 * rcv_idy / _size_y);
            yB = c_node(1.0 * (rcv_idy + 1) / _size_y);
        } else {
            xA = node(1.0 * rcv_idx / _size_x);
            xB = node(1.0 * (rcv_idx + 1) / _size_x);
            yA = node(1.0 * rcv_idy / _size_y);
            yB = node(1.0 * (rcv_idy + 1) / _size_y);
        }
        return {(xA + xB) * 0.5, (yA + yB) * 0.5};
    }

    /**
     * @brief Source point corresponding to a column index.
     *
     * @param[in] j Column index.
     * @return geom::Point.
     */
    obj_type col_point(uint64_t j) const {
        uint64_t src_idx{j / _size_y}, src_idy{j % _size_y};
        RealType xA, xB, yA, yB;
        if (_chebyshev) {
            xA = c_node(1.0 * src_idx / _size_x);
            xB = c_node(1.0 * (src_idx + 1) / _size_x);
            yA = c_node(1.0 * src_idy / _size_y);
            yB = c_node(1.0 * (src_idy + 1) / _size_y);
        } else {
            xA = node(1.0 * src_idx / _size_x);
            xB = node(1.0 * (src_idx + 1) / _size_x);
            yA = node(1.0 * src_idy / _size_y);
            yB = node(1.0 * (src_idy + 1) / _size_y);
        }
        return {(xA + xB) * 0.5, (yA + yB) * 0.5};
    }
};


#endif //MOSAIC_SKELETON_EXAMPLES_PRANDTL_KERNEL_HPP
