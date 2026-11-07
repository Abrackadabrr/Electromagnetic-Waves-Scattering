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
// Created by bulat_v on 20.08.23.
//


/**
 * @file
 *
 * @brief Defines template class geom::PointGrid.
 */


#ifndef MOSAIC_SKELETON_GEOMETRY_POINT_GRID_HPP
#define MOSAIC_SKELETON_GEOMETRY_POINT_GRID_HPP


#include <string>
#include <fstream>
#include <limits>

#include "smart_array/smart_array.hpp"

#include "msk/geometry/point.hpp"


namespace MSk {

namespace geom {


/**
 * @brief Arbitrary grid of points.
 *
 * @tparam DataType Points coordinates type.
 * @tparam dim Points dimension.
 */
template<typename DataType, uint64_t dim>
class PointGrid {
private:
    uint64_t _n_pts{0};
    SmartArray<Point<DataType, dim>> _own_pts; ///< Used if points are not given as an array and need to be stored
    const Point<DataType, dim> *_pts{nullptr};

public:
    using obj_type = Point<DataType, dim>;

private:
    /**
     * @brief Skip starting comment lines starting with '#' in a text file.
     *
     * @param[in, out] fin Open file object.
     */
    static void skip_comments(std::ifstream &fin) {
        for (auto c = fin.peek(); c == '#'; c = fin.peek()) {
            fin.ignore(std::numeric_limits<std::streamsize>::max(), fin.widen('\n'));
        }
    }

    /**
     * @brief Initialize the grid from an input text file.
     *
     * @details First reads the number of points, then sequentially reads points by coordinate.
     *
     * @param[in] fin Open file object.
     */
    void init_txt(std::ifstream &fin) {
        skip_comments(fin);
        uint32_t n_pts{};
        fin >> n_pts;
        _own_pts.resize(n_pts);
        DataType *ptr = reinterpret_cast<DataType *>(_own_pts.data());
        for (uint64_t k = 0; k < n_pts * dim; ++k) {
            fin >> ptr[k];
        }
    }

    /**
     * @brief Initialize the grid from an inputbinary file.
     *
     * @details First reads the number (uint32_t) of points, then sequentially reads points by coordinate.
     *
     * @param[in] fin Open file object.
     */
    void init_bin(std::ifstream &fin) {
        uint32_t n_pts{};
        fin.read(reinterpret_cast<char *>(&n_pts), sizeof(uint32_t));
        _own_pts.resize(n_pts);
        fin.read(reinterpret_cast<char *>(_own_pts.data()), sizeof(obj_type) * n_pts);
    }


public:
    /**
     * @brief Empty grid.
     */
    PointGrid() = default;

    /**
     * @brief Initialize with the uniform grid on the segment [0, 1].
     *
     * @param[in] num_points Number of points on the segment (boundaries included).
     */
    explicit PointGrid(uint64_t num_points) : _n_pts{num_points}, _own_pts{num_points}, _pts{_own_pts.data()} {
        if (num_points == 0) {
            _pts = nullptr;
        } else {
            generate_line(_own_pts.data(), _n_pts);
        }
    }

    /**
     * @brief Initialize the grid with array of points.
     *
     * @details Points are not copied, only pointer is stored.
     *
     * @param[in] num_points Number of points.
     * @param[in] points_ptr Points array
     */
    PointGrid(uint64_t num_points, const Point<DataType, dim> *points_ptr) : _n_pts{num_points}, _pts{points_ptr} {}

    /**
     * @brief Initialize the grid with array of points.
     *
     * @details Points are not copied directly as SmartArray does such work itself.
     *
     * @param[in] points SmartArray of points.
     */
    PointGrid(const SmartArray<Point<DataType, dim>> &points) : _n_pts{points.size()}, _own_pts{points},
                                                                _pts{_own_pts.data()} {}

    /**
     * @brief Initialize the grid with points from an input file.
     *
     * @details If file's name is ended with '.bin' then it is treated as binary, textual otherwise.
     * File should be formated in the follownig way. Textual format allows beginning comment lines
     * starting with '#' symbol. The first value is the number of points. Then there should come points
     * by their coordinates.
     *
     * @param[in] filename Name of a file.
     */
    PointGrid(const std::string &filename) {
        if (filename.size() > 4 && filename.compare(filename.size() - 4, 4, ".bin") == 0) {
            std::ifstream fin(filename, std::ios_base::in | std::ios_base::binary);
            init_bin(fin);
        } else {
            std::ifstream fin(filename, std::ios_base::in);
            init_txt(fin);
        }
        _n_pts = _own_pts.size();
        _pts = _own_pts.data();
    }

    PointGrid(const PointGrid &) = default;

    PointGrid(PointGrid &&) = default;

    PointGrid &operator =(const PointGrid &) = default;

    PointGrid &operator =(PointGrid&&) = default;

    const obj_type &operator[](uint64_t i) const {
        if (i >= _n_pts) {
            throw std::out_of_range("Geometry's Point index out of range");
        }
        return _pts[i];
    }

    /**
     * @brief Get the number of points in the grid.
     *
     * @return Number of points.
     */
    uint64_t size() const {
        return _n_pts;
    }
};


extern template class PointGrid<float, 2>;

extern template class PointGrid<double, 2>;


extern template class PointGrid<float, 3>;

extern template class PointGrid<double, 3>;


}
}
#endif //MOSAIC_SKELETON_GEOMETRY_POINT_GRID_HPP
