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
// Created by b.valiakhmetov on 17.05.23.
//


/**
 * @file
 *
 * @brief Contains magnetic matrix example.
 */


#ifndef MOSAIC_SKELETON_MATRIX_MAGNETIC_HPP
#define MOSAIC_SKELETON_MATRIX_MAGNETIC_HPP


#include <fstream>
#include <string>

#include "smart_array/smart_array.hpp"

#include "msk/matrix/matrix.hpp"
#include "msk/geometry/point.hpp"


namespace MSk {


/**
 * @brief Magnetic dipole field detection model matrix.
 *
 * @tparam DataType Element type.
 */
template<class DataType>
class Magnetic: public matrix::Matrix<DataType> {
private:
    uint64_t _n_rcv, _n_src; ///< Number of receivers and sources (points)
    bool _field, _grad; ///< If to evaluate field and gradient
    uint64_t _b_rows{1}, _b_cols{1}; ///< Sizes of the block of two points interaction
    SmartArray<geom::Point<double, 3>> _pts_rcv, _pts_src; ///< Arrays of points in space

public:
    using obj_type = geom::Point<double, 3>; ///< Spatial point type

private:
    /**
     * @brief Initialize array of spatial points from a textual file.
     *
     * @param[in] filename File name to read from.
     * @param[out] pts Output array.
     * @return Number of points.
     */
    static uint64_t init_points(const std::string &filename, SmartArray<obj_type> &pts) {
        std::ifstream fin(filename, std::ios_base::in | std::ios_base::binary);
        uint32_t n_pts{};
        fin.read(reinterpret_cast<char *>(&n_pts), sizeof(uint32_t));
        pts.resize(n_pts);
        fin.read(reinterpret_cast<char *>(pts.data()), 3 * sizeof(DataType) * n_pts);
        return n_pts;
    }

    /**
     * @brief Field compoment of two points interaction.
     *
     * @param[in] idx_rcv Receiver point number.
     * @param[in] idx_src Source point number.
     * @param[in] i_local Resulting field component index (0, or 1, or 2).
     * @param[in] j_local Source field component index (0, or 1, or 2).
     * @return Field component value.
     */
    double field(uint64_t idx_rcv, uint64_t idx_src, uint64_t i_local, uint64_t j_local) const {
        auto radius = _pts_src[idx_src] - _pts_rcv[idx_rcv];
        double rho = radius.norm();
        radius *= 1.0 / rho;
        double result = 3 * radius[i_local] * radius[j_local];
        if (i_local == j_local) {
            result -= 1.0;
        }
        return result / (rho * rho * rho);
    }

    /**
     * @brief Field gradient compoment of two points interaction.
     *
     * @param[in] idx_rcv Receiver point number.
     * @param[in] idx_src Source point number.
     * @param[in] i_local Resulting gradient component index (from 0 to 4).
     * @param[in] j_local Source field component index (0, or 1, or 2).
     * @return Gradient component value.
     */
    double grad(uint64_t idx_rcv, uint64_t idx_src, uint64_t i_local, uint64_t j_local) const {
        auto radius = _pts_src[idx_src] - _pts_rcv[idx_rcv];
        double rho = radius.norm();
        //radius *= 1.0 / rho;
        double dx{radius[0]}, dy{radius[1]}, dz{radius[2]};
        double result;
        if (i_local + j_local == 3) {
            result = -5 * dx * dy * dz;
        } else if (i_local == 0 || (i_local < 3 && j_local == 0)) {
            double first{radius[i_local + j_local]};
            double second{rho * rho};
            if (i_local == 0 && j_local == 0) {
                second *= 3;
            }
            second -= 5 * dx * dx;
            result = first * second;
        } else if (i_local == 4 || (i_local > 1 && j_local == 2)) {
            double first{radius[i_local + j_local - 4]};
            double second{rho * rho};
            if (i_local == 4 && j_local == 2) {
                second *= 3;
            }
            second -= 5 * dz * dz;
            result = first * second;
        } else {
            double first{radius[i_local - 1]};
            double second{rho * rho - 5 * dy * dy};
            result = first * second;
        }
        return 3 * result / (rho * rho * rho * rho/**/ * rho * rho * rho/**/);
    }

public:
    /**
     * @brief Initialize interaction type and geometry.
     *
     * @param[in] type Type of interaction to use: 0 - field, 1 - gradient, 2 - field & gradient.
     * @param[in] rows Row (receiving) points filename to init from.
     * @param[in] cols Column (source) points filename to init from.
     */
    Magnetic(uint64_t type, const std::string &rows, const std::string &cols) :
            matrix::Matrix<DataType>(0, 0), _field{type != 1}, _grad{type > 0} {
        _n_rcv = init_points(rows, _pts_rcv);
        _n_src = init_points(cols, _pts_src);
        _b_rows = static_cast<uint64_t>(_field) * 3 + static_cast<uint64_t>(_grad) * 5;
        _b_cols = 3;

        this->n_rows = _b_rows * _n_rcv;
        this->n_cols = _b_cols * _n_src;
    }

    ~Magnetic() override = default;

    DataType element(uint64_t i, uint64_t j) const override {
        uint64_t idx_rcv{i / _b_rows}, idx_src{j / _b_cols};
        uint64_t i_local{i % _b_rows}, j_local{j % _b_cols};
        // uint64_t idx_rcv{i % _pts_rcv.size()}, idx_src{j % _pts_src.size()};
        // uint64_t i_local{i / _pts_rcv.size()}, j_local{j / _pts_src.size()};
        if (_field && _grad) {
            if (i_local < 3) { // field
                return field(idx_rcv, idx_src, i_local, j_local);
            } else { // grad
                return grad(idx_rcv, idx_src, i_local - 3, j_local);
            }
        } else if (_field) {
            return field(idx_rcv, idx_src, i_local, j_local);
        } else {
            return grad(idx_rcv, idx_src, i_local, j_local);
        }
    }

    /**
     * @brief Receiver point corresponding to a row index.
     *
     * @param[in] i Row index.
     * @return geom::Point.
     */
    const geom::Point<double, 3> &row_point(uint64_t i) const {
        uint64_t idx_rcv{i / _b_rows};
        // uint64_t idx_rcv{i % _pts_rcv.size()};
        return _pts_rcv[idx_rcv];
    }

    /**
     * @brief Source point corresponding to a column index.
     *
     * @param[in] j Column index.
     * @return geom::Point.
     */
    const geom::Point<double, 3> &col_point(uint64_t j) const {
        uint64_t idx_src{j / _b_cols};
        // uint64_t idx_src{j % _pts_src.size()};
        return _pts_src[idx_src];
    }
};


template<class DataType, class RefMatrix>
class SlicedMatrix: public matrix::Matrix<DataType> {
private:
    RefMatrix &_matrix;
    const std::vector<uint64_t> _rows;
    const std::vector<uint64_t> _cols;

public:
    using obj_type = typename RefMatrix::obj_type;
public:
    SlicedMatrix(RefMatrix &matrix, std::vector<uint64_t> &&rows, std::vector<uint64_t> &&cols) :
            matrix::Matrix<DataType>(rows.empty() ? matrix.size().first : rows.size(),
                                     cols.empty() ? matrix.size().second : cols.size()),
            _matrix{matrix}, _rows{std::move(rows)}, _cols{std::move(cols)} {}

    ~SlicedMatrix() override = default;

    DataType element(uint64_t i, uint64_t j) const override {
        uint64_t s_i = _rows.empty() ? i : _rows[i];
        uint64_t s_j = _cols.empty() ? j : _cols[j];
        return _matrix.element(s_i, s_j);
    }

    /**
     * @brief Receiver point corresponding to a row index.
     *
     * @param[in] i Row index.
     * @return geom::Point.
     */
    obj_type row_point(uint64_t i) const {
        uint64_t s_i = _rows.empty() ? i : _rows[i];
        return _matrix.row_point(s_i);
    }

    /**
     * @brief Source point corresponding to a column index.
     *
     * @param[in] j Column index.
     * @return geom::Point.
     */
    obj_type col_point(uint64_t j) const {
        uint64_t s_j = _cols.empty() ? j : _cols[j];
        return _matrix.col_point(s_j);
    }
};

}
#endif //MOSAIC_SKELETON_MATRIX_MAGNETIC_HPP
