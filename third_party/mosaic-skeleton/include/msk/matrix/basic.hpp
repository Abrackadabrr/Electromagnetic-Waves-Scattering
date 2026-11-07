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
// Created by b.valiakhmetov on 9/9/22.
//


/**
 * @file
 *
 * @brief Contains some examples of matrix functions.
 */


#ifndef MOSAIC_SKELETON_MATRIX_EXAMPLES_HPP
#define MOSAIC_SKELETON_MATRIX_EXAMPLES_HPP


#include <cmath>
#include <complex>
#include <limits>

#include "msk/matrix/matrix.hpp"
#include "msk/geometry/point.hpp"



namespace MSk {

namespace matrix {


/**
 * @brief Basic functional matrix of integer arguments.
 *
 * @tparam DataType Element type.
 */
template<typename DataType>
class Elementary : public Matrix<DataType> {
private:
    std::function<DataType(uint64_t, uint64_t)> _function; ///< Element function

public:
    using obj_type = double;
public:
    /**
     * @brief Cunstruct matrix from fucntional object.
     *
     * @param[in] m Rows number.
     * @param[in] n Columns number.
     * @param[in] element_function Function of element evaluation.
     */
    Elementary(uint64_t m, uint64_t n,
               std::function<DataType(uint64_t, uint64_t)> element_function) :
               Matrix<DataType>(m, n), _function{element_function} {}

    DataType element(uint64_t i, uint64_t j) const override {
        return _function(i, j);
    }

    /**
     * @brief Space point corresponding to a row.
     *
     * @param[in] i Row index.
     * @return Point.
     */
    obj_type row_point(uint64_t i) const {
        return (1.0 * i) / this->n_rows;
    }

    /**
     * @brief Space point corresponding to a column.
     *
     * @param[in] j Column index.
     * @return Point.
     */
    obj_type col_point(uint64_t j) const {
        return (1.0 * j) / this->n_cols;
    }
};


/**
 * @brief Functional matrix of two points in space.
 *
 * @tparam DataType Element type.
 * @tparam PointType Point coordinate type.
 * @tparam dim Point dimension.
 */
template<class DataType, typename PointType, uint64_t dim = 3>
class BiPoint : public Matrix<DataType> {
public:
    using obj_type = geom::Point<PointType, dim>; ///< Spatial point type

private:
    std::function<DataType(obj_type, obj_type)> _function; ///< Element function
    const geom::PointGrid<PointType, dim> &rows, &cols; ///< Points corresponding to rows and column

public:
    /**
     * @brief Initialize the matrix with element fucntion spatial points.
     *
     * @param[in] element_function Element evaluation function.
     * @param[in] row_pts Points corresponding to rows indoces.
     * @param[in] col_pts Points corresponding to columns indoces.
     */
    BiPoint(std::function<DataType(const obj_type &, const obj_type &)> element_function,
            const geom::PointGrid<PointType, dim> &row_pts, const geom::PointGrid<PointType, dim> &col_pts) :
            Matrix<DataType>{row_pts.size(), col_pts.size()},
            _function{element_function}, rows{row_pts}, cols{col_pts} {}

    DataType element(uint64_t i, uint64_t j) const override {
        return _function(rows[i], cols[j]);
    }

    /**
     * @brief Space point corresponding to a row.
     *
     * @param[in] i Row index.
     * @return geom::Point.
     */
    const obj_type &row_point(uint64_t i) const {
        return rows[i];
    }

    /**
     * @brief Space point corresponding to a column.
     *
     * @param[in] j Column index.
     * @return geom::Point.
     */
    const obj_type &col_point(uint64_t j) const {
        return cols[j];
    }
};


/**
 * @brief Dense (dummy functional) matrix with precomputed elements corresponding to spatial points.
 *
 * @tparam DataType Matrix element type.
 * @tparam dim Space dimension.
 */
template<class DataType, uint64_t dim = 3>
class Dense: public Matrix<DataType> {
public:
    /**
     * @brief Dense matrix type.
     */
    enum class Type {
        rect, ///< General ectangular
        sym_up, ///< Symmetric stored as upper triangular
        sym_low, ///< Symmetric stored as lower triangular
        tri_up, ///< Upper triangular
        tri_low ///< Lower triangular
    };

    /**
     * @brief Elements storing mode.
     */
    enum class StoreMode {
        rows, ///< Rowwize
        columns ///< Columnwise
    };

private:
    using RealType = decltype(std::abs(DataType{}));
    Type type;
    StoreMode mode;
    SmartArray<DataType> data; ///< Matrix elements
    geom::PointGrid<RealType, dim> row_points, col_points; ///< Corresponding points

public:
    using obj_type = geom::Point<RealType, dim>;

private:
    /**
     * @brief Init matrix elements from a binary file.
     *
     * @param[in] fin Open file object.
     */
    void init_bin(std::ifstream &fin) {
        if (type == Type::rect) {
            uint32_t m{}, n{};
            fin.read(reinterpret_cast<char *>(&m), sizeof(uint32_t));
            fin.read(reinterpret_cast<char *>(&n), sizeof(uint32_t));
            this->n_rows = m;
            this->n_cols = n;
            data.resize(m * n);
            fin.read(reinterpret_cast<char *>(data.data()), sizeof(DataType) * m * n);
        } else {
            uint32_t n{};
            fin.read(reinterpret_cast<char *>(&n), sizeof(uint32_t));
            this->n_rows = this->n_cols = n;
            uint64_t n_elements{(n * (n + 1)) / 2};
            data.resize(n_elements);
            fin.read(reinterpret_cast<char *>(data.data()), sizeof(DataType) * n_elements);
        }
    }

    /**
     * @brief Skip comment line starting with '#' symbol in a textual file.
     *
     * @param[in] fin Open file object.
     */
    static void skip_comments(std::ifstream &fin) {
        for (auto c = fin.peek(); c == '#'; c = fin.peek()) {
            fin.ignore(std::numeric_limits<std::streamsize>::max(), fin.widen('\n'));
        }
    }

    /**
     * @brief Init matrix elements from a textual file.
     *
     * @param[in] fin Open file object.
     */
    void init_txt(std::ifstream &fin) {
        skip_comments(fin);
        uint64_t n_elements{};
        if (type == Type::rect) {
            uint32_t m{}, n{};
            fin >> m >> n;
            this->n_rows = m;
            this->n_cols = n;
            n_elements = m * n;
        } else {
            uint32_t n{};
            fin >> n;
            this->n_rows = this->n_cols = n;
            n_elements = (n * (n + 1)) / 2;
        }
        data.resize(n_elements);
        for (uint64_t k = 0; k < n_elements; ++k) {
            fin >> data[k];
        }
    }


public:
    /**
     * @brief Initialize dense matrix from files.
     *
     * @details If matrix is rectangular then file should contain its row and column sizes
     * of type uint32_t, then its elements. Otherwise (triangular or symmetric) a file
     * should start with a single size of type uint32_t, and then the elements should follow.
     *
     * @param[in] data_fname Elements file name.
     * @param[in] type Matrix type.
     * @param[in] mode Element storing mode.
     * @param[in] rows_fname Rows points file name.
     * @param[in] cols_fname Columns points file name.
     */
    Dense(const std::string &data_fname, Type type = Type::rect, StoreMode mode = StoreMode::columns,
          const std::string &rows_fname = "", const std::string &cols_fname = "") :
            Matrix<DataType>(0, 0), type{type}, mode{mode} {
        if (data_fname.size() > 4 && data_fname.compare(data_fname.size() - 4, 4, ".bin") == 0) {
            std::ifstream fin(data_fname, std::ios_base::in | std::ios_base::binary);
            init_bin(fin);
        } else {
            std::ifstream fin(data_fname, std::ios_base::in);
            init_txt(fin);
        }

        if (rows_fname.empty()) {
            row_points = std::move(geom::PointGrid<RealType, dim>{this->n_rows});
        } else {
            row_points = std::move(geom::PointGrid<RealType, dim>{rows_fname});
        }
        if (cols_fname == rows_fname) {
            col_points = row_points;
        } else if (cols_fname.empty()) {
            col_points = std::move(geom::PointGrid<RealType, dim>{this->n_cols});
        } else {
            col_points = std::move(geom::PointGrid<RealType, dim>{cols_fname});
        }
    }

    ~Dense() override = default;

    DataType element(uint64_t i, uint64_t j) const override {
        DataType val{};
        uint64_t m{this->n_rows}, n{this->n_cols};
        if (type == Type::rect) {
            if (mode == StoreMode::rows) {
                val = data[i * n + j];
            } else if (mode == StoreMode::columns) {
                val = data[i + j * m];
            }
        } else {
            if ((type == Type::tri_up && i > j) || (type == Type::tri_low && i < j)) {
                val = 0;
            } else {
                if (i > j) {
                    std::swap(i, j);
                }
                if ((mode == StoreMode::rows && (type == Type::sym_low || type == Type::tri_low)) ||
                    (mode == StoreMode::columns && (type == Type::sym_up || type == Type::tri_up))) {
                    val = data[i + (j * (j + 1)) / 2];
                    } else if ((mode == StoreMode::rows && (type == Type::sym_up || type == Type::tri_up)) ||
                               (mode == StoreMode::columns && (type == Type::sym_low || type == Type::tri_low))) {
                        val = data[(i * (2 * n - i + 1)) / 2 + (j - i)];
                               }
            }
        }
        return val;
    }

    /**
     * @brief Space point corresponding to a row.
     *
     * @param[in] i Row index.
     * @return geom::Point.
     */
    const obj_type &row_point(uint64_t i) const {
        return row_points[i];
    }

    /**
     * @brief Space point corresponding to a column.
     *
     * @param[in] j Column index.
     * @return geom::Point.
     */
    const obj_type &col_point(uint64_t j) const {
        return col_points[j];
    }
};


/**
 * @brief Bilinear functional matrix.
 *
 * @details Element value is Ax + By + C + Dxy, where x and y are row and column indices, respectively.
 *
 * @tparam DataType Matrix element type.
 */
template<class DataType>
class BiLinear : public Matrix<DataType> {
private:
    DataType _c0, _c1, _c2, _c12; ///< Coefficients.

public:
    /**
     * @brief Initialize function coefficients.
     *
     * @param[in] m Number of matrix rows.
     * @param[in] n Number of matrix columns.
     * @param[in] cx First linear coefficient.
     * @param[in] cy Second linear coefficient.
     * @param[in] c Free coefficient.
     * @param[in] cxy Bilinear coefficient.
     */
    BiLinear(uint64_t m, uint64_t n, 
            DataType cx = DataType{1.0}, DataType cy = DataType{1.0}, 
            DataType c = DataType{0.0}, DataType cxy = DataType{0.0}) :
            matrix::Matrix<DataType>{m, n}, _c0{c}, _c1{cx}, _c2{cy}, _c12{cxy} {}

    DataType element(uint64_t i, uint64_t j) const override {
        DataType x = i, y = j;
        return _c0 + _c1 * x + _c2 * y + _c12 * x * y;
    }

    /**
     * @brief Functor call type element getter.
     *
     * @param[in] i Row index.
     * @param[in] j Columns index.
     * @return Matrix element.
     */
    DataType operator()(uint64_t i, uint64_t j) const {
        return element(i, j);
    }
};


extern template class Elementary<float>;

extern template class Elementary<double>;

extern template class Elementary<std::complex<float>>;

extern template class Elementary<std::complex<double>>;


extern template class BiPoint<float, float, 2>;

extern template class BiPoint<double, double, 2>;

extern template class BiPoint<std::complex<float>, float, 2>;

extern template class BiPoint<std::complex<double>, double, 2>;


extern template class BiPoint<float, float, 3>;

extern template class BiPoint<double, double, 3>;

extern template class BiPoint<std::complex<float>, float, 3>;

extern template class BiPoint<std::complex<double>, double, 3>;


extern template class Dense<float, 2>;

extern template class Dense<double, 2>;

extern template class Dense<std::complex<float>, 2>;

extern template class Dense<std::complex<double>, 2>;


extern template class Dense<float, 3>;

extern template class Dense<double, 3>;

extern template class Dense<std::complex<float>, 3>;

extern template class Dense<std::complex<double>, 3>;


extern template class BiLinear<float>;

extern template class BiLinear<double>;

extern template class BiLinear<std::complex<float>>;

extern template class BiLinear<std::complex<double>>;


}
}

#endif //MOSAIC_SKELETON_MATRIX_EXAMPLES_HPP
