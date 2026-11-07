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
// created by bulat_v on 11.02.2022.
//


/**
 * @file
 *
 * @brief Defines basic classes representing functional matrices.
 */


#ifndef MOSAIC_SKELETON_MATRIX_MATRIX_HPP
#define MOSAIC_SKELETON_MATRIX_MATRIX_HPP


#include <omp.h>

#include <utility>
#include <functional>

#include "msk/geometry/point.hpp"
#include "msk/geometry/point_grid.hpp"


namespace MSk {

namespace matrix {


/**
 * @brief Base functional matrix class.
 *
 * @details Matrix represents the function of its arbitrary element evaluation.
 *
 * @tparam DataType Element type.
 */
template<typename DataType>
class Matrix {
protected:
    uint64_t n_rows{0}, n_cols{0}; ///< Matrix sizes

public:
    /**
     * @brief Abstract matrix of the given shape.
     *
     * @param[in] num_rows Number of rows.
     * @param[in] num_cols Number of columns.
     */
    Matrix(uint64_t num_rows, uint64_t num_cols) :
           n_rows{num_rows}, n_cols{num_cols}{}

    virtual ~Matrix() = default;

    /**
     * @brief Evaluate an element of the matrix.
     *
     * @param[in] i Row index.
     * @param[in] j Column index.
     * @return Matrix element value.
     */
    virtual DataType element(uint64_t i, uint64_t j) const = 0;

    /**
     * @brief Evaluate a row of the matrix.
     *
     * @tparam Idx Columns container type.
     * @param[in] i Row index.
     * @param[in] col_J Container of columns indices that supports iteration.
     * @param[out] arr Array to write the evaluated elements.
     */
    template<class Idx>
    void row(uint64_t i, const Idx &col_J, DataType *arr) const {
        if (arr == nullptr) {
            return;
        }
        /*
        for (auto it_j = col_J.begin(); it_j != col_J.end(); ++it_j, ++arr) {
            *arr = element(i, *it_j);
        }
        */

        uint64_t n{col_J.size()};
        auto it_j = col_J.begin();
#pragma omp parallel for
        for (uint64_t j = 0; j < n; ++j) {
            arr[j] = element(i, *(it_j + j));
        }
    }

    /**
     * @brief Evaluate a column of the matrix.
     *
     * @tparam Idx Columns container type.
     * @param[in] row_I Container of rows indices that supports iteration.
     * @param[in] j Column index.
     * @param[out] arr Array to write the evaluated elements.
     */
    template<class Idx>
    void column(const Idx &row_I, uint64_t j, DataType *arr) const {
        if (arr == nullptr) {
            return;
        }
        /*
        for (auto it_i = row_I.begin(); it_i != row_I.end(); ++it_i, ++arr) {
            *arr = element(*it_i, j);
        }
        */
        uint64_t m{row_I.size()};
        auto it_i = row_I.begin();
#pragma omp parallel for
        for (uint64_t i = 0; i < m; ++i) {
            arr[i] = element(*(it_i + i), j);
        }
    }

    /**
     * @brief Evaluate a subblock of the matrix.
     *
     * @tparam Idx Columns container type.
     * @param[in] row_I Container of rows indices that supports iteration.
     * @param[in] col_J Container of columns indices that supports iteration.
     * @param[out] arr Array to write the evaluated elements.
     */
    template<class Idx>
    void block(const Idx &row_I, const Idx &col_J, DataType *arr) const {
        if (arr == nullptr) {
            return;
        }
        /*
        for (auto it_j = col_J.begin(); it_j != col_J.end(); ++it_j) {
            for (auto it_i = row_I.begin(); it_i != row_I.end(); ++it_i, ++arr) {
                *arr = element(*it_i, *it_j);
            }
        }
        */
        uint64_t m{row_I.size()}, n{col_J.size()};
        auto it_i = row_I.begin(), it_j = col_J.begin();
#pragma omp parallel for collapse(2)
        for (uint64_t j = 0; j < n; ++j) {
            for (uint64_t i = 0; i < m; ++i) {
                arr[i + m * j] = element(*(it_i + i), *(it_j + j));
            }
        }

    }

    /**
     * @brief Evaluate an arbitrary list of the matrix elements.
     *
     * @tparam Idx Columns container type.
     * @param[in] row_I Container of rows indices that supports iteration.
     * @param[in] col_J Container of columns indices that supports iteration.
     * @param[out] arr Array to write the evaluated elements.
     */
    template<class Idx>
    void list(const Idx &row_I, const Idx &col_J, DataType *arr) const {
        if (arr == nullptr) {
            return;
        }
        if (row_I.size() != col_J.size()) {
            return;
        }
        /*
        auto it_i = row_I.begin();
        auto it_j = col_J.begin();
        for (; it_i != row_I.end(); ++it_i, ++it_j, ++arr) {
            *arr = element(*it_i, *it_j);
        }
        */
        uint64_t n{col_J.size()};
        auto it_i = row_I.begin(), it_j = col_J.begin();
#pragma omp parallel for
        for (uint64_t k = 0; k < n; ++k) {
            arr[k] = element(*(it_i + k), *(it_j + k));
        }
    }

    /**
     * @brief Matrix size.
     *
     * @return Rows number, columns number.
     */
    std::pair<uint64_t, uint64_t> size() const {
        return std::make_pair(n_rows, n_cols);
    }

};


extern template class Matrix<float>;

extern template class Matrix<double>;

extern template class Matrix<std::complex<float>>;

extern template class Matrix<std::complex<double>>;


}
} //namespace MSk


#endif //MOSAIC_SKELETON_MATRIX_MATRIX_HPP

