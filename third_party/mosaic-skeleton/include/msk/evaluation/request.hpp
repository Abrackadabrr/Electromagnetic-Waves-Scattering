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
// Created by b.valiakhmetov on 23.03.23.
//


/**
 * @file
 *
 * @brief Defines class Request.
 */


#ifndef MOSAIC_SKELETON_EVALUATION_REQUEST_HPP
#define MOSAIC_SKELETON_EVALUATION_REQUEST_HPP


#include <algorithm>
#include <cstdint>
#include <complex>
#include <utility>
#include <vector>
#include <iostream>

#include "msk/index/index.hpp"


namespace MSk {

namespace eval {
/**
 * @brief Request of matrix elements for future evaluation.
 *
 * @tparam DataType Elements type.
 */
template<typename DataType>
class Request {
public:
    /**
     * @brief Request type
     */
    enum class Type {
        none, ///< Empty
        element, ///< Single element
        row, ///< Single row (not all columns possibly)
        column, ///< Single column (not all rows possibly)
        block, ///< Subblock
        list, ///< Arbitrary list
    };

private:
    Type _type{Type::none}; ///< Request type
    uint64_t _n_elems{0}; ///< Number of containing elements
    DataType *_ptr{nullptr}; ///< Pointer to write the result
    uint64_t _block_id{0}; ///< Identificator of a host Block
    std::vector<uint64_t> _row_idx; ///< Row indices
    std::vector<uint64_t> _col_idx; ///< Column indoces

private:
    /**
     * @brief Convert Type to integer for transfering reasons.
     *
     * @param[in] t Type value.
     * @return Integer equivalent.
     */
    uint8_t type_to_int(Type t) const {
        switch (t) {
            case Type::none:
                return 0;
            case Type::element:
                return 1;
            case Type::row:
                return 2;
            case Type::column:
                return 3;
            case Type::block:
                return 4;
            case Type::list:
                return 5;
            default:
                return 0;
        }
    }

    /**
     * @brief Convert integer to Type for transfering reasons.
     *
     * @param[in] i Integer value.
     * @return Type equivalent.
     */
    Type int_to_type(uint8_t i) const {
        switch (i) {
            case 0:
                return Type::none;
            case 1:
                return Type::element;
            case 2:
                return Type::row;
            case 3:
                return Type::column;
            case 4:
                return Type::block;
            case 5:
                return Type::list;
            default:
                return Type::none;
        }
    }

public:
    /**
     * @brief Empty request.
     */
    Request() : _type{Type::none}, _n_elems{0} {}

    ~Request() = default;

    Request(const Request &other) :
            _type{other._type}, _n_elems{other._n_elems}, _ptr{other._ptr},
            _block_id{other._block_id}, _row_idx{other._row_idx}, _col_idx{other._col_idx} {}

    Request &operator=(const Request &other) {
        if (&other == this) {
            return *this;
        }
        _type = other._type;
        _n_elems = other._n_elems;
        _ptr = other._ptr;
        _block_id = other._block_id;
        _row_idx = other._row_idx;
        _col_idx = other._col_idx;
        return *this;
    }

    Request(Request &&other) noexcept {
        _type = std::move(other._type);
        other._type = Type::none;
        _n_elems = std::move(other._n_elems);
        other._n_elems = 0;
        _ptr = std::move(other._ptr);
        other._ptr = nullptr;
        _block_id = std::move(other._block_id);
        other._block_id = 0;
        _row_idx = std::move(other._row_idx);
        other._row_idx.clear();
        _col_idx = std::move(other._col_idx);
        other._col_idx.clear();
    }

    Request &operator=(Request &&other) noexcept {
        if (&other == this) {
            return *this;
        }
        _type = std::move(other._type);
        other._type = Type::none;
        _n_elems = std::move(other._n_elems);
        other._n_elems = 0;
        _ptr = std::move(other._ptr);
        other._ptr = nullptr;
        _block_id = std::move(other._block_id);
        other._block_id = 0;
        _row_idx = std::move(other._row_idx);
        other._row_idx.clear();
        _col_idx = std::move(other._col_idx);
        other._col_idx.clear();
        return *this;
    }

    /**
     * @brief Get the id if the Block that requested elements.
     *
     * @return Block identificator.
     */
    uint64_t block_id() const {
        return _block_id;
    }
    /**
     * @brief Set the id if the Block that requested elements.
     *
     * @param[in] b_id Block identificator.
     */
    void set_block_id(uint64_t b_id) {
        _block_id = b_id;
    }

    /**
     * @brief Fill request with one element.
     *
     * @param[in] i Row index.
     * @param[in] j Column index.
     * @param[in] ptr Pointer to write the evaluated element.
     */
    void element(uint64_t i, uint64_t j, DataType *ptr) {
        _type = Type::element;
        _n_elems = 1;
        _row_idx.resize(1);
        _row_idx[0] = i;
        _col_idx.resize(1);
        _col_idx[0] = j;
        _ptr = ptr;
    }

    /**
     * @brief Fill request with one row (not all columns possibly).
     *
     * @tparam Idx Columns container type.
     * @param[in] r Row index.
     * @param[in] columns Container of columns indices that supports iteration.
     * @param[in] ptr Pointer to write the evaluated elements.
     */
    template<class Idx>
    void row(uint64_t r, const Idx &columns, DataType *ptr) {
        _type = Type::row;
        _n_elems = columns.size();
        _row_idx.resize(1);
        _row_idx[0] = r;
        _col_idx.resize(_n_elems);
        uint64_t j{0};
        for (auto it = columns.begin(); it != columns.end(); ++it, ++j) {
            _col_idx[j] = *it;
        }
        _ptr = ptr;
    }

    /**
     * @brief Fill request with one column (not all rows possibly).
     *
     * @tparam Idx Rows container type.
     * @param[in] rows Container of rows indices that supports iteration.
     * @param[in] c Column index.
     * @param[in] ptr Pointer to write the evaluated elements.
     */
    template<class Idx>
    void column(const Idx &rows, uint64_t c, DataType *ptr) {
        _type = Type::column;
        _n_elems = rows.size();
        _col_idx.resize(1);
        _col_idx[0] = c;
        _row_idx.resize(_n_elems);
        uint64_t i{0};
        for (auto it = rows.begin(); it != rows.end(); ++it, ++i) {
            _row_idx[i] = *it;
        }
        _ptr = ptr;
    }

    /**
     * @brief Fill request with matrix subblock.
     *
     * @tparam Idx Rows and columns container type.
     * @param[in] rows Container of block rows indices that supports iteration.
     * @param[in] columns Container of block columns indices that supports iteration.
     * @param[in] ptr Pointer to write the evaluated elements.
     */
    template<class Idx>
    void block(const Idx &rows, const Idx &columns, DataType *ptr) {
        _type = Type::block;
        uint64_t m{rows.size()}, n{columns.size()};
        _n_elems = m * n;
        _row_idx.resize(m);
        _col_idx.resize(n);
        uint64_t i{0};
        for (auto it = rows.begin(); it != rows.end(); ++it, ++i) {
            _row_idx[i] = *it;
        }
        uint64_t j{0};
        for (auto it = columns.begin(); it != columns.end(); ++it, ++j) {
            _col_idx[j] = *it;
        }
        _ptr = ptr;
    }

    /**
     * @brief Fill request with an arbitrary list of matrix elements.
     *
     * @tparam Idx Rows and columns container type.
     * @param[in] rows Container of rows indices that supports iteration.
     * @param[in] columns Container of columns indices that supports iteration.
     * @param[in] ptr Pointer to write the evaluated elements.
     */
    template<class Idx>
    void list(const Idx &rows, const Idx &columns, DataType *ptr) {
        _type = Type::list;
        _n_elems = std::min(rows.size(), columns.size());
        _row_idx.resize(_n_elems);
        _col_idx.resize(_n_elems);
        uint64_t i{0};
        for (auto it = rows.begin(); it != rows.end() && i < _n_elems; ++it, ++i) {
            _row_idx[i] = *it;
        }
        uint64_t j{0};
        for (auto it = columns.begin(); it != columns.end() && j < _n_elems; ++it, ++j) {
            _col_idx[j] = *it;
        }
        _ptr = ptr;
    }

    Type type() const {
        return _type;
    }

    uint64_t size() const {
        return _n_elems;
    }

    bool empty() const {
        return _n_elems == 0;
    }

    /**
     * @brief Pointer to put evaluated elements.
     *
     * @return Pointer.
     */
    DataType *result_ptr() const {
        return _ptr;
    }

    /**
     * @brief Set the pointer for putting evaluated elements.
     *
     * @param[in] ptr Pointer.
     */
    void set_result_ptr(DataType *ptr) {
        _ptr = ptr;
    }

    /**
     * @brief Redefine block-related local indices if elements to be evaluated
     * according to matrix-related global indices.
     * @param[in] row_idx Rows indices.
     * @param[in] col_idx Columns indices.
     */
    void at_index(const index::Index<uint64_t> &row_idx, const index::Index<uint64_t> &col_idx) {
        for (auto &r: _row_idx) {
            r = row_idx[r];
        }
        for (auto &c: _col_idx) {
            c = col_idx[c];
        }
    }

    /**
     * @brief Extract elements' indices as lists of rows and columns indices.
     *
     * @param[out] rows Vector of rows indices.
     * @param[out] columns Vector of columns indices.
     * @return Number of elements to evaluate.
     */
    uint64_t extract_list(std::vector<uint64_t> &rows, std::vector<uint64_t> &columns) const {
        rows.resize(_n_elems);
        columns.resize(_n_elems);
        if (_type == Type::element) {
            rows[0] = _row_idx[0];
            columns[0] = _col_idx[0];
        } else if (_type == Type::row) {
            for (uint64_t i = 0; i < _n_elems; ++i) {
                rows[i] = _row_idx[0];
            }
            std::copy(_col_idx.cbegin(), _col_idx.cend(), columns.begin());
        } else if (_type == Type::column) {
            std::copy(_row_idx.cbegin(), _row_idx.cend(), rows.begin());
            for (uint64_t i = 0; i < _n_elems; ++i) {
                columns[i] = _col_idx[0];
            }
        } else if (_type == Type::block) {
            const uint64_t n_rows{_row_idx.size()}, n_cols{_col_idx.size()};
            for (uint64_t j = 0; j < n_cols; ++j) {
                std::copy(_row_idx.cbegin(), _row_idx.cend(), rows.begin() + j * n_rows);
                for (uint64_t i = 0; i < n_rows; ++i) {
                    columns[i + n_rows * j] = _col_idx[j];
                }
            }
        } else if (_type == Type::list) {
            std::copy(_row_idx.cbegin(), _row_idx.cend(), rows.begin());
            std::copy(_col_idx.cbegin(), _col_idx.cend(), columns.begin());
        }
        return _n_elems;
    }

    /**
     * @brief Extract elements' indices in the way they are stored inside Request.
     *
     * @param[out] rows Vector of rows indices.
     * @param[out] columns Vector of columns indices.
     * @return Number of elements to evaluate.
     */
    uint64_t extract_compact(std::vector<uint64_t> &rows, std::vector<uint64_t> &columns) const {
        rows = _row_idx;
        columns = _col_idx;
        return _n_elems;
    }

    /**
     * @brief Split the request into two (cut from the end).
     *
     * @details The procedure tries to preserve type of the request if possible.
     *
     * @param[in] n_split Number of elements to cut.
     * @return New request.
     */
    Request split(uint64_t n_split) {
        if (n_split >= _n_elems || _type == Type::none || _type == Type::element) {
            return std::move(*this);
        }

        Request result;
        const uint64_t n_remain = _n_elems - n_split;
        result._n_elems = n_split;
        result._ptr = _ptr + n_remain;
        result._block_id = _block_id;

        if (_type == Type::row) {
            result._type = _type;
            result._row_idx.resize(1, _row_idx[0]);
            result._col_idx = std::vector<uint64_t>(_col_idx.cbegin() + n_remain, _col_idx.cend());
            _col_idx.resize(n_remain);
        } else if (_type == Type::column) {
            result._type = _type;
            result._col_idx.resize(1, _col_idx[0]);
            result._row_idx = std::vector<uint64_t>(_row_idx.cbegin() + n_remain, _row_idx.cend());
            _row_idx.resize(n_remain);
        } else if (_type == Type::block) {
            if (n_split % _row_idx.size() == 0) { // Cut several full columns if possible
                result._type = _type;
                result._row_idx = _row_idx;
                const uint64_t remain_cols = _col_idx.size() - n_split / _row_idx.size();
                result._col_idx = std::vector<uint64_t>(_col_idx.cbegin() + remain_cols, _col_idx.cend());
                _col_idx.resize(remain_cols);
            } else { // Otherwise restore the request as List
                std::vector<uint64_t> new_rows, new_cols;
                extract_list(new_rows, new_cols);
                _row_idx = std::move(new_rows);
                _col_idx = std::move(new_cols);
                _type = Type::list;
            }
        }
        if (_type == Type::list) {
            result._type = _type;
            result._row_idx.insert(result._row_idx.cbegin(), _row_idx.begin() + n_remain, _row_idx.end());
            result._col_idx.insert(result._col_idx.cbegin(), _col_idx.begin() + n_remain, _col_idx.end());
            _row_idx.resize(n_remain);
            _col_idx.resize(n_remain);
        }
        _n_elems = n_remain;
        return result;
    }

    /**
     * @brief Clear the request.
     */
    void free() {
        _type = Type::none;
        _n_elems = 0;
        _row_idx.clear();
        _col_idx.clear();
        _ptr = nullptr;
    }

    /**
     * @brief Size of the request packed into bytes.
     *
     * @return Number of bytes.
     */
    uint64_t packed_size() const {
        uint64_t sz = 1 + sizeof(_n_elems) + sizeof(_ptr) + sizeof(_block_id);
        sz += 2 * sizeof(uint64_t) + _row_idx.size() * sizeof(_row_idx[0]) + _col_idx.size() * sizeof(_col_idx[0]);
        return sz;
    }

    /**
     * @brief Pack the request in bytes.
     *
     * @param[out] ptr Array to pack to.
     * @return Packed size (number of bytes).
     */
    uint64_t pack(uint8_t *ptr) {
        // Pack metadata
        ptr[0] = type_to_int(_type);
        ptr += 1;
        reinterpret_cast<uint64_t *>(ptr)[0] = _n_elems;
        ptr += sizeof(_n_elems);
        reinterpret_cast<DataType **>(ptr)[0] = _ptr;
        ptr += sizeof(_ptr);
        reinterpret_cast<uint64_t *>(ptr)[0] = _block_id;
        ptr += sizeof(_block_id);

        // Pack rows
        uint64_t n_rows = _row_idx.size();
        reinterpret_cast<uint64_t *>(ptr)[0] = n_rows;
        ptr += sizeof(n_rows);
        std::memcpy(ptr, _row_idx.data(), n_rows * sizeof(_row_idx[0]));
        ptr += _row_idx.size() * sizeof(_row_idx[0]);

        // Pack columns
        uint64_t n_cols = _col_idx.size();
        reinterpret_cast<uint64_t *>(ptr)[0] = n_cols;
        ptr += sizeof(n_cols);
        std::memcpy(ptr, _col_idx.data(), n_cols * sizeof(_col_idx[0]));
        return packed_size();
    }

    /**
     * @brief Unpack the request from bytes sequence.
     *
     * @param[in] ptr Array of bytes.
     * @return Packed size (number of bytes).
     */
    uint64_t unpack(const uint8_t *ptr) {
        // Unpack metadata
        _type = int_to_type(ptr[0]);
        ptr += 1;
        _n_elems = reinterpret_cast<const uint64_t *>(ptr)[0];
        ptr += sizeof(_n_elems);
        _ptr = reinterpret_cast<DataType * const*>(ptr)[0];
        ptr += sizeof(_ptr);
        _block_id = reinterpret_cast<const uint64_t *>(ptr)[0];
        ptr += sizeof(_block_id);

        // Unpack rows
        uint64_t n_rows = reinterpret_cast<const uint64_t *>(ptr)[0];
        ptr += sizeof(n_rows);
        _row_idx.resize(n_rows);
        std::memcpy(_row_idx.data(), ptr, n_rows * sizeof(_row_idx[0]));
        ptr += _row_idx.size() * sizeof(_row_idx[0]);

        // Unpack columns
        uint64_t n_cols = reinterpret_cast<const uint64_t *>(ptr)[0];
        ptr += sizeof(n_cols);
        _col_idx.resize(n_cols);
        std::memcpy(_col_idx.data(), ptr, n_cols * sizeof(_col_idx[0]));
        return packed_size();
    }
};


extern template class Request<float>;

extern template class Request<double>;

extern template class Request<std::complex<float>>;

extern template class Request<std::complex<double>>;


}
}

#endif //MOSAIC_SKELETON_EVALUATION_REQUEST_HPP
