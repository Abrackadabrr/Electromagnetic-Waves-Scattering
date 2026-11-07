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
// Created by bulat_v on 22.01.24.
//


/**
 * @file
 *
 * @brief Contains template classes for Toeplitz and circulant matrix data format.
 */


#ifndef MOSAIC_SKELETON_BLOCK_TOEPLITZ_HPP
#define MOSAIC_SKELETON_BLOCK_TOEPLITZ_HPP


#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"

#include "msk/block/data.hpp"


namespace MSk {

namespace block {


/**
 * @brief Data representation for diagonal block format.
 *
 * @tparam DataType Block element type.
 */
template<class DataType>
class Toeplitz : public Data<DataType> {
private:
    using RealType = decltype(std::abs(DataType()));
    uint64_t _memory{0};
    double _norm{0};
    SmartArray<DataType> _data;

private:
    /**
     * @brief Number of bytes required to store and recover the representation.
     *
     * @return uint64_t bytes.
     */
    uint64_t _local_packed_size() const {
        return sizeof(uint64_t) + sizeof(double) + _memory * sizeof(DataType);
    }

public:
    Toeplitz(uint64_t m, uint64_t n) : Data<DataType>(m, n) {}

    Toeplitz(const Toeplitz &) = default;

    Toeplitz(Toeplitz &&other) noexcept : Data<DataType>(std::move(other)) {
        _memory = other._memory;
        other._memory = 0;
        _norm = other._norm;
        other._norm = 0.0;
        _data = std::move(other._data);
    }

    Toeplitz &operator=(const Toeplitz &) = default;

    Toeplitz &operator=(Toeplitz &&other) noexcept {
        Data<DataType>::operator=(std::move(other));
        _memory = other._memory;
        other._memory = 0;
        _norm = other._norm;
        other._norm = 0.0;
        _data = std::move(other._data);
        return *this;
    }

    ~Toeplitz() override = default;

    /**
     * @brief Copy block elements from an array.
     *
     * @param[in] elements SmartArray of elements, starting from first column's last to first row's last.
     * @param[in] norm Precomputed norm of the block.
     *
     * @note If norm is not provided (equals zero) then it will be computed.
     */
    void set_data(const SmartArray<DataType> &elements, double norm = 0.0) {
        uint64_t sz{std::min(this->n_rows, this->n_cols)};
        _data = elements;
        _data.resize(sz);
        if (norm == 0.0) {
            _norm = BLAS::nrm2(sz, _data.data(), 1);
        } else {
            _norm = norm;
        }
        _memory = sz;
    }

    /**
     * @brief Move block elements from an array.
     *
     * @param[in, out] elements SmartArray of elements, starting from first column's last to first row's last.
     * @param[in] norm Precomputed norm of the block.
     *
     * @note If norm is not provided (equals zero) then it will be computed.
     */
    void set_data(SmartArray<DataType> &&elements, double norm = 0.0) {
        uint64_t sz{std::min(this->n_rows, this->n_cols)};
        _data = std::move(elements);
        _data.resize(sz);
        if (norm == 0.0) {
            _norm = 0;
            // for (uint64_t i = 0; i < )
            // _norm = BLAS::nrm2(sz, _data.data(), 1);
        } else {
            _norm = norm;
        }
        _memory = sz;
    }

    uint64_t memory() const override {
        return _memory;
    }

    RealType frob_norm() const override {
        return _norm;
    }

    bool exact() const override {
        return true;
    }

    void matvec(const DataType *x, DataType *y) const override {
        if (_data.size() == 0) {
            return;
        }
        for (uint64_t i = 0; i < std::min(this->n_rows, this->n_cols); ++i) {
            y[i] += _data[i] * x[i];
        }
    }

    void transp_matvec(const DataType *x, DataType *y) const override {
        if (_data.size() == 0) {
            return;
        }
        for (uint64_t i = 0; i < std::min(this->n_rows, this->n_cols); ++i) {
            y[i] += _data[i] * x[i];
        }
    }

    void convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size) const override {
        if (_data.size() == 0) {
            return;
        }
        const uint64_t m{this->n_rows}, n{this->n_cols};
        if (size == 0) {
            size = m + n - 1;
        }
        size = std::min(size, m + n - 1);
        for (uint64_t i = 0; i < std::min(this->n_rows, this->n_cols) && 2 * i < size; ++i) {
            z[2 * i] += y[i] * _data[i] * x[i];
        }
    }

    void conj_convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size) const override {
        if (_data.size() == 0) {
            return;
        }
        const uint64_t m{this->n_rows}, n{this->n_cols};
        if (size == 0) {
            size = m + n - 1;
        }
        size = std::min(size, m + n - 1);
        for (uint64_t i = 0; i < std::min(this->n_rows, this->n_cols) && 2 * i < size; ++i) {
            z[2 * i] += y[i] * conjugate(_data[i]) * x[i];
        }
    }

    void to_dense(DataType *ptr, uint64_t ldim) const override {
        const DataType zero{0};
        if (ldim == this->n_rows) {
            BLAS::copy(this->n_rows * this->n_cols, &zero, 0, ptr, 1);
        } else {
            for (uint64_t j = 0; j < this->n_cols; ++j) {
                BLAS::copy(this->n_rows, &zero, 0, ptr + j * ldim, 1);
            }
        }
        BLAS::copy(std::min(this->n_rows, this->n_cols), _data.data(), 1, ptr, ldim + 1);
    }

    uint64_t pack(uint8_t *ptr) const override {
        uint64_t data_off{Data<DataType>::pack(ptr)};
        if (ptr != nullptr) {
            ptr += data_off;
            reinterpret_cast<uint64_t *>(ptr)[0] = _memory;
            ptr += sizeof(uint64_t);
            reinterpret_cast<double *>(ptr)[0] = _norm;
            ptr += sizeof(double);
            std::memcpy(ptr, _data.data(), _memory * sizeof(DataType));
        }
        return data_off + _local_packed_size();
    }

    uint64_t unpack(uint8_t *ptr) override {
        uint64_t data_off{Data<DataType>::unpack(ptr)};
        if (ptr != nullptr) {
            ptr += data_off;
            _memory = reinterpret_cast<uint64_t *>(ptr)[0];
            ptr += sizeof(uint64_t);
            _norm = reinterpret_cast<double *>(ptr)[0];
            ptr += sizeof(double);
            _data.resize(_memory);
            std::memcpy(_data.data(), ptr, _memory * sizeof(DataType));
        }
        return data_off + _local_packed_size();
    }

    void free() override {
        Data<DataType>::free();
        _data.free();
        _norm = 0;
        _memory = 0;
    }

    comp::Result recompress(comp::Parameters recomp_params) override {
        (void) recomp_params;
        return {_norm, 0, 1, true};
    }
};
}
}


#endif //MOSAIC_SKELETON_BLOCK_TOEPLITZ_HPP
