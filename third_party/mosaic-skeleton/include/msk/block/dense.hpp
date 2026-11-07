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
// Created by b.valiakhmetov on 16.02.23.
//


/**
 * @file
 *
 * @brief Contains template class block::Dense, which represents dense block data format.
 */


#ifndef MOSAIC_SKELETON_BLOCK_DENSE_HPP
#define MOSAIC_SKELETON_BLOCK_DENSE_HPP


#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"

#include "msk/block/data.hpp"
#include "msk/utility/convolution.hpp"


namespace MSk {

namespace block {


/**
 * @brief Data representation for dense block format.
 *
 * @tparam DataType Block element type.
 */
template<class DataType>
class Dense : public Data<DataType> {
private:
    using RealType = decltype(std::abs(DataType()));
    double _norm{0};
    SmartArray<DataType> _data;

private:
    /**
     * @brief Number of bytes required to store and recover the representation.
     *
     * @return uint64_t bytes.
     */
    uint64_t _local_packed_size() const {
        return sizeof(double) + + sizeof(uint64_t) + memory() * sizeof(DataType);
    }

public:
    Dense(uint64_t m, uint64_t n) : Data<DataType>(m, n) {}

    Dense(const Dense &) = default;

    Dense(Dense &&other) noexcept : Data<DataType>(std::move(other)) {
        _norm = other._norm;
        other._norm = 0.0;
        _data = std::move(other._data);
    }

    Dense &operator=(const Dense &) = default;

    Dense &operator=(Dense &&other) noexcept {
        Data<DataType>::operator=(std::move(other));
        _norm = other._norm;
        other._norm = 0.0;
        _data = std::move(other._data);
        return *this;
    }

    ~Dense() override = default;

    const DataType *raw_data() const {
        return _data.data();
    }

    /**
     * @brief Copy block elements from an array.
     *
     * @param[in] data SmartArray of elements.
     * @param[in] norm Precomputed norm of the block.
     *
     * @note If norm is not provided (equals zero) then it will be computed.
     */
    void set_data(const SmartArray<DataType> &data, double norm = 0.0) {
        _data = data;
        if (norm == 0.0) {
            _norm = BLAS::nrm2(this->n_rows * this->n_cols, _data.data(), 1);
        } else {
            _norm = norm;
        }
    }

    /**
     * @brief Move block elements from an array.
     *
     * @param[in, out] data SmartArray of elements.
     * @param[in] norm Precomputed norm of the block.
     *
     * @note If norm is not provided (equals zero) then it will be computed.
     */
    void set_data(SmartArray<DataType> &&data, double norm = 0.0) {
        _data = std::move(data);
        if (norm == 0.0) {
            _norm = BLAS::nrm2(this->n_rows * this->n_cols, _data.data(), 1);
        } else {
            _norm = norm;
        }
    }

    uint64_t memory() const override {
        return _data.size();
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
        BLAS::gemv('N', this->n_rows, this->n_cols, DataType{1.0}, _data.data(), this->n_rows,
                   x, 1, DataType{1.0}, y, 1);
    }

    void matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const override {
        if (_data.size() == 0) {
            return;
        }
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemm('N', 'N', m, block_size, n, DataType{1.0}, _data.data(), m, X, ldim_X, DataType{1.0}, Y, ldim_Y);
    }

    void compl_conj_matvec(const DataType *x, DataType *y) const override {
        if (_data.size() == 0) {
            return;
        }
        BLAS::gemv('Z', this->n_rows, this->n_cols, DataType{1.0}, _data.data(), this->n_rows,
                   x, 1, DataType{1.0}, y, 1);
    }

    void compl_conj_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const override {
        if (_data.size() == 0) {
            return;
        }
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemm('Z', 'N', m, block_size, n, DataType{1.0}, _data.data(), m, X, ldim_X, DataType{1.0}, Y, ldim_Y);
    }

    void transp_matvec(const DataType *x, DataType *y) const override {
        if (_data.size() == 0) {
            return;
        }
        BLAS::gemv('T', this->n_rows, this->n_cols, DataType{1.0}, _data.data(), this->n_rows,
                   x, 1, DataType{1.0}, y, 1);
    }

    void transp_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const override {
        if (_data.size() == 0) {
            return;
        }
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemm('T', 'N', n, block_size, m, DataType{1.0}, _data.data(), m, X, ldim_X, DataType{1.0}, Y, ldim_Y);
    }

    void conj_matvec(const DataType *x, DataType *y) const override {
        if (_data.size() == 0) {
            return;
        }
        BLAS::gemv('C', this->n_rows, this->n_cols, DataType{1.0}, _data.data(), this->n_rows,
                   x, 1, DataType{1.0}, y, 1);
    }

    void conj_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const override {
        if (_data.size() == 0) {
            return;
        }
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemm('C', 'N', n, block_size, m, DataType{1.0}, _data.data(), m, X, ldim_X, DataType{1.0}, Y, ldim_Y);
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
        convolve_add_dense(m, n, _data.data(), y, x, z, size);
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
        conj_convolve_add_dense(m, n, _data.data(), y, x, z, size);
    }

    void to_dense(DataType *ptr, uint64_t ldim) const override {
        if (ldim == this->n_rows) {
            BLAS::copy(this->n_rows * this->n_cols, _data.data(), 1, ptr, 1);
        } else {
            for (uint64_t j = 0; j < this->n_cols; ++j) {
                BLAS::copy(this->n_rows, _data.data() + j * this->n_rows, 1, ptr + j * ldim, 1);
            }
        }
    }

    uint64_t pack(uint8_t *ptr) const override {
        uint64_t data_off{Data<DataType>::pack(ptr)};
        if (ptr != nullptr) {
            ptr += data_off;
            reinterpret_cast<double *>(ptr)[0] = _norm;
            ptr += sizeof(double);
            uint64_t mem{memory()};
            reinterpret_cast<uint64_t *>(ptr)[0] = mem;
            ptr += sizeof(uint64_t);
            if (mem > 0) {
                std::memcpy(ptr, _data.data(), mem * sizeof(DataType));
            }
        }
        return data_off + _local_packed_size();
    }

    uint64_t unpack(uint8_t *ptr) override {
        uint64_t data_off{Data<DataType>::unpack(ptr)};
        if (ptr != nullptr) {
            ptr += data_off;
            _norm = reinterpret_cast<double *>(ptr)[0];
            ptr += sizeof(double);
            uint64_t mem = reinterpret_cast<uint64_t *>(ptr)[0];
            ptr += sizeof(uint64_t);
            if (mem > 0) {
                _data.resize(mem);
                std::memcpy(_data.data(), ptr, mem * sizeof(DataType));
            }
        }
        return data_off + _local_packed_size();
    }

    void free() override {
        Data<DataType>::free();
        _data.free();
        _norm = 0;
    }

    comp::Result recompress(comp::Parameters recomp_params) override {
        (void) recomp_params;
        return {_norm, 0, memory(), true};
    }
};


extern template class Dense<float>;

extern template class Dense<double>;

extern template class Dense<std::complex<float>>;

extern template class Dense<std::complex<double>>;


}
}


#endif //MOSAIC_SKELETON_BLOCK_DENSE_HPP
