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
// Created by bulat_v on 21.01.24.
//


/**
 * @file
 *
 * @brief Contains template classes for diagonal, scalar and zeroed matrix data format.
 */


#ifndef MOSAIC_SKELETON_BLOCK_DIAGONAL_HPP
#define MOSAIC_SKELETON_BLOCK_DIAGONAL_HPP


#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"

#include "msk/block/data.hpp"
#include "msk/utility/common.hpp"


namespace MSk {

namespace block {


/**
 * @brief Data representation for diagonal block format.
 *
 * @tparam DataType Block element type.
 */
template<class DataType>
class Diag : public Data<DataType> {
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
        return sizeof(double) + sizeof(uint64_t) + memory() * sizeof(DataType);
    }

public:
    Diag(uint64_t m, uint64_t n) : Data<DataType>(m, n) {}

    Diag(const Diag &) = default;

    Diag(Diag &&other) noexcept : Data<DataType>(std::move(other)) {
        _norm = other._norm;
        other._norm = 0.0;
        _data = std::move(other._data);
    }

    Diag &operator=(const Diag &) = default;

    Diag &operator=(Diag &&other) noexcept {
        Data<DataType>::operator=(std::move(other));
        _norm = other._norm;
        other._norm = 0.0;
        _data = std::move(other._data);
        return *this;
    }

    ~Diag() override = default;

    const DataType *raw_data() const {
        return _data.data();
    }

    /**
     * @brief Copy block elements from an array.
     *
     * @param[in] diagonal SmartArray of elements.
     * @param[in] norm Precomputed norm of the block.
     *
     * @note If norm is not provided (equals zero) then it will be computed.
     */
    void set_data(const SmartArray<DataType> &diagonal, double norm = 0.0) {
        uint64_t sz{std::min(this->n_rows, this->n_cols)};
        _data = diagonal;
        _data.resize(sz);
        if (norm == 0.0) {
            _norm = BLAS::nrm2(sz, _data.data(), 1);
        } else {
            _norm = norm;
        }
    }

    /**
     * @brief Move block elements from an array.
     *
     * @param[in, out] diagonal SmartArray of elements.
     * @param[in] norm Precomputed norm of the block.
     *
     * @note If norm is not provided (equals zero) then it will be computed.
     */
    void set_data(SmartArray<DataType> &&diagonal, double norm = 0.0) {
        uint64_t sz{std::min(this->n_rows, this->n_cols)};
        _data = std::move(diagonal);
        _data.resize(sz);
        if (norm == 0.0) {
            _norm = BLAS::nrm2(sz, _data.data(), 1);
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
        for (uint64_t i = 0; i < std::min(this->n_rows, this->n_cols); ++i) {
            y[i] += _data[i] * x[i];
        }
    }

    void transp_matvec(const DataType *x, DataType *y) const override {
        matvec(x, y);
    }

    void conj_matvec(const DataType *x, DataType *y) const override {
        if (_data.size() == 0) {
            return;
        }
        for (uint64_t i = 0; i < std::min(this->n_rows, this->n_cols); ++i) {
            y[i] += conjugate(_data[i]) * x[i];
        }
    }

    void compl_conj_matvec(const DataType *x, DataType *y) const override {
        conj_matvec(x, y);
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



/**
 * @brief Data representation for scalar block format.
 *
 * @tparam DataType Block element type.
 */
template<class DataType>
class Scalar : public Data<DataType> {
private:
    using RealType = decltype(std::abs(DataType()));
    DataType _element{};

public:
    Scalar(uint64_t m, uint64_t n, DataType elem = {}) : Data<DataType>(m, n), _element{elem} {}

    Scalar(const Scalar &) = default;

    Scalar(Scalar &&other) noexcept = default;

    Scalar &operator=(const Scalar &) = default;

    Scalar &operator=(Scalar &&other) noexcept = default;

    ~Scalar() override = default;

    /**
     * @brief Copy block elements from an array.
     *
     * @param[in] elem Diagonal element of the block.
     */
    void set_data(const DataType &elem) {
        _element = elem;
    }

    uint64_t memory() const override {
        return 1;
    }

    RealType frob_norm() const override {
        return std::sqrt(std::min(this->n_rows, this->n_cols)) * std::abs(_element);
    }

    bool exact() const override {
        return true;
    }

    void matvec(const DataType *x, DataType *y) const override {
        for (uint64_t i = 0; i < std::min(this->n_rows, this->n_cols); ++i) {
            y[i] += _element * x[i];
        }
    }

    void transp_matvec(const DataType *x, DataType *y) const override {
        matvec(x, y);
    }

    void conj_matvec(const DataType *x, DataType *y) const override {
        const DataType C{conjugate(_element)};
        for (uint64_t i = 0; i < std::min(this->n_rows, this->n_cols); ++i) {
            y[i] += C * x[i];
        }
    }

    void compl_conj_matvec(const DataType *x, DataType *y) const override {
        conj_matvec(x, y);
    }    

    void convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size) const override {
        const uint64_t m{this->n_rows}, n{this->n_cols};
        if (size == 0) {
            size = m + n - 1;
        }
        size = std::min(size, m + n - 1);
        for (uint64_t i = 0; i < std::min(this->n_rows, this->n_cols) && 2 * i < size; ++i) {
            z[2 * i] += y[i] * _element * x[i];
        }
    }

    void conj_convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size) const override {
        const uint64_t m{this->n_rows}, n{this->n_cols};
        if (size == 0) {
            size = m + n - 1;
        }
        size = std::min(size, m + n - 1);
        for (uint64_t i = 0; i < std::min(this->n_rows, this->n_cols) && 2 * i < size; ++i) {
            z[2 * i] += y[i] * conjugate(_element) * x[i];
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
        BLAS::copy(std::min(this->n_rows, this->n_cols), &_element, 0, ptr, ldim + 1);
    }

    uint64_t pack(uint8_t *ptr) const override {
        uint64_t data_off{Data<DataType>::pack(ptr)};
        if (ptr != nullptr) {
            ptr += data_off;
            std::memcpy(ptr, &_element, sizeof(DataType));
        }
        return data_off + sizeof(DataType);
    }

    uint64_t unpack(uint8_t *ptr) override {
        uint64_t data_off{Data<DataType>::unpack(ptr)};
        if (ptr != nullptr) {
            ptr += data_off;
            std::memcpy(&_element, ptr, sizeof(DataType));
        }
        return data_off + sizeof(DataType);
    }

    void free() override {
        Data<DataType>::free();
    }

    comp::Result recompress(comp::Parameters recomp_params) override {
        (void) recomp_params;
        return {frob_norm(), 0, 1, true};
    }
};


/**
 * @brief Data representation for zero block format.
 *
 * @tparam DataType Block element type.
 */
template<class DataType>
class Zero : public Scalar<DataType> {
public:
    Zero(uint64_t m, uint64_t n) : Scalar<DataType>(m, n, 0) {}

    void matvec(const DataType *x, DataType *y) const override {
        return;
    }

    void compl_conj_matvec(const DataType *x, DataType *y) const override {
        return;
    }

    void transp_matvec(const DataType *x, DataType *y) const override {
        return;
    }

    void conj_matvec(const DataType *x, DataType *y) const override {
        return;
    }

    void convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size) const override {
        return;
    }

    void conj_convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size) const override {
        return;
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
    }
};


extern template class Diag<float>;

extern template class Diag<double>;

extern template class Diag<std::complex<float>>;

extern template class Diag<std::complex<double>>;


extern template class Scalar<float>;

extern template class Scalar<double>;

extern template class Scalar<std::complex<float>>;

extern template class Scalar<std::complex<double>>;


extern template class Zero<float>;

extern template class Zero<double>;

extern template class Zero<std::complex<float>>;

extern template class Zero<std::complex<double>>;


}
}


#endif //MOSAIC_SKELETON_BLOCK_DIAGONAL_HPP
