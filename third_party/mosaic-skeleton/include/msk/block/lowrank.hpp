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
// Created by b.valiakhmetov on 17.02.23.
//


/**
 * @file
 *
 * @brief Contains template class block::LowRank, which represents low-rank block data format.
 */


#ifndef MOSAIC_SKELETON_BLOCK_LOWRANK_HPP
#define MOSAIC_SKELETON_BLOCK_LOWRANK_HPP


#include <cstdint>
#include <type_traits>

#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"
#include "cpp_lapack/LAPACK.hpp"

#include "msk/block/data.hpp"
#include "msk/block/dense.hpp"
#include "msk/utility/svd.hpp"
#include "msk/utility/convolution.hpp"


namespace MSk {

namespace block {


/**
 * @brief Data representation for low-rank block format.
 *
 * @details Block is stored as the product of two factors: A = U V^*.
 *
 * @tparam DataType Block element type.
 */
template<class DataType>
class LowRank : public Data<DataType> {
private:
    using RealType = decltype(std::abs(DataType()));
    double _norm{0};
    uint64_t _rank{0};

    bool _orth_cols{false}, _orth_rows{false};
    SmartArray<DataType> col_factor; /**< Columns' basis U */
    SmartArray<DataType> row_factor; /**< Rows' conjugated basis V */

private:
    /**
     * @brief Number of bytes required to store and recover the representation.
     *
     * @return uint64_t bytes.
     */
    uint64_t _local_packed_size() const {
        return sizeof(double) + sizeof(uint64_t) + memory() * sizeof(DataType);
    }

    // Transposed matvec for complex type
    template<typename DummyVoid = void>
    typename std::enable_if<!std::is_same<DataType, RealType>::value, DummyVoid>::type
    _transp_matvec(const DataType *x, DataType *y) const {
        uint64_t m{this->n_rows}, n{this->n_cols};
        SmartArray<DataType> tmp(_rank), x_conj(m), y_conj(n);

        BLAS::copy(m, x, 1, x_conj.data(), 1);
        LAPACK::lacgv(m, x_conj.data(), 1);

        BLAS::gemv('C', m, _rank, DataType{1.0}, col_factor.data(), m,
                   x_conj.data(), 1, DataType{0.0}, tmp.data(), 1);
        BLAS::gemv('N', n, _rank, DataType{1.0}, row_factor.data(), n,
                   tmp.data(), 1, DataType{0.0}, y_conj.data(), 1);
                   
        LAPACK::lacgv(n, y_conj.data(), 1);
        BLAS::axpy(n, DataType{1.0}, y_conj.data(), 1, y, 1);
    }

    // Transposed matvec for real type
    template<typename DummyVoid = void>
    typename std::enable_if<std::is_same<DataType, RealType>::value, DummyVoid>::type
    _transp_matvec(const DataType *x, DataType *y) const {
        conj_matvec(x, y);
    }

    // Block transposed matvec for complex type
    template<typename DummyVoid = void>
    typename std::enable_if<!std::is_same<DataType, RealType>::value, DummyVoid>::type
    _transp_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const {
        uint64_t m{this->n_rows}, n{this->n_cols};
        SmartArray<DataType> tmp(_rank * block_size), X_conj(m * block_size), Y_conj(n * block_size);

        if (ldim_X == m) {
            BLAS::copy(m * block_size, X, 1, X_conj.data(), 1);
        } else {
            for (uint64_t k = 0; k < block_size; ++k) {
                BLAS::copy(m, X + k * ldim_X, 1, X_conj.data() + k * m, 1);
            }
        }

        LAPACK::lacgv(m * block_size, X_conj.data(), 1);

        BLAS::gemm('C', 'N', _rank, block_size, m, DataType{1.0}, col_factor.data(), m, 
                   X_conj.data(), m, DataType{0.0}, tmp.data(), _rank);
        BLAS::gemm('N', 'N', n, block_size, _rank, DataType{1.0}, row_factor.data(), n, 
                   tmp.data(), _rank, DataType{0.0}, Y_conj.data(), n);

        LAPACK::lacgv(n * block_size, Y_conj.data(), 1);
        
        if (ldim_Y == n) {
            BLAS::axpy(n * block_size, DataType{1.0}, Y_conj.data(), 1, Y, 1);
        } else {
            for (uint64_t k = 0; k < block_size; ++k) {
                BLAS::axpy(n, DataType{1.0}, Y_conj.data() + k * n, 1, Y + k * ldim_Y, 1);
            }
        }
    }

    // Block transposed matvec for real type
    template<typename DummyVoid = void>
    typename std::enable_if<std::is_same<DataType, RealType>::value, DummyVoid>::type
    _transp_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const {
        conj_matvec(block_size, X, ldim_X, Y, ldim_Y);
    }

    void recalculate_norm() {
        uint64_t m{this->n_rows}, n{this->n_cols};
        SmartArray<DataType> uu(_rank * _rank), vv(_rank * _rank);
        BLAS::gemm('C', 'N', _rank, _rank, m, DataType{1.0},
                    col_factor.data(), m, col_factor.data(), m,
                    DataType{0.0}, uu.data(), _rank);
        BLAS::gemm('C', 'N', _rank, _rank, n, DataType{1.0},
                    row_factor.data(), n, row_factor.data(), n,
                    DataType{0.0}, vv.data(), _rank);
        DataType norm{0};
        for (uint64_t r = 0; r < _rank; ++r) {
            norm += BLAS::dotc(_rank, uu.data() + r, _rank, vv.data() + r * _rank, 1);
        }
        _norm = std::sqrt(std::abs(norm));
    }

public:
    LowRank(uint64_t m, uint64_t n) : Data<DataType>(m, n) {}

    LowRank(const LowRank &) = default;

    LowRank(LowRank &&) = default;

    LowRank &operator=(const LowRank &) = default;

    LowRank &operator=(LowRank &&) = default;

    ~LowRank() override = default;

    const DataType *raw_col_factor() const {
        return col_factor.data();
    }

    const DataType *raw_row_factor() const {
        return row_factor.data();
    }

    /**
     * @brief Copy block representation.
     *
     * @param[in] cols Basis of columns.
     * @param[in] rows Conjugated basis of rows.
     * @param[in] rank Number of columns/rows in bases.
     * @param[in] norm Precomputed norm of the block.
     *
     * @note If norm is not provided (equals zero) then it will be computed.
     */
    void set_data(const SmartArray<DataType> &cols, const SmartArray<DataType> &rows,
                  uint64_t rank, double norm = 0.0) {
        uint64_t m{this->n_rows}, n{this->n_cols};
        col_factor = cols;
        row_factor = rows;
        _rank = rank;
        if (norm == 0.0) {
            recalculate_norm();
        } else {
            _norm = norm;
        }
    }

    /**
     * @brief Move block representation.
     *
     * @param[in, out] cols Basis of columns.
     * @param[in, out] rows Conjugated basis of rows.
     * @param[in] rank Number of columns/rows in bases.
     * @param[in] norm Precomputed norm of the block.
     *
     * @note If norm is not provided (equals zero) then it will be computed.
     */
    void set_data(SmartArray<DataType> &&cols, SmartArray<DataType> &&rows,
                  uint64_t rank, double norm = 0.0) {
        uint64_t m{this->n_rows}, n{this->n_cols};
        col_factor = std::move(cols);
        row_factor = std::move(rows);
        _rank = rank;
        if (norm == 0.0) {
            recalculate_norm();
        } else {
            _norm = norm;
        }
    }

    uint64_t memory() const override {
        return col_factor.size() + row_factor.size();
    }

    RealType frob_norm() const override {
        return _norm;
    }

    bool exact() const override {
        return false;
    }

    void matvec(const DataType *x, DataType *y) const override {
        SmartArray<DataType> tmp(_rank);
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemv('C', n, _rank, DataType{1.0}, row_factor.data(), n,
                   x, 1, DataType{0.0}, tmp.data(), 1);
        BLAS::gemv('N', m, _rank, DataType{1.0}, col_factor.data(), m,
                   tmp.data(), 1, DataType{1.0}, y, 1);
    }

    void matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const override {
        SmartArray<DataType> tmp(_rank * block_size);
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemm('C', 'N', _rank, block_size, n, DataType{1.0}, row_factor.data(), n, 
                   X, ldim_X, DataType{0.0}, tmp.data(), _rank);
        BLAS::gemm('N', 'N', m, block_size, _rank, DataType{1.0}, col_factor.data(), m, 
                   tmp.data(), _rank, DataType{1.0}, Y, ldim_Y);
    }

    void transp_matvec(const DataType *x, DataType *y) const override {
        _transp_matvec(x, y);
    }

    void transp_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const override {
        _transp_matvec(block_size, X, ldim_X, Y, ldim_Y);
    }

    void conj_matvec(const DataType *x, DataType *y) const override {
        SmartArray<DataType> tmp(_rank);
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemv('C', m, _rank, DataType{1.0}, col_factor.data(), m,
                   x, 1, DataType{0.0}, tmp.data(), 1);
        BLAS::gemv('N', n, _rank, DataType{1.0}, row_factor.data(), n,
                   tmp.data(), 1, DataType{1.0}, y, 1);
    }

    void conj_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const override {
        SmartArray<DataType> tmp(_rank * block_size);
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemm('C', 'N', _rank, block_size, m, DataType{1.0}, col_factor.data(), m, 
                   X, ldim_X, DataType{0.0}, tmp.data(), _rank);
        BLAS::gemm('N', 'N', n, block_size, _rank, DataType{1.0}, row_factor.data(), n, 
                   tmp.data(), _rank, DataType{1.0}, Y, ldim_Y);
    }

    void compl_conj_matvec(const DataType *x, DataType *y) const override {
        SmartArray<DataType> tmp(_rank);
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemv('T', n, _rank, DataType{1.0}, row_factor.data(), n,
                   x, 1, DataType{0.0}, tmp.data(), 1);
        BLAS::gemv('Z', m, _rank, DataType{1.0}, col_factor.data(), m,
                   tmp.data(), 1, DataType{1.0}, y, 1);
    }

    void compl_conj_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const override {
        SmartArray<DataType> tmp(_rank * block_size);
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemm('T', 'N', _rank, block_size, n, DataType{1.0}, row_factor.data(), n, 
                   X, ldim_X, DataType{0.0}, tmp.data(), _rank);
        BLAS::gemm('Z', 'N', m, block_size, _rank, DataType{1.0}, col_factor.data(), m, 
                   tmp.data(), _rank, DataType{1.0}, Y, ldim_Y);
    }

    void convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size) const override {
        const uint64_t m{this->n_rows}, n{this->n_cols};
        if (size == 0) {
            size = m + n - 1;
        }
        size = std::min(size, m + n - 1);

        // call fft-based convolution
        convolve_add_lowrank(m, n, _rank, col_factor.data(), row_factor.data(), y, x, z, size);
    }

    void conj_convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size) const override {
        const uint64_t m{this->n_rows}, n{this->n_cols};
        if (size == 0) {
            size = m + n - 1;
        }
        size = std::min(size, m + n - 1);

        // call fft-based convolution
        convolve_add_lowrank(n, m, _rank, row_factor.data(), col_factor.data(), y, x, z, size);
    }

    void to_dense(DataType *ptr, uint64_t ldim) const override {
        uint64_t m{this->n_rows}, n{this->n_cols};
        BLAS::gemm('N', 'C', m, n, _rank, DataType{1.0},
                   col_factor.data(), m, row_factor.data(), n,
                   DataType{0.0}, ptr, ldim);
    }

    /**
     * @brief Get the block rank.
     *
     * @return uint64_t rank.
     */
    uint64_t rank() const {
        return _rank;
    }

    /**
     * @brief Check if rows basis is orthogonal.
     *
     * @return Orthogonality flag.
     */
    bool orth_rows() const {
        return _orth_rows;
    }

    /**
     * @brief Check if columns basis is orthogonal.
     *
     * @return Orthogonality flag.
     */
    uint64_t orth_cols() const {
        return _orth_cols;
    }

    uint64_t pack(uint8_t *ptr) const override {
        uint64_t data_off{Data<DataType>::pack(ptr)};
        if (ptr != nullptr) {
            ptr += data_off;
            reinterpret_cast<double *>(ptr)[0] = _norm;
            ptr += sizeof(double);
            reinterpret_cast<uint64_t *>(ptr)[0] = _rank;
            ptr += sizeof(uint64_t);
            if (_rank > 0) {
                uint64_t col_factor_size{col_factor.size() * sizeof(DataType)};
                std::memcpy(ptr, col_factor.data(), col_factor_size);
                ptr += col_factor_size;
                uint64_t row_factor_size{row_factor.size() * sizeof(DataType)};
                std::memcpy(ptr, row_factor.data(), row_factor_size);
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
            _rank = reinterpret_cast<uint64_t *>(ptr)[0];
            ptr += sizeof(uint64_t);
            if (_rank > 0) {
                col_factor.resize(this->n_rows * _rank);
                uint64_t col_factor_size{this->n_rows * _rank * sizeof(DataType)};
                std::memcpy(col_factor.data(), ptr, col_factor_size);
                ptr += col_factor_size;
                row_factor.resize(this->n_cols * _rank);
                uint64_t row_factor_size{this->n_cols * _rank * sizeof(DataType)};
                std::memcpy(row_factor.data(), ptr, row_factor_size);
            }
        }
        return data_off + _local_packed_size();
    }

    void free() override {
        Data<DataType>::free();
        col_factor.free();
        row_factor.free();
        _norm = 0;
        _rank = 0;
    }

    comp::Result recompress(comp::Parameters recomp_params) override {
        uint64_t m{this->n_rows}, n{this->n_cols};
        if (m == 0 || n == 0 || col_factor.size() == 0 || row_factor.size() == 0) {
            return {_norm, 0, _rank * (m + n), true};
        }
        auto res = recompress_lowrank(m, n, _rank, col_factor, row_factor, recomp_params, _orth_cols, _orth_rows);
        _rank = res.memory / (m + n);
        _norm = res.matr_norm;
        _orth_cols = true;
        _orth_rows = true;
        col_factor.resize(m * _rank);
        row_factor.resize(n * _rank);
        return res;
    }

    /**
     * @brief Retreive columns basis U to an array.
     *
     * @param[in] ptr Pointer to write to.
     * @param[in] ldim Leading dimension of a resulting array.
     */
    void extract_col_factor(DataType *ptr, uint64_t ldim) const {
        uint64_t m{this->n_rows};
        if (ldim == m) {
            BLAS::copy(m * _rank, col_factor.data(), 1, ptr, 1);
        } else {
            for (uint64_t r = 0; r < _rank; ++r) {
                BLAS::copy(m, col_factor.data() + r * m, 1, ptr + r * ldim, 1);
            }
        }
    }

    /**
     * @brief Retreive conjugated rows basis V to an array.
     *
     * @param[in] ptr Pointer to write to.
     * @param[in] ldim Leading dimension of a resulting array.
     */
    void extract_row_factor(DataType *ptr, uint64_t ldim) const {
        uint64_t n{this->n_cols};
        if (ldim == n) {
            BLAS::copy(n * _rank, row_factor.data(), 1, ptr, 1);
        } else {
            for (uint64_t r = 0; r < _rank; ++r) {
                BLAS::copy(n, row_factor.data() + r * n, 1, ptr + r * ldim, 1);
            }
        }
    }

    explicit operator Dense<DataType>() const {
        Dense<DataType> dns{this->n_rows, this->n_cols};
        SmartArray<DataType> data{this->n_rows * this->n_cols};
        to_dense(data.data(), this->n_rows);
        dns.set_data(std::move(data), _norm);
        return dns;
    }
};


extern template class LowRank<float>;

extern template class LowRank<double>;

extern template class LowRank<std::complex<float>>;

extern template class LowRank<std::complex<double>>;


}
}


#endif //MOSAIC_SKELETON_BLOCK_LOWRANK_HPP
