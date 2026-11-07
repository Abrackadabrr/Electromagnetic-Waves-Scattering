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
 * @brief Contains template abstract class block::Data, which represents an arbitrary structured matrix.
 */


#ifndef MOSAIC_SKELETON_BLOCK_DATA_HPP
#define MOSAIC_SKELETON_BLOCK_DATA_HPP


#include <cmath>
#include <complex>
#include <cstdint>
#include <utility>

#include "msk/computation/computer.hpp"


namespace MSk {

namespace block {


/**
 * @brief Abstract rectangular block representation.
 *
 * Interface to a block data representation for an arbitrary matrix format.
 * It requires some statistics, matrix-vector product, convertion to dense form
 * and packing/unpacking to bytes sequence for transfering reasons.
 *
 * @tparam DataType Block element type.
 */
template<class DataType>
class Data {
    /**
     * @brief Number of bytes used to store the representation.
     *
     * @return uint64_t bytes.
     */
    uint64_t _local_packed_size() const {
        return 2 * sizeof(uint64_t);
    }

protected:
    using RealType = decltype(std::abs(DataType()));
    uint64_t n_rows{0}, n_cols{0};

public:
    /**
     * @brief Basic block Data requires only size.
     *
     * @param[in] m Number of rows.
     * @param[in] n Number of columns.
     */
    Data(uint64_t m, uint64_t n) : n_rows{m}, n_cols{n} {}

    virtual ~Data() = default;

    /**
     * @brief Memory amount used to store the representation.
     *
     * @return Number of DataType memory cells.
     */
    virtual uint64_t memory() const = 0;

    /**
     * @brief Frobenius norm of the block.
     *
     * @return RealType.
     */
    virtual RealType frob_norm() const = 0;

    /**
     * @brief Additive matrix-vector product: y += A x.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     */
    virtual void matvec(const DataType *x, DataType *y) const = 0;

    /**
     * @brief Additive matrix-block-vector product: Y += A X.
     *
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     */
    virtual void matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const {
        for (uint64_t k = 0; k < block_size; ++k) {
            matvec(X + k * ldim_X, Y + k * ldim_Y);
        }
    }

    /**
     * @brief Additive matrix-vector product: y += A x.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     */
    virtual void compl_conj_matvec(const DataType *x, DataType *y) const = 0;

    /**
     * @brief Additive matrix-block-vector product: Y += A X.
     *
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     */
    virtual void compl_conj_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const {
        for (uint64_t k = 0; k < block_size; ++k) {
            compl_conj_matvec(X + k * ldim_X, Y + k * ldim_Y);
        }
    }

    /**
     * @brief Additive transposed matrix-vector product: y += A^* x.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     */
    virtual void transp_matvec(const DataType *x, DataType *y) const = 0;

    /**
     * @brief Additive transposed matrix-block-vector product: Y += A^* X.
     *
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     */
    virtual void transp_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const {
        for (uint64_t k = 0; k < block_size; ++k) {
            transp_matvec(X + k * ldim_X, Y + k * ldim_Y);
        }
    }

    /**
     * @brief Additive conjugated matrix-vector product: y += A^* x.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     */
    virtual void conj_matvec(const DataType *x, DataType *y) const = 0;

    /**
     * @brief Additive block conjugated matrix-block-vector product: Y += A^* X.
     *
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     */
    virtual void conj_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const {
        for (uint64_t k = 0; k < block_size; ++k) {
            conj_matvec(X + k * ldim_X, Y + k * ldim_Y);
        }
    }

    /**
     * @brief Convolve two vectors using the block as the kernel and add to a vector.
     *
     * @details Convolution result contains first size entries of 'full' convolution.
     * So z_k += SUM[over j] y_j Data_{j, k - j} x_{k - j} is computed for k < size
     * (indexing starts from zero). Parameter size=0 is also treated as 'full'
     * of size M + N - 1 for Data of size (M, N).
     *
     * @param[in] y Left-hand-side vector to convolve.
     * @param[in] x Right-hand-side vector to convolve.
     * @param[in, out] z Resulting convolution vector.
     * @param[in] size Required size of result.
     */
    virtual void convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size) const = 0;

    /**
     * @brief Convolve two vectors using the block as the conjugated kernel and add to a vector.
     *
     * @details Convolution result contains first size entries of 'full' convolution.
     * So z_k += SUM[over j] y_j conjugate(Data_{k - j, j}) x_{k - j} is computed for k < size
     * (indexing starts from zero). Parameter size=0 is also treated as 'full'
     * of size M + N - 1 for Data of size (M, N).
     *
     * @param[in] y Left-hand-side vector to convolve.
     * @param[in] x Right-hand-side vector to convolve.
     * @param[in, out] z Resulting convolution vector.
     * @param[in] size Required size of result.
     */
    virtual void conj_convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size) const = 0;

    /**
     * @brief Store block in dense form.
     *
     * @param[out] ptr Pointer to the block write the result to.
     * @param[in] ldim Leading dimension of resulting block.
     */
    virtual void to_dense(DataType *ptr, uint64_t ldim) const = 0;

    /**
     * @brief Flag that shows if data representation stores exact or approximated block.
     *
     * @return Boolean flag.
     */
    virtual bool exact() const = 0;

    /**
     * @brief Size of the block.
     *
     * @return Pair M, N - number of rows and columns respectively.
     */
    std::pair<uint64_t, uint64_t> size() const {
        return {n_rows, n_cols};
    }

    /**
     * @brief Pack the representation as a sequence of bytes, or get the size is requires.
     *
     * @param ptr[out] Pointer to array to write the representation, or nullptr if the size is requested.
     * @return Number of bytes required to store the representation.
     */
    virtual uint64_t pack(uint8_t *ptr) const {
        if (ptr != nullptr) {
            uint64_t *put_ptr = reinterpret_cast<uint64_t *>(ptr);
            put_ptr[0] = n_rows;
            put_ptr[1] = n_cols;
        }
        return _local_packed_size();
    }

    /**
     * @brief Unpack the representation from a sequence of bytes, or get the size is requires.
     *
     * @param ptr[in] Pointer to array to read the representation, or nullptr if the size is requested.
     * @return Number of bytes required to store the representation.
     */
    virtual uint64_t unpack(uint8_t *ptr) {
        if (ptr != nullptr) {
            uint64_t *get_ptr = reinterpret_cast<uint64_t *>(ptr);
            n_rows = get_ptr[0];
            n_cols = get_ptr[1];
        }
        return _local_packed_size();
    }

    /**
     * @brief Free the representation: destroy the content, deallocate memory.
     */
    virtual void free() {}

    /**
     * @brief Recompress data representation if possible.
     *
     * @param[in] recomp_params Recompression parameters.
     * @return Recompression statistics.
     */
    virtual comp::Result recompress(comp::Parameters recomp_params) = 0;
};


extern template class Data<float>;

extern template class Data<double>;

extern template class Data<std::complex<float>>;

extern template class Data<std::complex<double>>;


}

}


#endif //MOSAIC_SKELETON_BLOCK_DATA_HPP
