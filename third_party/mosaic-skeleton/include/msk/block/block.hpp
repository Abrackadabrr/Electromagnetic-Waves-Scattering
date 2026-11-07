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
// created by bulat_v on 03.02.2022.
//


/**
 * @file
 *
 * @brief Contains template class block::Block, which is used to operate on a block of the matrix.
 */


#ifndef MOSAIC_SKELETON_BLOCK_BLOCK_HPP
#define MOSAIC_SKELETON_BLOCK_BLOCK_HPP


#include "msk_config.h"

#include <cassert>
#include <cstdint>
#include <utility>

#ifdef MSK_USE_MPI
#include "cpp_mpi/MPI.hpp"
#endif

#include "msk/index/index.hpp"
#include "msk/block/data.hpp"
#include "msk/block/dense.hpp"
#include "msk/matrix/matrix.hpp"
#include "msk/evaluation/request.hpp"
#include "msk/computation/computer.hpp"


namespace MSk {

namespace block {


/**
 * @brief Subblock of a matrix.
 *
 * Block provides indexing over rows and columns.
 * It is utilized as the interface to an abstract data representation format
 * and an abstract method of its computation.
 *
 * @tparam DataType Matrix element type.
 * @tparam Idx Indexing interface.
 */
template<class DataType, class Idx = index::Index<uint64_t>>
class Block {
protected:
    uint64_t _id{0};
    Idx row_idx, col_idx;
    using RealType = decltype(std::abs(DataType()));

private:
    comp::Computer<DataType> *_computer{nullptr};
    Data<DataType> *_data{nullptr};

#ifdef MSK_USE_MPI
    MPI_Comm _comm{MPI_COMM_NULL};
#endif

public:
    Block() : row_idx{}, col_idx{} {}

    /**
     * @brief Construct Block over a given pair of row and column indices with optional data and its computer.
     *
     * @param[in] row Row indices.
     * @param[in] col Column indices.
     * @param[in] data Data representation.
     * @param[in] computer Computer that can evaluate the data given.
     */
    Block(const Idx &row, const Idx &col,
          Data<DataType> *data = nullptr, comp::Computer<DataType> *computer = nullptr) :
          row_idx{row}, col_idx{col}, _computer{computer}, _data{data} {}

    Block(const Block &other) = delete;

    Block(Block &&other) noexcept {
        _id = std::move(other._id);
        row_idx = std::move(other.row_idx);
        col_idx = std::move(other.col_idx);
        _computer = std::move(other._computer);
        other._computer = nullptr;
        _data = std::move(other._data);
        other._data = nullptr;
#ifdef MSK_USE_MPI
        _comm = other._comm;
        other._comm = MPI_COMM_NULL;
#endif
    }

    Block &operator=(const Block &other) = delete;

    Block &operator=(Block &&other) noexcept {
        if (this == &other) {
            return *this;
        }
        _id = std::move(other._id);
        row_idx = std::move(other.row_idx);
        col_idx = std::move(other.col_idx);
        _computer = std::move(other._computer);
        other._computer = nullptr;
        _data = std::move(other._data);
        other._data = nullptr;
#ifdef MSK_USE_MPI
        _comm = other._comm;
        other._comm = MPI_COMM_NULL;
#endif
        return *this;
    }

    /**
     * @brief Destruct internal data and computer.
     */
    ~Block() {
        delete _data;
        delete _computer;
    }

    /**
     * @brief Set data representation.
     *
     * @param[in] ptr Pointer to allocated data.
     */
    void set_data(Data<DataType> *ptr) {
        if (_data != nullptr) {
            delete _data;
        }
        _data = ptr;
    }

    /**
     * @brief Get the data representation.
     */
    Data<DataType> *data() const {
        return _data;
    }

    /**
     * @brief Set the data evaluation method.
     *
     * @param[in] ptr Pointer to allocated computer.
     */
    void set_computer(comp::Computer<DataType> *ptr) {
        if (_computer != nullptr) {
            delete _computer;
        }
        _computer = ptr;
    }

    /**
     * @brief Get the data computer.
     */
    comp::Computer<DataType> *computer() const {
        return _computer;
    }

    void set_id(uint64_t new_id) {
        _id = new_id;
    }

    /**
     * @brief Get user defined ID if a Block.
     *
     * @return id
     */
    uint64_t id() const {
        return _id;
    }

#ifdef MSK_USE_MPI
    /**
     * @brief Get MPI communicator.
     *
     * @return MPI_Comm
     */
    MPI_Comm comm() const {
        return _comm;
    }

    void set_comm(MPI_Comm comm) {
        _comm = comm;
    }
#endif

    /**
     * @brief Perform matrix-vector product in form: y += A x.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     */
    void matvec(const DataType *x, DataType *y) const {
        if (_data != nullptr) {
            _data->matvec(x, y);
        }
    }

    /**
     * @brief Perform matrix-block-vector product in form: Y += A X.
     *
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     */
    void matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const {
        if (_data != nullptr) {
            _data->matvec(block_size, X, ldim_X, Y, ldim_Y);
        }
    }

    /**
     * @brief Perform matrix-vector product in form: y += A x.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     */
    void compl_conj_matvec(const DataType *x, DataType *y) const {
        if (_data != nullptr) {
            _data->compl_conj_matvec(x, y);
        }
    }

    /**
     * @brief Perform matrix-block-vector product in form: Y += A X.
     *
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     */
    void compl_conj_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const {
        if (_data != nullptr) {
            _data->compl_conj_matvec(block_size, X, ldim_X, Y, ldim_Y);
        }
    }

    /**
     * @brief Perform transposed matrix-vector product in form: y += A^T x.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     */
    void transp_matvec(const DataType *x, DataType *y) const {
        if (_data != nullptr) {
            _data->transp_matvec(x, y);
        }
    }

    /**
     * @brief Perform transposed matrix-block-vector product in form: Y += A^T X.
     *
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     */
    void transp_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const {
        if (_data != nullptr) {
            _data->transp_matvec(block_size, X, ldim_X, Y, ldim_Y);
        }
    }

    /**
     * @brief Perform conjugated matrix-vector product in form: y += A^* x.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     */
     void conj_matvec(const DataType *x, DataType *y) const {
        if (_data != nullptr) {
            _data->conj_matvec(x, y);
        }
    }
    
    /**
     * @brief Perform conjugated matrix-block-vector product in form: Y += A^* X.
     *
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     */
    void conj_matvec(const uint64_t block_size, const DataType *X, const uint64_t ldim_X, DataType *Y, const uint64_t ldim_Y) const {
        if (_data != nullptr) {
            _data->conj_matvec(block_size, X, ldim_X, Y, ldim_Y);
        }
    }

    /**
     * @brief Convolve two vectors using the block as the kernel.
     *
     * @details For details see Data::convolve.
     *
     * @param[in] y Left-hand-side vector to convolve.
     * @param[in] x Right-hand-side vector to convolve.
     * @param[out] z Resulting convolution vector.
     * @param[in] size Required size of result.
     */
    void convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size = 0) const {
        const uint64_t m = row_idx.size(), n = col_idx.size();
        size = std::min(size, m + n - 1);
        if (size == 0) {
            size = m + n - 1;
        }
        if (_data != nullptr) {
            _data->convolve(y, x, z, size);
        }
    }

    /**
     * @brief Convolve two vectors using the block as the kernel.
     *
     * @details For details see Data::convolve.
     *
     * @param[in] y Left-hand-side vector to convolve.
     * @param[in] x Right-hand-side vector to convolve.
     * @param[out] z Resulting convolution vector.
     * @param[in] size Required size of result.
     */
    void conj_convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size = 0) const {
        const uint64_t m = row_idx.size(), n = col_idx.size();
        size = std::min(size, m + n - 1);
        if (size == 0) {
            size = m + n - 1;
        }
        if (_data != nullptr) {
            _data->conj_convolve(y, x, z, size);
        }
    }

    /**
     * @brief Store a block in dense form.
     *
     * @param[out] ptr Pointer to place data.
     * @param[in] ldim Leading dimension.
     */
    void to_dense(DataType *ptr, uint64_t ldim) const {
        if (_data != nullptr) {
            _data->to_dense(ptr, ldim);
        }
    }

    /**
     * @brief Frobenius norm of a block.
     *
     * @return Frobenius norm (real value).
     */
    RealType frob_norm() const {
        if (_data == nullptr) {
            return 0;
        }
        return _data->frob_norm();
    }

    /**
     * @brief Number of DataType values used to store data representation of a Block.
     *
     * @return Integer number of memory cells.
     */
    uint64_t memory() const {
        if (_data == nullptr) {
            return 0;
        }
        return _data->memory();
    }

    /**
     * @brief Flag that shows if block's data representation is exact or approximated.
     *
     * @return Boolean flag.
     */
    bool exact() const {
        if (_data == nullptr) {
            return true;
        }
        return _data->exact();
    }

    /**
     * @brief Compute defined data representation of a Block via computer.
     *
     * @param[in] A Matrix from which to evaluate elements.
     * @return Computation result.
     */
    comp::Result compute(matrix::Matrix<DataType> &A) {
        if (_data == nullptr || _computer == nullptr) {
            return {0, 0, 0};
        }
        return _computer->compute(A, row_idx, col_idx);
    }

    /**
     * @brief Perform one step of comutation and request additional elements of a block.
     *
     * @param[out] request Request to put elements to be evaluated.
     * @return Intermediate computation result.
     */
    comp::Result step_enqueue(eval::Request<DataType> &request) {
        if (_data == nullptr || _computer == nullptr) {
            return {0, 0, 0};
        }
        comp::Result result = _computer->step_enqueue(request);
        request.at_index(row_idx, col_idx);
        request.set_block_id(_id);
        return result;
    }

#ifdef MSK_USE_MPI
    /**
     * @brief Move block to another process (MPI rank).
     *
     * @param[in] to Flag: to send (true) or receive (false) data.
     * @param[in] partner_id Process MPI rank to migrate to.
     * @param[in] migr_comm MPI communicator to work on.
     *
     * @note Operation is blocking.
     */
    void migrate(bool to, int partner_id, MPI_Comm migr_comm = MPI_COMM_NULL) {
        if (_data == nullptr) {
            throw std::runtime_error("Sending/receiving empty block");
        }
        // Use block's own communicator in case provided one is undefined (default)
        if (migr_comm == MPI_COMM_NULL) {
            migr_comm = _comm;
        }
        SmartArray<uint8_t> packed_data;
        if (to) { // send
            uint64_t packed_size{_data->pack(nullptr)};

            packed_data.resize(packed_size);
            _data->pack(packed_data.data());
            MPI_Send(packed_data.data(), packed_size, MPI_UINT8_T, partner_id, _id, migr_comm);
            _data->free();
        } else { // receive
            MPI_Status stat;
            MPI_Probe(partner_id, _id, migr_comm, &stat);
            int packed_size{};
            MPI_Get_count(&stat, MPI_UINT8_T, &packed_size);

            packed_data.resize(packed_size);
            MPI_Recv(packed_data.data(), packed_size, MPI_UINT8_T, partner_id, _id, migr_comm, &stat);
            _data->unpack(packed_data.data());
        }
    }
#endif

    /**
     * @brief Perform recompression procedure.
     *
     * @param[in] params Recompression parameters.
     * @return Recomputation result.
     */
    comp::Result recompress(comp::Parameters params) {
        if (_data == nullptr) {
            return {0, 0, 0, true};
        }
        return _data->recompress(params);
    }

    void densify() {
        if (_data == nullptr || dynamic_cast<Dense<DataType> *>(_data) != nullptr) {
            return;
        }
        uint64_t m{row_idx.size()}, n{col_idx.size()};
        SmartArray<DataType> elements(m * n);
        _data->to_dense(elements.data(), m);
        auto new_data = new Dense<DataType>(m, n);
        new_data->set_data(std::move(elements), _data->frob_norm());
        delete _data;
        _data = new_data;
    }

    /**
     * @brief Size of a block: M x N.
     *
     * @return Pair M, N.
     */
    std::pair<typename Idx::idx_type, typename Idx::idx_type> size() const {
        return std::make_pair(row_idx.size(), col_idx.size());
    }

    /**
     * @brief Corner of a block inside the whole matrix.
     *
     * @return Pair i0, j0.
     */
    std::pair<typename Idx::idx_type, typename Idx::idx_type> corner() const {
        return std::make_pair(*row_idx.begin(false), *col_idx.begin(false));
    }

    /**
     * @brief Beginning and ending iterators over block's rows.
     * .
     * @param[in] use_perm Flag of using permutation.
     * @return Pair if iterators.
     */
    std::pair<typename Idx::iterator_type, typename Idx::iterator_type> row_range(bool use_perm = false) const {
        return std::make_pair(row_idx.begin(use_perm), row_idx.end(use_perm));
    }

    /**
     * @brief Beginning and ending iterators over block's columns.
     * .
     * @param[in] use_perm Flag of using permutation.
     * @return Pair if iterators.
     */
    std::pair<typename Idx::iterator_type, typename Idx::iterator_type> col_range(bool use_perm = false) const {
        return std::make_pair(col_idx.begin(use_perm), col_idx.end(use_perm));
    }

    friend bool operator==(const Block &lhs, const Block &rhs) {
        return lhs.row_idx == rhs.row_idx && lhs.col_idx == rhs.col_idx;
    }

    /**
     * @brief Row-major comparison of blocks positions in the matrix.
     *
     * @param[in] other Block to compare to.
     * @return True (integer > 0) if current block is positioned before another one.
     */
    int is_row_less(const Block &other) const {
        uint64_t my_i0, my_j0, ot_i0, ot_j0;
        std::tie(my_i0, my_j0) = corner();
        std::tie(ot_i0, ot_j0) = other.corner();
        if (my_i0 > ot_i0) {
            return 0;
        }
        if (my_i0 < ot_i0) {
            return 2;
        }
        if (my_j0 > ot_j0) {
            return 0;
        }
        if (my_j0 < ot_j0) {
            return 3;
        }
        return 0;
    }

    /**
     * @brief Columns-major comparison of blocks positions in the matrix.
     *
     * @param[in] other Block to compare to.
     * @return True (integer > 0) if current block is positioned before another one.
     */
    int is_col_less(const Block &other) const {
        uint64_t my_i0, my_j0, ot_i0, ot_j0;
        std::tie(my_i0, my_j0) = corner();
        std::tie(ot_i0, ot_j0) = other.corner();
        if (my_j0 > ot_j0) {
            return 0;
        }
        if (my_j0 < ot_j0) {
            return 2;
        }
        if (my_i0 > ot_i0) {
            return 0;
        }
        if (my_i0 < ot_i0) {
            return 3;
        }
        return 0;
    }

    /**
     * @brief Index over rows of a block.
     *
     * @return Index.
     */
    const Idx &rows() const {
        return row_idx;
    }

    /**
     * @brief Index over rows of a block.
     *
     * @return Index.
     */
    Idx &rows() {
        return row_idx;
    }

    /**
     * @brief Index over columns of a block.
     *
     * @return Index.
     */
    const Idx &columns() const {
        return col_idx;
    }

    /**
     * @brief Index over columns of a block.
     *
     * @return Index.
     */
    Idx &columns() {
        return col_idx;
    }

    /**
     * @brief Index over rows and columns of a block.
     *
     * @return Pair Index, Index.
     */
    std::pair<Idx, Idx> get_idx() const {
        return std::make_pair(row_idx, col_idx);
    }
};


extern template class Block<float>;

extern template class Block<double>;

extern template class Block<std::complex<float>>;

extern template class Block<std::complex<double>>;


}

}


#endif //MOSAIC_SKELETON_BLOCK_BLOCK_HPP
