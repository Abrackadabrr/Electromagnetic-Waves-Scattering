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
// Created by b.valiakhmetov on 21.01.24.
//


/**
 * @file
 *
 * @brief Defines dense block computation procedure.
 */


#ifndef MOSAIC_SKELETON_COMPUTATION_DAIGONAL_HPP
#define MOSAIC_SKELETON_COMPUTATION_DAIGONAL_HPP


#include <cmath>

#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"

#include "msk/computation/computer.hpp"
#include "msk/block/diagonal.hpp"
#include "msk/utility/logging.hpp"


namespace MSk {

namespace comp {


/**
 * @brief Basic procedure of a block evaluation in diagonal format.
 *
 * @tparam DataType Block element type.
 */
template<typename DataType>
class Diag: public Computer<DataType> {
private:
    using RealType = decltype(std::abs(DataType()));

    block::Diag<DataType> *_block{nullptr}; ///< Diagonal block to approximate

    Result approx_result{0.0, 0.0, 0, false}; ///< Current epproximation state

    SmartArray<DataType> _data{0};

    std::shared_ptr<log::Logger> log;

public:
    /**
     * @brief Initialize data to approximate.
     *
     * @param[in] data Diagonal block to compute.
     * @param[in] log_ptr Logger for computation procedure.
     */
    Diag(block::Diag<DataType> *data = nullptr,
         std::shared_ptr<log::Logger> log_ptr = nullptr) : _block{data}, log{std::move(log_ptr)} {
    }

    ~Diag() override = default;

    Result compute(matrix::Matrix<DataType> &A,
                   const index::Index<uint64_t> &row_idx,
                   const index::Index<uint64_t> &col_idx) override {
        uint64_t m{row_idx.size()}, n{col_idx.size()};
        if (_block == nullptr) {
            _block = new block::Diag<DataType>(m, n);
        } else if (std::make_pair(m, n) != _block->size()) {
            return approx_result = {0, 0, 0, false};
        }
        uint64_t sz{std::min(m, n)};
        _data.resize(sz);
        if (m == n) {
            A.list(row_idx, col_idx, _data.data());
        } else {
            std::vector<uint64_t> true_row_idx(sz), true_col_idx(sz);
            for (uint64_t i = 0; i < sz; ++i) {
                true_row_idx[i] = row_idx[i];
                true_col_idx[i] = col_idx[i];
            }
            A.list(true_row_idx, true_col_idx, _data.data());
        }
        RealType f_norm = BLAS::nrm2(sz, _data.data(), 1);
        _block->set_data(std::move(_data), f_norm);
        return approx_result = {f_norm, 0.0, sz, true};
    }

    Result step_enqueue(eval::Request<DataType> &request) override {
        if (approx_result.finished) { // If already approximated do nothing
            request.free();
            return approx_result;
        }

        uint64_t m{}, n{};
        std::tie(m, n) = _block->size();
        uint64_t sz{std::min(m, n)};
        if (approx_result.memory == 0) {
            _data.resize(sz);

            // Request all elements
            std::vector<uint64_t> row_idx(sz), col_idx(sz);
            for (uint64_t i = 0; i < sz; ++i) {
                row_idx[i] = col_idx[i] = i;
            }
            request.list(row_idx, col_idx, _data.data());
            return approx_result = {0, 0, sz, false};
        }

        // All elements are received
        approx_result.matr_norm = BLAS::nrm2(sz, _data.data(), 1);
        _block->set_data(std::move(_data), approx_result.matr_norm);
        return approx_result;
    }

    /**
     * @brief Get current computation result.
     *
     * @return Intermediate Result.
     */
    Result result() const {
        return approx_result;
    }

    /**
     * @brief Check if computation has finished.
     *
     * @return Computation finality.
     */
    bool finished() const {
        return approx_result.finished;
    }
};


/**
 * @brief Basic procedure of a block evaluation in scalar block format.
 *
 * @tparam DataType Block element type.
 */
template<typename DataType>
class Scalar: public Computer<DataType> {
private:
    using RealType = decltype(std::abs(DataType()));

    block::Scalar<DataType> *_block{nullptr}; ///< Scalar block to approximate

    DataType element{};
    Result approx_result{0.0, 0.0, 0, false}; ///< Current approximation state

    std::shared_ptr<log::Logger> log;

public:
    /**
     * @brief Initialize data to approximate.
     *
     * @param[in] data Diagonal block to compute.
     * @param[in] log_ptr Logger for computation procedure.
     */
    Scalar(block::Scalar<DataType> *data = nullptr,
           std::shared_ptr<log::Logger> log_ptr = nullptr) : _block{data}, log{std::move(log_ptr)} {}

    ~Scalar() override = default;

    Result compute(matrix::Matrix<DataType> &A,
                   const index::Index<uint64_t> &row_idx,
                   const index::Index<uint64_t> &col_idx) override {
        uint64_t m{row_idx.size()}, n{col_idx.size()};
        if (_block == nullptr) {
            _block = new block::Scalar<DataType>(m, n);
        } else if (std::make_pair(m, n) != _block->size()) {
            return approx_result = {0, 0, 0, false};
        }
        _block->set_data(A.element(row_idx[0], col_idx[0]));
        return approx_result = {_block->frob_norm(), 0.0, 1, true};
    }

    Result step_enqueue(eval::Request<DataType> &request) override {
        if (approx_result.finished) { // If already approximated do nothing
            request.free();
            return approx_result;
        }

        uint64_t m{}, n{};
        std::tie(m, n) = _block->size();
        if (approx_result.memory == 0) {
            // Request the single element
            request.element(0, 0, &element);
            return approx_result = {0, 0, 1, false};
        }

        // All elements are received
        _block->set_data(element);
        approx_result.matr_norm = _block->frob_norm();
        return approx_result;
    }

    /**
     * @brief Get current computation result.
     *
     * @return Intermediate Result.
     */
    Result result() const {
        return approx_result;
    }

    /**
     * @brief Check if computation has finished.
     *
     * @return Computation finality.
     */
    bool finished() const {
        return approx_result.finished;
    }
};


/**
 * @brief Basic procedure of a block evaluation in scalar block format.
 *
 * @tparam DataType Block element type.
 */
template<typename DataType>
class None: public Computer<DataType> {
private:
    using RealType = decltype(std::abs(DataType()));

    block::Data<DataType> *_block{nullptr}; ///< Diagonal block to approximate

    Result approx_result{0.0, 0.0, 0, true}; ///< Current epproximation state

    std::shared_ptr<log::Logger> log;

public:
    /**
     * @brief Initialize data to approximate.
     *
     * @param[in] data Any block that does not require computation.
     * @param[in] log_ptr Logger for computation procedure.
     */
    None(block::Data<DataType> *data = nullptr,
         std::shared_ptr<log::Logger> log_ptr = nullptr) : _block{data}, log{std::move(log_ptr)} {
        if (_block != nullptr) {
            uint64_t m, n;
            std::tie(m, n) = _block->size();
            approx_result.memory = _block->memory();
            approx_result.matr_norm = _block->frob_norm();
        }
    }

    ~None() override = default;

    Result compute(matrix::Matrix<DataType> &A,
                   const index::Index<uint64_t> &row_idx,
                   const index::Index<uint64_t> &col_idx) override {
        (void) A;
        (void) row_idx;
        (void) col_idx;
        return approx_result;
    }

    Result step_enqueue(eval::Request<DataType> &request) override {
        // Do nothing
        (void) request;
        return approx_result;
    }

    /**
     * @brief Get current computation result.
     *
     * @return Intermediate Result.
     */
    Result result() const {
        return approx_result;
    }

    /**
     * @brief Check if computation has finished.
     *
     * @return Computation finality.
     */
    bool finished() const {
        return approx_result.finished;
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


extern template class None<float>;

extern template class None<double>;

extern template class None<std::complex<float>>;

extern template class None<std::complex<double>>;


}
}


#endif //MOSAIC_SKELETON_COMPUTATION_DAIGONAL_HPP
