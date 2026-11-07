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
 * @brief Defines dense block computation procedure.
 */


#ifndef MOSAIC_SKELETON_COMPUTATION_DENSE_HPP
#define MOSAIC_SKELETON_COMPUTATION_DENSE_HPP


#include <cmath>

#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"

#include "msk/computation/computer.hpp"
#include "msk/block/dense.hpp"
#include "msk/utility/logging.hpp"


namespace MSk {

namespace comp {


/**
 * @brief Basic procedure of a block evaluation in dense format.
 *
 * @tparam DataType Block element type.
 */
template<typename DataType>
class Dense: public Computer<DataType> {
private:
    using RealType = decltype(std::abs(DataType()));

    block::Dense<DataType> *_block{nullptr}; ///< Dense block to approximate

    Result approx_result{0.0, 0.0, 0, false}; ///< Current epproximation state

    SmartArray<DataType> _data{0};
    bool _by_col{true}; ///< Flag: evaluate column-by-column or at one step (in case of step-by-step approximation)
    uint64_t next_col{0}; ///< Next column to evaluate in this case

    std::shared_ptr<log::Logger> log;

public:
    /**
     * @brief Initialize data to approximate.
     *
     * @param[in] data Dense block to compute.
     * @param[in] by_col Flag: to evaluate column-by-column or not.
     * @param[in] log_ptr Logger for computation procedure.
     */
    Dense(block::Dense<DataType> *data = nullptr, bool by_col = true,
          std::shared_ptr<log::Logger> log_ptr = nullptr) : _block{data}, _by_col{by_col},
                                                               log{std::move(log_ptr)} {}

    ~Dense() override = default;

    Result compute(matrix::Matrix<DataType> &A,
                   const index::Index<uint64_t> &row_idx,
                   const index::Index<uint64_t> &col_idx) override {
        uint64_t m{row_idx.size()}, n{col_idx.size()};
        if (_block == nullptr) {
            _block = new block::Dense<DataType>(m, n);
        } else if (std::make_pair(m, n) != _block->size()) {
            return approx_result = {0, 0, 0, false};
        }
        _data.resize(m * n);
        A.block(row_idx, col_idx, _data.data());
        RealType f_norm = BLAS::nrm2(m * n, _data.data(), 1);
        _block->set_data(std::move(_data), f_norm);
        return approx_result = {f_norm, 0.0, m * n, true};
    }

    Result step_enqueue(eval::Request<DataType> &request) override {
        if (approx_result.finished) { // If already approximated do nothing
            request.free();
            return approx_result;
        }

        uint64_t m{}, n{};
        std::tie(m, n) = _block->size();
        if (next_col == 0) {
            _data.resize(m * n);
        }

        // Request another bunch of elements
        if (next_col < n) { // If there are elements to evaluate
            if (!_by_col) { // Request all at once
                std::vector<uint64_t> row_idx(m), col_idx(n);
                for (uint64_t i = 0; i < m; ++i) {
                    row_idx[i] = i;
                }
                for (uint64_t j = 0; j < n; ++j) {
                    col_idx[j] = j;
                }
                request.block(row_idx, col_idx, _data.data());
                next_col = n;
                return approx_result = {0, 0, m * next_col, false};
            } else { // Request next column
                std::vector<uint64_t> row_idx(m);
                for (uint64_t i = 0; i < m; ++i) {
                    row_idx[i] = i;
                }
                request.column(row_idx, next_col, _data.data() + next_col * m);
                ++next_col;
                return approx_result = {0, 0, 0, false};
            }
        }

        // All elements are received
        RealType f_norm = BLAS::nrm2(m * n, _data.data(), 1);
        _block->set_data(std::move(_data), f_norm);
        return approx_result = {f_norm, 0.0, m * n, true};
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


extern template class Dense<float>;

extern template class Dense<double>;

extern template class Dense<std::complex<float>>;

extern template class Dense<std::complex<double>>;


}
}


#endif //MOSAIC_SKELETON_COMPUTATION_DENSE_HPP
