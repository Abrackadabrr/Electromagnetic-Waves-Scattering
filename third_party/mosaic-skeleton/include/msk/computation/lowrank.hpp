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


/**
 * @file
 *
 * @brief Defines low-rank block approximation procedures.
 */


#ifndef MOSAIC_SKELETON_COMPUTATION_LOWRANK_HPP
#define MOSAIC_SKELETON_COMPUTATION_LOWRANK_HPP


#include <cstdint>
#include <cstdlib>
#include <memory>
#include <random>
#include <vector>

#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"
#include "cpp_lapack/LAPACK.hpp"

#include "msk/computation/computer.hpp"
#include "msk/block/lowrank.hpp"
#include "msk/utility/logging.hpp"
#include "msk/utility/svd.hpp"
#include "msk/evaluation/request.hpp"


namespace MSk {

namespace comp {


/**
 * @brief Low-rank approximation of a block using singular value decomposition.
 *
 * @tparam DataType Block element type.
 */
template<typename DataType>
class SVD : public Computer<DataType> {
private:
    using RealType = decltype(std::abs(DataType()));

    const Parameters params; ///< Approximation parameters
    const bool _by_col{true}; ///< If to evaluate matrix element column by column

    block::LowRank<DataType> *_block{nullptr}; ///< Block to approximate
    Result approx_result{0.0, 0.0, 0, false}; ///< Approximation result

    // Temporary approximation data
    SmartArray<DataType> _data{0}; ///< Block in dense form
    uint64_t next_col{0}; ///< Next column to evaluate, if _by_col == true

    std::shared_ptr<log::Logger> log{nullptr}; ///< Logger

public:
    /**
     * @brief Initialize data to approximate and parameters.
     *
     * @param[in] data Low-rank block to approximate.
     * @param[in] parameters Approximation parameters.
     * @param[in] by_col Flag: request block's elements column-by-column or all at once.
     * @param[in] log_ptr Logger.
     */
    SVD(block::LowRank<DataType> *data = nullptr,
        const Parameters &parameters = {0, 1e-3, 100}, const bool by_col = true, 
        std::shared_ptr<log::Logger> log_ptr = nullptr) :
            Computer<DataType>{}, _block{data},
            params{parameters}, _by_col{by_col}, log{std::move(log_ptr)} {
        uint64_t m{}, n{};
        std::tie(m, n) = _block->size();
        _data.resize(m * n);
    }

    ~SVD() override = default;

    Result compute(matrix::Matrix<DataType> &A,
                   const index::Index<uint64_t> &row_idx,
                   const index::Index<uint64_t> &col_idx) override {
        uint64_t m{}, n{};
        std::tie(m, n) = _block->size();
        A.block(row_idx, col_idx, _data.data());
        SmartArray<DataType> col_factor, row_factor;
        // Perform SVD on the data computed
        approx_result = perform_svd(m, n, _data, col_factor, row_factor, params);
        _data.free();
        uint64_t r{approx_result.memory / (m + n)};
        _block->set_data(std::move(col_factor), std::move(row_factor), r, approx_result.matr_norm);
        return approx_result;
    }

    Result step_enqueue(eval::Request<DataType> &request) override {
        if (approx_result.finished) { // If already approximated do nothing
            request.free();
            return approx_result;
        }
        uint64_t m{}, n{};
        std::tie(m, n) = _block->size();

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
                return approx_result = {0, 0, next_col * m, false};
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
        SmartArray<DataType> col_factor, row_factor;
        approx_result = perform_svd(m, n, _data, col_factor, row_factor, params);
        _data.free();
        uint64_t r{approx_result.memory / (m + n)};
        _block->set_data(std::move(col_factor), std::move(row_factor), r, approx_result.matr_norm);

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

    /**
     * @brief Get current approximation result.
     *
     * @return Intermediate Result.
     */
    Result result() const {
        return approx_result;
    }
};


/**
 * @brief Low-rank approximation of a block using the matrix cross method.
 *
 * @tparam DataType Block element type.
 */
template<typename DataType>
class Cross : public Computer<DataType> {
private:
    using RealType = decltype(std::abs(DataType()));

    const Parameters params; ///< Parameters of approximation
    const uint64_t n_max_lookups{2}; ///< Number of row/columns lookups for maximum
    const RealType rho{1.0 + 1e-6}; ///< Look for rho-local maximal element at each row/column
    const bool _recomp{false}; ///< Perform recompression at the end of approximation

    uint64_t n_test_cols{4};
    SmartArray<DataType> test_cols; ///< Testing columns for 1st maximal element lookup
    SmartArray<uint64_t> test_cols_ids; ///< Testing columns indices sequence
    uint64_t new_tcol_id_idx{0}; ///< Index of next test column id to take
    bool no_more_test_candidates{false};
    uint64_t place_tcol_idx{n_test_cols}; ///< Index of new test column to be placed

    bool init{false}; ///< Current state
    uint64_t lookup{0}; ///< Lookup number
    RealType elem_max{0}; ///< Current local maximal element
    RealType prev_elem_max{0}; ///< Previous local maximal element
    uint64_t i_max{0}, j_max{0}; ///< Position of the maximal element

    uint64_t _rank{0}; ///< Current approximant rank
    SmartArray<DataType> U_factor, V_factor; ///< Approximation factors
    std::vector<uint64_t> rows_used; ///< Rows used for approximation
    std::vector<uint64_t> columns_used; ///< Columns used for approximation
    RealType UV_norm_sq{0}; ///< Squared norm of current approximation U V^T
    RealType est_err_norm{0}; ///< Current error norm estimation

    uint64_t n_rows, n_cols;
    block::LowRank<DataType> *_block{nullptr};
    Result approx_result{0.0, 0.0, 0, false};
    std::shared_ptr<log::Logger> log{nullptr};

private:
    void init_test_columns(eval::Request<DataType> &request) {
        test_cols_ids.resize(n_cols);
        for (uint64_t i = 0; i < n_cols; ++i) {
            test_cols_ids[i] = i;
        }
        std::random_device device;
        std::mt19937_64 engine(device());
        std::shuffle(test_cols_ids.data(), test_cols_ids.data() + n_cols, engine);
        new_tcol_id_idx = n_cols - 1;
        no_more_test_candidates = false;
        test_cols.resize(n_rows * n_test_cols);
        // Enqueue test cols
        std::vector<uint64_t> row_idx(n_rows), col_idx(n_test_cols);
        for (uint64_t i = 0; i < n_rows; ++i) {
            row_idx[i] = i;
        }
        for (uint64_t j = 0; j < n_test_cols; ++j) {
            col_idx[j] = test_cols_ids[j];
        }
        request.block(row_idx, col_idx, test_cols.data());
        place_tcol_idx = n_test_cols;
    }

    /**
     * @brief Enqueue new bunch of block elements to evaluate.
     *
     * @param[out] request Request to enqueue elements.
     */
    void enqueue(eval::Request<DataType> &request) {
        if (lookup % 2 == 1) {
            std::vector<uint64_t> col_idx(n_cols);
            for (uint64_t j = 0; j < n_cols; ++j) {
                col_idx[j] = j;
            }
            V_factor.resize((_rank + 1) * n_cols);
            request.row(i_max, col_idx, &V_factor[_rank * n_cols]);
        } else {
            std::vector<uint64_t> row_idx(n_rows);
            for (uint64_t i = 0; i < n_rows; ++i) {
                row_idx[i] = i;
            }
            U_factor.resize((_rank + 1) * n_rows);
            request.column(row_idx, j_max, &U_factor[_rank * n_rows]);
        }
    }

    /**
     * @brief Add new skeleton to the approximation.
     *
     */
    void add_skeleton() {
        const DataType elem{V_factor[j_max + _rank * n_cols]}; // maximal element
        const DataType col_coef{DataType{1.0} / std::sqrt(std::abs(elem))};
        const DataType row_coef{std::sqrt(std::abs(elem)) / elem};
        BLAS::scal(n_rows, col_coef, &U_factor[_rank * n_rows], 1);
        BLAS::scal(n_cols, row_coef, &V_factor[_rank * n_cols], 1);
        
        SmartArray<DataType> u_dot(_rank), v_dot(_rank);
        BLAS::gemv('C', n_rows, _rank, DataType{1.0}, U_factor.data(), n_rows,
                   &U_factor[_rank * n_rows], 1, DataType{0.0}, u_dot.data(), 1);
        BLAS::gemv('C', n_cols, _rank, DataType{1.0}, V_factor.data(), n_cols,
                   &V_factor[_rank * n_cols], 1, DataType{0.0}, v_dot.data(), 1);
        const DataType uv_new_dot = BLAS::dotu(_rank, u_dot.data(), 1, v_dot.data(), 1);
        const RealType u_new = BLAS::nrm2(n_rows, &U_factor[_rank * n_rows], 1);
        const RealType v_new = BLAS::nrm2(n_cols, &V_factor[_rank * n_cols], 1);
        UV_norm_sq += 2 * std::real(uv_new_dot) + u_new * u_new * v_new * v_new;

        // Perform rank-1 update of test columns with the last added skeleton
        SmartArray<DataType> coefs(n_test_cols);
        for (uint64_t j = 0; j < n_test_cols; ++j) {
            coefs[j] = V_factor[test_cols_ids[j] + _rank * n_cols];
        }
        BLAS::geru(n_rows, n_test_cols, DataType{-1.0},
                &U_factor[_rank * n_rows], 1,
                coefs.data(), 1, test_cols.data(), n_rows);

        ++_rank;
        approx_result.matr_norm = std::sqrt(UV_norm_sq);
        approx_result.memory = _rank * (n_rows + n_cols);

        rows_used.push_back(i_max);
        columns_used.push_back(j_max);
    }

    /**
     * @brief Finalize the approximation algorithm.
     */
    void finish() {
        U_factor.resize(n_rows * _rank);
        V_factor.resize(n_cols * _rank);
        // Conjugate rows basis V for complex data type
        if (!std::is_same<DataType, RealType>::value) {
            LAPACK::lacgv(V_factor.size(), V_factor.data(), 1);
        }
        _block->set_data(std::move(U_factor), std::move(V_factor), _rank, approx_result.matr_norm);
        // Perform recompression if required
        if (_recomp) {
            approx_result = _block->recompress(params);
        }
        approx_result.finished = true;
    }

    static std::tuple<uint64_t, uint64_t, RealType> find_maximal(const uint64_t line_size, const uint64_t n_lines, const DataType *lines) {
        RealType e_max{std::abs(lines[0])};
        uint64_t k_max{0};
        for (uint64_t k = 1; k < line_size * n_lines; ++k) {
            if (std::abs(lines[k]) > e_max) {
                k_max = k;
                e_max = std::abs(lines[k]);
            }
        }
        return {k_max % line_size, k_max / line_size, e_max};
    }

public:
    /**
     * @brief Initialize data to approximate and parameters.
     *
     * @param[in] data Low-rank block to approximate.
     * @param[in] parameters Approximation parameters.
     * @param[in] recompress Flag: perform recompression at the end of approximation or not.
     * @param[in] log_ptr Logger.
     */
    Cross(block::LowRank<DataType> *data = nullptr, Parameters parameters = {0, 1e-3, 100},
          bool recompress = false, std::shared_ptr<log::Logger> log_ptr = nullptr) :
            Computer<DataType>(), _block{data}, log{std::move(log_ptr)},
            params{parameters}, _recomp{recompress} {
        std::tie(n_rows, n_cols) = _block->size();
    }

    ~Cross() override = default;

    Result compute(matrix::Matrix<DataType> &A,
                   const index::Index<uint64_t> &row_idx,
                   const index::Index<uint64_t> &col_idx) override {
        uint64_t m{row_idx.size()}, n{col_idx.size()};
        SmartArray<DataType> block(m * n);
        A.block(row_idx, col_idx, block.data());
        approx_result = perform_svd(m, n, block, U_factor, V_factor, params);
        block.free();
        _rank = approx_result.memory / (m + n);
        _block->set_data(std::move(U_factor), std::move(V_factor), _rank, approx_result.matr_norm);
        return approx_result;
    }

    Result step_enqueue(eval::Request<DataType> &request) override {
        if (approx_result.finished) { // If already approximated do nothing
            request.free();
            return approx_result;
        }
        if (!init) {
            init_test_columns(request);
            init = true;
            return approx_result;
        }

        while (true) {
            if (lookup == 0) { // Check test columns
                // If new test column was evaluated, update it with current approximation
                if (_rank > 0 && place_tcol_idx < n_test_cols) {
                    BLAS::gemv('N', n_rows, _rank,
                            DataType{-1.0}, &U_factor[0], n_rows,
                            V_factor.data() + j_max, n_cols,
                            DataType{1.0}, &U_factor[_rank * n_rows], 1);
                    BLAS::copy(n_rows, &U_factor[_rank * n_rows], 1, &test_cols[place_tcol_idx * n_rows], 1);
                    place_tcol_idx = n_test_cols;
                }
                // Find the maximal element in testing columns
                prev_elem_max = 0;
                uint64_t j_test_idx;
                std::tie(i_max, j_test_idx, elem_max) = find_maximal(n_rows, n_test_cols, test_cols.data());
                j_max = test_cols_ids[j_test_idx];
                lookup = 1;
                enqueue(request); // Row evaluation
                return approx_result;
            }
            // Update newly evaluated row or column: a - U v or a^T - u^T V^T
            if (_rank > 0) {
                if (lookup % 2 == 1) { // Check row i_max
                    BLAS::gemv('N', n_cols, _rank,
                            DataType{-1.0}, V_factor.data(), n_cols,
                            U_factor.data() + i_max, n_rows,
                            DataType{1.0}, V_factor.data() + _rank * n_cols, 1);
                } else { // Check column j_max
                    BLAS::gemv('N', n_rows, _rank,
                            DataType{-1.0}, U_factor.data(), n_rows,
                            V_factor.data() + j_max, n_cols,
                            DataType{1.0}, U_factor.data() + _rank * n_rows, 1);
                }
            }
            // Look for local maximal element
            if (lookup < n_max_lookups && elem_max >= rho * prev_elem_max) {
                prev_elem_max = elem_max;
                uint64_t j_idx{0};
                if (lookup % 2 == 1) { // In row
                    std::tie(j_max, j_idx, elem_max) = find_maximal(n_cols, 1, &V_factor[_rank * n_cols]);
                    ++lookup;
                    // Check if column is in test list to skip evaluation
                    for (j_idx = 0; j_idx < n_test_cols && test_cols_ids[j_idx] != j_max; ++j_idx) {}
                    if (j_idx >= n_test_cols) {
                        enqueue(request); // Column evaluation
                        return approx_result;
                    }
                    // Required column need not to be evaluated
                    U_factor.resize((_rank + 1) * n_rows);
                    BLAS::copy(n_rows, &test_cols[j_idx * n_rows], 1, &U_factor[_rank * n_rows], 1);
                    if (lookup < n_max_lookups && elem_max >= rho * prev_elem_max) {
                        prev_elem_max = elem_max;
                        std::tie(i_max, j_idx, elem_max) = find_maximal(n_rows, 1, &U_factor[_rank * n_rows]);
                        ++lookup;
                        enqueue(request); // Row evaluation
                        return approx_result;
                    }
                } else { // In column
                    std::tie(i_max, j_idx, elem_max) = find_maximal(n_rows, 1, &U_factor[_rank * n_rows]);
                    ++lookup;
                    enqueue(request); // Row evaluation
                    return approx_result;
                }
            }
            // If no lookups needed, check for convergence
            RealType error{0}; // Approximation error estimate
            if (lookup % 2 == 1) {
                error = BLAS::nrm2(n_cols, &V_factor[_rank * n_cols], 1) * std::sqrt(n_rows - _rank);
            } else {
                error = BLAS::nrm2(n_rows, &U_factor[_rank * n_rows], 1) * std::sqrt(n_cols - _rank);
            }
            approx_result.err_norm = error;
            const RealType err_bound = params.abs_tol + params.rel_tol * std::sqrt(UV_norm_sq);
            const uint64_t rank_bound = std::min((params.max_rank > 0) ? params.max_rank : n_rows, std::min(n_rows, n_cols));
            // Finish approximation
            if (error < err_bound || _rank >= rank_bound) {
                finish();
                return approx_result;
            }
            // Add another cross skeleton
            add_skeleton();
            if (_rank >= rank_bound) {
                finish();
                return approx_result;
            }

            // Restart lookups from test columns
            lookup = 0;
            
            // Check: last added column can't be used for testing            
            for (place_tcol_idx = 0; place_tcol_idx < n_test_cols && test_cols_ids[place_tcol_idx] != j_max; ++place_tcol_idx) {}
            if (place_tcol_idx < n_test_cols) {  // Need to find new test column to replace id[place_tcol_idx]
                uint64_t replacement_idx = n_cols;
                while (!no_more_test_candidates && new_tcol_id_idx >= n_test_cols) {
                    const uint64_t candidate_idx = new_tcol_id_idx;
                    if (new_tcol_id_idx == 0) no_more_test_candidates = true;
                    else --new_tcol_id_idx;
                    if (std::find(columns_used.begin(), columns_used.end(),
                                  test_cols_ids[candidate_idx]) == columns_used.end()) {
                        replacement_idx = candidate_idx;
                        break;
                    }
                }
                if (replacement_idx < n_cols) {
                    std::swap(test_cols_ids[place_tcol_idx], test_cols_ids[replacement_idx]);
                    j_max = test_cols_ids[place_tcol_idx];  // New test column to be evaluated
                    enqueue(request);
                    return approx_result;
                } else { // No new test columns (all are already used)
                    // Reduce test columns
                    --n_test_cols;
                    std::swap(test_cols_ids[place_tcol_idx], test_cols_ids[n_test_cols]);
                    // Move spoiled column to the end
                    if (place_tcol_idx < n_test_cols) {
                        BLAS::swap(n_rows, &test_cols[place_tcol_idx * n_rows], 1, &test_cols[n_test_cols * n_rows], 1);
                    }
                    place_tcol_idx = n_test_cols;
                    if (n_test_cols == 0) {
                        finish();
                        return approx_result;
                    }
                }
            }
        }
    }
};


extern template class SVD<float>;

extern template class SVD<double>;

extern template class SVD<std::complex<float>>;

extern template class SVD<std::complex<double>>;


extern template class Cross<float>;

extern template class Cross<double>;

extern template class Cross<std::complex<float>>;

extern template class Cross<std::complex<double>>;


}
}


#endif //MOSAIC_SKELETON_COMPUTATION_LOWRANK_HPP

