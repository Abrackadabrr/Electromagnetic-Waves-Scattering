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
 * @brief Defines blocks computation basics.
 */


#ifndef MOSAIC_SKELETON_COMPUTATION_COMPUTER_HPP
#define MOSAIC_SKELETON_COMPUTATION_COMPUTER_HPP


#include "msk/matrix/matrix.hpp"
#include "msk/index/index.hpp"
#include "msk/evaluation/request.hpp"


namespace MSk {

namespace comp {


/**
 * @brief Parameters of a Block approximation.
 */
struct Parameters {
    double abs_tol{0}; ///< Absolute error tolerance.
    double rel_tol{0}; ///< Relative error tolerance.
    uint64_t max_rank{0}; ///< Rank bound.

    Parameters(double at = 0, double rt = 0, uint64_t mr = 0) : abs_tol{at}, rel_tol{rt}, max_rank{mr} {}
};


/**
 * @brief Result of a Block computation (approximation) procedure.
 */
struct Result {
    double matr_norm{0}; ///< Norm of approximated block.
    double err_norm{0}; ///< Absolute error norm estimate.
    uint64_t memory{0}; ///< Number of elements stored for approximation.
    bool finished{false}; ///< Flag of procedure finish.

    Result(double nrm = 0.0, double err = 0.0, uint64_t mem = 0, bool finish = false) :
            matr_norm{nrm}, err_norm{err}, memory{mem}, finished{finish} {}
};


/**
 * @brief Abstract computation procedure of an arbitrary matrix block in some format.
 *
 * @details Computer defines the interface to design methods of evaluation of blocks.
 *
 * @tparam DataType Element type of block.
 */
template<class DataType>
class Computer {
public:
    Computer() = default;

    virtual ~Computer() = default;

    /**
     * @brief Performs computation (approximation) procedure of block.
     *
     * @param[in] A Matrix that provides elements evaluation.
     * @param[in] row_idx Row indices of a block inside the whole matrix.
     * @param[in] col_idx Column indices of a block inside the whole matrix.
     * @return Approximation result.
     */
    virtual Result compute(matrix::Matrix<DataType> &A,
                           const index::Index<uint64_t> &row_idx,
                           const index::Index<uint64_t> &col_idx) = 0;

    /**
     * @brief Performs one step of computation (approximation) procedure.
     *
     * @details The procedure should be designed as step-by-step iterative method.
     * At each step some computational work is done, and then some matrix elements are requested.
     *
     * @param[out] request Request to put matrix elements that are required for the next step.
     * @return Intermediate computation result.
     */
    virtual Result step_enqueue(eval::Request<DataType> &request) = 0;

};


extern template class Computer<float>;

extern template class Computer<double>;

extern template class Computer<std::complex<float>>;

extern template class Computer<std::complex<double>>;


}
}


#endif //MOSAIC_SKELETON_COMPUTATION_COMPUTER_HPP
