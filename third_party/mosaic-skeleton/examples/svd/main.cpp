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
// Created by b.valiakhmetov on 5/26/23.
//


#include <iostream>
#include <vector>
#include <random>

#include "toml.hpp"
#include "common.hpp"

#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"

#include "msk/block/lowrank.hpp"
#include "msk/utility/testing.hpp"


int main(int argc, char *argv[]) {
    std::string problem_file;
    if (argc > 1) {
        problem_file = argv[1];
    } else {
        problem_file = std::string(MSK_EXAMPLES_DIRECTORY) + "/svd/parameters.toml";
    }
    auto params = toml::parse(problem_file);
    auto &matrix_params = toml::find(params, "matrix");
    uint64_t n_rows = toml::find<toml::integer>(matrix_params, "n_rows");
    uint64_t n_cols = toml::find<toml::integer>(matrix_params, "n_cols");
    uint64_t n_blocks = toml::find<toml::integer>(matrix_params, "n_blocks");
    uint64_t rank = toml::find<toml::integer>(matrix_params, "rank");

    if (rank < 1 || rank > n_rows || rank > n_cols) {
        rank = std::min(2ul, std::min(n_rows, n_cols));
    }
    if (n_blocks < 2) {
        n_blocks = 2;
    }

    using DataType = double;
    using RealType = decltype(std::abs(DataType{1.0}));

    constexpr DataType zero{0};

    std::normal_distribution<double> normal(0, 1);
    std::default_random_engine eng(42);

    uint64_t b_n_rows{n_rows / n_blocks};
    std::vector<const MSk::block::LowRank<DataType>*> lrank_blocks;
    SmartArray<DataType> true_matrix(n_rows * n_cols);
    BLAS::copy(n_rows * rank, &zero, 0, true_matrix.data(), 1);

    // Generate block matrix
    for (uint64_t i = 0; i < n_blocks; ++i) {
        auto b = new MSk::block::LowRank<DataType>(b_n_rows, n_cols);
        SmartArray<DataType> u(b_n_rows * rank), v(n_cols * rank);

        for (uint64_t j = 0; j < rank; ++j) {
            for (uint64_t k = 0; k < n_cols; ++k) {
                v[k + n_cols * j] = normal(eng);
            }
        }
        for (uint64_t j = 0; j < rank; ++j) {
            for (uint64_t k = 0; k < b_n_rows; ++k) {
                u[k + b_n_rows * j] = normal(eng);
            }
        }
        for (uint64_t j = 0; j < rank; ++j) {
            for (uint64_t k = 0; k < j; ++k) {
                DataType proj = BLAS::dotc(b_n_rows, &u[k * b_n_rows], 1, &u[j * b_n_rows], 1);
                BLAS::axpy(b_n_rows, -proj, &u[k * b_n_rows], 1, &u[j * b_n_rows], 1);
            }
            RealType nrm = BLAS::nrm2(b_n_rows, &u[j * b_n_rows], 1);
            BLAS::rscal(b_n_rows, 1.0 / nrm, &u[j * b_n_rows], 1);
        }

        BLAS::gemm('N', 'C', b_n_rows, n_cols, rank,
                   DataType{1.0},u.data(), b_n_rows, v.data(), n_cols,
                   DataType{0.0}, true_matrix.data() + b_n_rows * i, n_rows);

        b->set_data(std::move(u), std::move(v), rank);
        lrank_blocks.push_back(b);
    }

    // Get full SVD (true)
    SmartArray<DataType> Ut(n_rows * n_rows);
    SmartArray<DataType> Vt(n_cols * n_cols);
    SmartArray<RealType> St(std::min(n_rows, n_cols));
    SmartArray<DataType> copy_matrix = true_matrix;
    LAPACK::gesdd('S', n_rows, n_cols, copy_matrix.data(), n_rows,
                  St.data(), Ut.data(), n_rows, Vt.data(), n_cols);
    for (uint64_t i = 0; i < n_cols; ++i) {
        for (uint64_t j = i + 1; j < n_cols; ++j) {
            std::swap(Vt[i + n_cols * j], Vt[j + n_cols * i]);
        }
    }

    SmartArray<DataType> U;
    SmartArray<DataType> V;
    SmartArray<RealType> S;
    uint64_t new_rank = MSk::block_svd(lrank_blocks.size(), lrank_blocks, U, S, V);
    std::cout << "Rank = " << new_rank << '\n';

    // Check factors error
    SmartArray<DataType> UU = Ut, VV = Vt;
    SmartArray<RealType> SS = St;
    for (uint64_t i = 0; i < rank; ++i) {
        DataType factor = UU[i * n_rows] / U[i * n_rows];
        factor /= std::abs(factor);
        BLAS::axpy(n_rows, DataType(-1.0) * factor, &U[i * n_rows], 1, &UU[i * n_rows], 1);
    }
    double u_err = BLAS::nrm2(n_rows * rank, UU.data(), 1);
    double u_nrm = BLAS::nrm2(n_rows * rank, Ut.data(), 1);
    std::cout << "U_err = " << u_err / u_nrm << std::endl;

    for (uint64_t i = 0; i < rank; ++i) {
        DataType factor = Vt[i * n_cols] / V[i * n_cols];
        factor /= std::abs(factor);
        BLAS::axpy(n_cols, DataType(-1.0) * factor, &V[i * n_cols], 1, &VV[i * n_cols], 1);
    }
    double v_err = BLAS::nrm2(n_cols * rank, VV.data(), 1);
    double v_nrm = BLAS::nrm2(n_cols * rank, Vt.data(), 1);
    std::cout << "V_err = " << v_err / v_nrm << std::endl;

    BLAS::axpy(rank, RealType(-1.0), S.data(), 1, SS.data(), 1);

    double s_err = BLAS::nrm2(rank, SS.data(), 1);
    double s_nrm = BLAS::nrm2(rank, St.data(), 1);
    std::cout << "S_err = " << s_err / s_nrm << std::endl;

    // Check full error
    SmartArray<DataType> check(n_rows * n_cols);
    for (uint64_t i = 0; i < new_rank; ++i) {
        BLAS::rscal(n_rows, S[i], U.data() + n_rows * i, 1);
    }
    BLAS::gemm('N', 'C', n_rows, n_cols, new_rank, DataType{1.0}, U.data(), n_rows,
               V.data(), n_cols, DataType{0.0}, check.data(), n_rows);

    BLAS::axpy(n_rows * n_cols, DataType(-1.0), true_matrix.data(), 1, check.data(), 1);

    RealType err = BLAS::nrm2(n_rows * n_cols, check.data(), 1);
    RealType nrm = BLAS::nrm2(n_rows * n_cols, true_matrix.data(), 1);

    std::cout << "Block-SVD err = " << err / nrm << std::endl;

    for (uint64_t i = 0; i < n_blocks; ++i) {
        delete lrank_blocks[i];
    }
}
