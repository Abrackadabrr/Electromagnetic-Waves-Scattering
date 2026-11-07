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
 * @brief Defines Conjugate Gradient method and its utilities.
 */


 #ifndef MOSAIC_SKELETON_SOLVERS_LINSYS_LEAST_SQUARES_HPP
 #define MOSAIC_SKELETON_SOLVERS_LINSYS_LEAST_SQUARES_HPP
 
 
 #include <cstdint>
 
 #include "smart_array/smart_array.hpp"
 #include "cpp_blas/BLAS.hpp"
 #include "cpp_lapack/LAPACK.hpp"

 #include "msk/utility/svd.hpp"
 
 
namespace MSk {
 
namespace solvers {
 
namespace linsys {


template<typename DataType>
void block_to_block(uint64_t n_rows, uint64_t n_cols, 
                    const DataType *block_in, uint64_t ld_block_in,
                    DataType *block_out, uint64_t ld_block_out) {
    if (n_rows == ld_block_in && n_rows == ld_block_out) {
        BLAS::copy(n_rows * n_cols, block_in, 1, block_out, 1);
    } else {
        for (uint64_t k = 0; k < n_cols; ++k) {
            BLAS::copy(n_rows, block_in + k * ld_block_in, 1, block_out + k * ld_block_out, 1);
        }
    }
}
 
// A X = B
// Q R X = B
// R X= Q^* B
// B = R^-1 Q^* B
template<typename DataType>
void lsq_qr(uint64_t n_rows, uint64_t n_cols, uint64_t n_rhs, 
            const DataType *matrix, uint64_t ld_matrix, 
            const DataType *rhs, uint64_t ld_rhs, 
            DataType *sol, uint64_t ld_sol) {
    if (n_rows < n_cols) {
        return;
    }
    uint64_t lwork = -1;
    int info{0};
    SmartArray<DataType> work(1);
    work[0] = 0;
    DataType *tmp{nullptr};
    LAPACK::geqrf(n_rows, n_cols, tmp, ld_matrix, tmp, work.data(), lwork, info);
    lwork = std::real(work[0]);
    work.resize(lwork);

    SmartArray<DataType> tau(std::min(n_rows, n_cols));
    SmartArray<DataType> ref(n_rows * n_cols);
    block_to_block(n_rows, n_cols, matrix, ld_matrix, ref.data(), n_rows);
    // Compute decomposition: A = Q R
    LAPACK::geqrf(n_rows, n_cols, ref.data(), n_rows, tau.data(), work.data(), lwork, info);

    SmartArray<DataType> rhs_proj(n_rows * n_rhs);
    block_to_block(n_rows, n_rhs, rhs, ld_rhs, rhs_proj.data(), n_rows);

    // Calculate Q^* B
    LAPACK::unmqr('L', 'C', n_rows, n_rhs, n_cols, ref.data(), n_rows, tau.data(), rhs_proj.data(), n_rows);

    // Solve R X = Q^* B
    BLAS::trsm('L', 'U', 'N', 'N', n_cols, n_rhs, DataType{1.0}, ref.data(), n_rows, rhs_proj.data(), n_rows);
    block_to_block(n_rows, n_cols, rhs_proj.data(), n_rows, sol, ld_sol);
}
 
 
template<typename DataType>
void lsq_svd(uint64_t n_rows, uint64_t n_cols, uint64_t n_rhs,
                const DataType *matrix, uint64_t ld_matrix,
                const DataType *rhs, uint64_t ld_rhs,
                DataType *sol, uint64_t ld_sol) {
    using RealType = decltype(std::abs(DataType()));
    SmartArray<DataType> U, V;
    SmartArray<RealType> S;
    SmartArray<DataType> raw_matrix(n_rows * n_cols);
    copy_matrix(n_rows, n_cols, matrix, ld_matrix, raw_matrix.data());
    auto res = perform_svd(n_rows, n_cols, raw_matrix, U, S, V, {0.0, 1e-8});

    uint64_t rank = res.memory / (n_rows + n_cols);

    SmartArray<DataType> UcB(rank * n_rhs);
    BLAS::gemm('C', 'N', rank, n_rhs, n_rows, DataType{1.0}, U.data(), n_rows, rhs, ld_rhs, DataType{0.0}, UcB.data(), rank);
    
    for (uint64_t i = 0; i < rank; ++i) {
        for (uint64_t k = 0; k < n_rhs; ++k) {
            UcB[i + k * rank] /= S[i];
        }
    }
    BLAS::gemm('N', 'N', n_cols, n_rhs, rank, DataType{1.0}, V.data(), n_cols, UcB.data(), rank, DataType{0}, sol, ld_sol);
}


 
}
}
}
 
 
#endif //MOSAIC_SKELETON_SOLVERS_LINSYS_LEAST_SQUARES_HPP
 