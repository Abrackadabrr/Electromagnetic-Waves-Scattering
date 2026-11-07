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
// Created by b.valiakhmetov on 5/25/23.
//


/**
 * @file
 *
 * @brief Contains extended SVD routines .
 */


#ifndef MOSAIC_SKELETON_UTILITY_SVD_HPP
#define MOSAIC_SKELETON_UTILITY_SVD_HPP


#include <cstdint>
#include <cassert>

#include "cpp_lapack/LAPACK.hpp"
#include "cpp_blas/BLAS.hpp"

#include "msk/computation/computer.hpp"


namespace MSk {


template<typename RealType>
uint64_t truncate_singular_values(uint64_t n, const RealType *S, const comp::Parameters &params, const bool frob = true) {
    const RealType norm = frob ? BLAS::nrm2(n, S, 1) : S[0];
    const RealType bound = norm * params.rel_tol + params.abs_tol;
    const uint64_t max_rank{std::min(n, params.max_rank > 0 ? params.max_rank : n)};
    uint64_t rank{n};
    if (frob) {
        RealType next_err_sq = S[rank - 1] * S[rank - 1];
        while (rank > max_rank || std::sqrt(next_err_sq) < bound) {
            --rank;
            if (rank == 0) {
                break;
            }
            next_err_sq += S[rank - 1] * S[rank - 1];
        }
    } else {
        for (rank = n; rank > 0 && (rank > max_rank || S[rank - 1] < bound); --rank) {}
    }
    return rank;
}


/**
 * @brief Perform truncated SVD over a rectangular matrix.
 *
 * @details Procedure uses LAPACK::gesdd subroutine.
 * The decomposition is written in form A = U S V^*, and U, S, V are returned.
 *
 * @tparam DataType Matrix element type.
 * @param[in] m Number of matrix rows.
 * @param[in] n Number of matrix columns.
 * @param[in, out] data Matrix to process, on output its contents are destroyed.
 * @param[out] U Left singular vectors.
 * @param[out] S Singular values.
 * @param[out] V Right singular vectors.
 * @param[in] params Prameters of decomposition.
 * @return comp::Result decomposition information.
 */
template<typename DataType>
comp::Result perform_svd(uint64_t m, uint64_t n, SmartArray<DataType> &data,
                         SmartArray<DataType> &U, SmartArray<decltype(std::real(DataType{}))> &S, SmartArray<DataType> &V,
                         comp::Parameters params, bool frob_truncate = true) {
    using RealType = decltype(std::abs(DataType()));
    uint64_t min_mn = std::min(m, n);
    U.resize(m * min_mn);
    S.resize(min_mn);
    SmartArray<DataType> VT(n * min_mn);

    // Perform SVD
    LAPACK::gesdd<DataType>('S', m, n, data.data(), m, S.data(), U.data(), m, VT.data(), min_mn);

    // Do truncation by Frobenius norm error
    uint64_t rank = truncate_singular_values(min_mn, S.data(), params, frob_truncate);

    // Save the result
    U.resize(m * rank);
    V.resize(n * rank);

    RealType appr_norm = BLAS::nrm2(rank, S.data(), 1);
    RealType appr_err = BLAS::nrm2(min_mn - rank, S.data() + rank, 1);
    S.resize(rank);

    // Restore V by columns and conjugate
    for (uint64_t i = 0; i < rank; ++i) {
        BLAS::copy(n, VT.data() + i, min_mn, V.data() + i * n, 1);
    }
    LAPACK::lacgv(n * rank, V.data(), 1);

    return {appr_norm, appr_err, rank * (m + n), true};
}


/**
 * @brief Perform low-rank approximation of a rectangular matrix using truncated SVD.
 *
 * @details Procedure uses LAPACK::gesdd subroutine.
 * The approximation is written in form A = U V^*, and U, V are returned.
 *
 * @tparam DataType Matrix element type.
 * @param[in] m Number of matrix rows.
 * @param[in] n Number of matrix columns.
 * @param[in, out] data Matrix to process, on output its contents are destroyed.
 * @param[out] col_factor Basis of columns, U.
 * @param[out] row_factor Conjugated basis of rows, V.
 * @param[in] params Prameters of approximation.
 * @return comp::Result approximation information.
 */
template<typename DataType>
comp::Result perform_svd(uint64_t m, uint64_t n, SmartArray<DataType> &data,
                         SmartArray<DataType> &col_factor, SmartArray<DataType> &row_factor,
                         comp::Parameters params) {
    using RealType = decltype(std::abs(DataType()));
    uint64_t min_mn = std::min(m, n);
    SmartArray<RealType> sing(min_mn);
    comp::Result result = perform_svd<DataType>(m, n, data, col_factor, sing, row_factor, params);
    uint64_t rank = result.memory / (m + n);

    for (uint64_t i = 0; i < rank; ++i) {
        BLAS::rscal(n, sing[i], row_factor.data() + i * n, 1);
    }
    return result;
}

/**
 * @brief Recompress matrix of low-rank (skeleton) format.
 *
 * @details Input matricex is in form A = U V^*. Recompression is performed using
 * QR-decompositions of bases and truncated SVD for R_u R_v^*.
 *
 * @tparam DataType Matrix element type.
 * @param[in] m Number of rows.
 * @param[in] n Number of columns.
 * @param[in] rank Input number of columns/rows in bases.
 * @param[in, out] col_factor Basis of columns.
 * @param[in, out] row_factor Conjugated basis of rows.
 * @param[in] params Recompression parameters.
 * @param[in] orth_cols If basis of columns is already arthogonal.
 * @param[in] orth_rows If basis of rows is already arthogonal.
 * @return comp::Result recompression result.
 */
template<typename DataType>
comp::Result recompress_lowrank(uint64_t m, uint64_t n, uint64_t rank,
                                SmartArray<DataType> &col_factor, SmartArray<DataType> &row_factor,
                                const comp::Parameters &params,
                                const bool orth_cols = false, const bool orth_rows = false) {
    if (orth_cols) {
        assert(rank <= m && "Number of orthogonal columns can't exceed space dimension");
    }
    if (orth_rows) {
        assert(rank <= n && "Number of orthogonal rows can't exceed space dimension");
    }
    using RealType = decltype(std::abs(DataType()));
    // Do only truncation if both factors are already orthogonal
    if (orth_rows && orth_cols) {
        SmartArray<RealType> norms(rank);
        SmartArray<uint64_t> index(rank);
        for (uint64_t i = 0; i < rank; ++i) {
            index[i] = i;
            norms[i] = BLAS::nrm2(m, col_factor.data() + i * m, 1) * BLAS::nrm2(n, row_factor.data() + i * n, 1);
        }
        std::stable_sort(index.data(), index.data() + rank, [&norms](const uint64_t &i, const uint64_t &j) {return norms[i] > norms[j];});
        SmartArray<RealType> sing(rank);
        for (uint64_t i = 0; i < rank; ++i) {
            sing[i] = norms[index[i]];
        }
        uint64_t new_rank = truncate_singular_values(rank, sing.data(), params, true);
        RealType new_nrm = BLAS::nrm2(new_rank, sing.data(), 1);
        RealType new_err = BLAS::nrm2(rank - new_rank, sing.data() + new_rank, 1);
        // Gather new basis vectors to factors' front

        for (uint64_t i_put{new_rank}, i_get{0}; i_put < rank && i_get < new_rank; ++i_put, ++i_get) {
            while (i_put < rank && index[i_put] >= new_rank) {
                ++i_put;
            }
            if (i_put >= rank) {
                break;
            }
            // (i_put >= new_rank && index[i_put] < new_rank) => then vector at index[i_put] can be overwritten
            while (i_get < new_rank && index[i_get] < new_rank) {
                ++i_get;
            }
            if (i_get >= new_rank) {
                break;
            }
            // (i_get < new_rank && index[i_get] >= new_rank) => vector at index[i_get] should be moved to new basis
            BLAS::copy(m, col_factor.data() + index[i_get] * m, 1, col_factor.data() + index[i_put] * m, 1);
            BLAS::copy(n, row_factor.data() + index[i_get] * n, 1, row_factor.data() + index[i_put] * n, 1);
        }
        return {new_nrm, new_err, new_rank * (m + n), true};
    }

    constexpr DataType zero{0.0};
    const uint64_t min_mr{std::min(m, rank)}, min_nr{std::min(n, rank)};

    // Compute work size for qr factorization
    uint64_t lwork = -1;
    int info{0};
    SmartArray<DataType> work(2);
    work[0] = work[1] = 0;
    DataType *tmp{nullptr};
    SmartArray<DataType> tau_u, tau_v;
    if (!orth_cols) {
        LAPACK::geqrf(m, rank, tmp, m, tmp, &work[0], lwork, info);
        tau_u.resize(min_mr);
    }
    if (!orth_rows) {
        LAPACK::geqrf(n, rank, tmp, n, tmp, &work[1], lwork, info);
        tau_v.resize(min_nr);
    }
    lwork = std::max(std::real(work[0]), std::real(work[1]));
    work.resize(lwork);

    // Calculate Ru Rv^*
    SmartArray<DataType> RuRv(min_mr * min_nr);
    BLAS::copy(min_mr * min_nr, &zero, 0, RuRv.data(), 1);
    if (!orth_cols && orth_rows) {
        // Compute QR for U: U = Qu Ru
        LAPACK::geqrf(m, rank, col_factor.data(), m, tau_u.data(), work.data(), lwork, info);
        // Extract Ru and scale its columns by norms of V'th columns
        for (uint64_t j = 0; j < rank; ++j) { // rank <= n is guaranteed
            BLAS::copy(std::min(j + 1, m), col_factor.data() + j * m, 1, RuRv.data() + j * min_mr, 1);
            RealType norm = BLAS::nrm2(n, row_factor.data() + j * n, 1);
            if (std::abs(norm - 1.0) > 5 * std::numeric_limits<RealType>::epsilon()) {
                BLAS::rscal(std::min(j + 1, m), norm, RuRv.data() + j * min_mr, 1);
                BLAS::rscal(n, RealType{1.0} / norm, row_factor.data() + j * n, 1);
            }
        }
    } else if (orth_cols && !orth_rows) {
        // Compute QR for V: V = Qv Rv
        LAPACK::geqrf(n, rank, row_factor.data(), n, tau_v.data(), work.data(), lwork, info);
        // Extract Rv* and multiply it by norms of U'th columns
        for (uint64_t j = 0; j < rank; ++j) { // rank <= m is guaranteed
            BLAS::copy(std::min(j + 1, n), row_factor.data() + j * n, 1, RuRv.data() + j, min_mr);
            RealType norm = BLAS::nrm2(m, col_factor.data() + j * m, 1);
            if (std::abs(norm - 1.0) > 5 * std::numeric_limits<RealType>::epsilon()) {
                BLAS::rscal(std::min(j + 1, n), norm, RuRv.data() + j, min_mr);
                BLAS::rscal(m, RealType{1.0} / norm, col_factor.data() + j * m, 1);
            }
        }
        LAPACK::lacgv(min_mr * min_nr, RuRv.data(), 1);
    } else { // !orth_cols && !orth_rows
        // Compute QR for V: V = Qv Rv
        LAPACK::geqrf(n, rank, row_factor.data(), n, tau_v.data(), work.data(), lwork, info);
        // Compute QR for U: U = Qu Ru
        LAPACK::geqrf(m, rank, col_factor.data(), m, tau_u.data(), work.data(), lwork, info);

        if (n >= rank) { // Rv is triangular, not trapezoidal
            // Extract Ru
            for (uint64_t j = 0; j < rank; ++j) {
                BLAS::copy(std::min(j + 1, m), col_factor.data() + j * m, 1, RuRv.data() + j * min_mr, 1);
            }
            // Calculate Ru Rv^*
            BLAS::trmm('R', 'U', 'C', 'N', min_mr, min_nr, DataType{1.0},
                       row_factor.data(), n, RuRv.data(), min_mr);
        } else if (m >= rank) {
            // Extract Rv^*
            for (uint64_t j = 0; j < rank; ++j) {
                BLAS::copy(std::min(j + 1, n), row_factor.data() + j * n, 1, RuRv.data() + j, min_mr);
            }
            LAPACK::lacgv(min_mr * min_nr, RuRv.data(), 1);
            BLAS::trmm('L', 'U', 'N', 'N', min_mr, min_nr, DataType{1.0},
                       col_factor.data(), m, RuRv.data(), min_mr);
        } else { // m < rank && n < rank
            if (m >= n) {
                // Extract part (m rows of r) of Rv^*
                for (uint64_t j = 0; j < m; ++j) {
                    BLAS::copy(std::min(j + 1, n), row_factor.data() + j * n, 1, RuRv.data() + j, m);
                }
                LAPACK::lacgv(m * n, RuRv.data(), 1);
                // Multiply by triangular part of Ru
                BLAS::trmm('L', 'U', 'N', 'N', m, n, DataType{1.0},
                           col_factor.data(), m, RuRv.data(), m);
            } else { // n > m
                // Extract part (n columns of r) of Ru
                for (uint64_t j = 0; j < n; ++j) {
                    BLAS::copy(std::min(j + 1, m), col_factor.data() + j * m, 1, RuRv.data() + j * m, 1);
                }
                // Multiply by triangular part of Rv^*
                BLAS::trmm('R', 'U', 'C', 'N', m, n, DataType{1.0},
                           row_factor.data(), n, RuRv.data(), m);
            }
            // Multiply remaining rectangular parts of Ru and Rv^*
            uint64_t k{std::max(m, n)};
            BLAS::gemm('N', 'C', m, n, rank - k, DataType{1.0}, col_factor.data() + k * m, m,
                       row_factor.data() + k * n, n, DataType{1.0}, RuRv.data(), m);
        }
    }

    // Apply SVD to Ru Rv^* and get new factors: Ru Rv^* = Xu Xv^*
    SmartArray<DataType> Xu, Xv;
    comp::Result result = perform_svd(min_mr, min_nr, RuRv, Xu, Xv, params);

    // Build new factors: U = Qu Xu, V = Qv Xv
    uint64_t new_rank{result.memory / (min_mr + min_nr)};
    SmartArray<DataType> new_col_factor(m * new_rank), new_row_factor(n * new_rank);

    lwork = std::max(m, n) * std::max(min_mr, min_nr);
    if (work.size() < lwork) {
        work.resize(lwork);
    }
    if (!orth_cols) { // Qu consists of reflections, apply them to Xu
        BLAS::copy(m * new_rank, &zero, 0, new_col_factor.data(), 1);
        for (uint64_t j = 0; j < new_rank; ++j) {
            BLAS::copy(min_mr, Xu.data() + j * min_mr, 1, new_col_factor.data() + j * m, 1);
        }
        LAPACK::unmqr('L', 'N', m, new_rank, min_mr, col_factor.data(), m,
                      tau_u.data(), new_col_factor.data(), m);
    } else { // Qu = U, do multiplication
        BLAS::gemm('N', 'N', m, new_rank, min_mr, DataType{1.0}, col_factor.data(), m,
                   Xu.data(), min_mr, DataType{0.0}, new_col_factor.data(), m);
    }
    if (!orth_rows) {
        BLAS::copy(n * new_rank, &zero, 0, new_row_factor.data(), 1);
        for (uint64_t j = 0; j < new_rank; ++j) {
            BLAS::copy(min_nr, Xv.data() + j * min_nr, 1, new_row_factor.data() + j * n, 1);
        }
        LAPACK::unmqr('L', 'N', n, new_rank, min_nr, row_factor.data(), n,
                      tau_v.data(), new_row_factor.data(), n);
    } else { // Qu = U, do multiplication
        BLAS::gemm('N', 'N', n, new_rank, min_nr, DataType{1.0}, row_factor.data(), n,
                   Xv.data(), min_nr, DataType{0.0}, new_row_factor.data(), n);
    }

    col_factor = std::move(new_col_factor);
    row_factor = std::move(new_row_factor);

    result.memory = new_rank * (m + n);
    return result;
}


extern template uint64_t truncate_singular_values(uint64_t n, const float *S, const comp::Parameters &params, const bool frob);

extern template uint64_t truncate_singular_values(uint64_t n, const double *S, const comp::Parameters &params, const bool frob);


extern template comp::Result perform_svd(uint64_t m, uint64_t n, SmartArray<float> &data,
         SmartArray<float> &U, SmartArray<float> &S, SmartArray<float> &V,
         comp::Parameters params, bool frob_truncate);

extern template comp::Result perform_svd(uint64_t m, uint64_t n, SmartArray<double> &data,
         SmartArray<double> &U, SmartArray<double> &S, SmartArray<double> &V,
         comp::Parameters params, bool frob_truncate);

extern template comp::Result perform_svd(uint64_t m, uint64_t n, SmartArray<std::complex<float>> &data,
         SmartArray<std::complex<float>> &U, SmartArray<float> &S, SmartArray<std::complex<float>> &V,
         comp::Parameters params, bool frob_truncate);

extern template comp::Result perform_svd(uint64_t m, uint64_t n, SmartArray<std::complex<double>> &data,
         SmartArray<std::complex<double>> &U, SmartArray<double> &S, SmartArray<std::complex<double>> &V,
         comp::Parameters params, bool frob_truncate);


extern template comp::Result perform_svd(uint64_t m, uint64_t n, SmartArray<float> &data,
         SmartArray<float> &col_factor, SmartArray<float> &row_factor,
         comp::Parameters params);

extern template comp::Result perform_svd(uint64_t m, uint64_t n, SmartArray<double> &data,
         SmartArray<double> &col_factor, SmartArray<double> &row_factor,
         comp::Parameters params);

extern template comp::Result perform_svd(uint64_t m, uint64_t n, SmartArray<std::complex<float>> &data,
         SmartArray<std::complex<float>> &col_factor, SmartArray<std::complex<float>> &row_factor,
         comp::Parameters params);

extern template comp::Result perform_svd(uint64_t m, uint64_t n, SmartArray<std::complex<double>> &data,
         SmartArray<std::complex<double>> &col_factor, SmartArray<std::complex<double>> &row_factor,
         comp::Parameters params);


extern template comp::Result recompress_lowrank(uint64_t m, uint64_t n, uint64_t rank,
        SmartArray<float> &col_factor, SmartArray<float> &row_factor,
        const comp::Parameters &params, const bool orth_cols, const bool orth_rows);

extern template comp::Result recompress_lowrank(uint64_t m, uint64_t n, uint64_t rank,
        SmartArray<double> &col_factor, SmartArray<double> &row_factor,
        const comp::Parameters &params, const bool orth_cols, const bool orth_rows);

extern template comp::Result recompress_lowrank(uint64_t m, uint64_t n, uint64_t rank,
        SmartArray<std::complex<float>> &col_factor, SmartArray<std::complex<float>> &row_factor,
        const comp::Parameters &params, const bool orth_cols, const bool orth_rows);

extern template comp::Result recompress_lowrank(uint64_t m, uint64_t n, uint64_t rank,
        SmartArray<std::complex<double>> &col_factor, SmartArray<std::complex<double>> &row_factor,
        const comp::Parameters &params, const bool orth_cols, const bool orth_rows);

}


#endif //MOSAIC_SKELETON_UTILITY_SVD_HPP
