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
// Created by bulat_v on 25.06.23.
//


/**
 * @file
 *
 * @brief Defines procedures for mosaic-skeleton approximation post-processing.
 */


#ifndef MOSAIC_SKELETON_UTILITY_TESTING_HPP
#define MOSAIC_SKELETON_UTILITY_TESTING_HPP


#include "msk_config.h"

#ifdef MSK_USE_MPI
#include "cpp_mpi/MPI.hpp"
#endif


#include <cstdint>
#include <iostream>
#include <string>
#include <chrono>
#include <random>


#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"
#include "cpp_lapack/LAPACK.hpp"

#include "msk/block/block.hpp"
#include "msk/computation/dense.hpp"
#include "msk/block/dense.hpp"
#include "msk/block/lowrank.hpp"
#include "msk/utility/logging.hpp"
#include "msk/mosaic_skeleton.hpp"
#include "msk/solvers/linsys/GMRES.hpp"
#include "msk/utility/convolution.hpp"


namespace MSk {


/**
 * @brief Basic rectangular matrix block printing to standard output.
 *
 * @details Matrix is expected to be stored by columns.
 *
 * @tparam DataType Matrix element type.
 * @param[in] m Number of rows.
 * @param[in] n Number of columns.
 * @param[in] ptr Matrix elements pointer.
 * @param[in] ldim Leading dimension of the block.
 * @param[in] pre String to be printed right before the matrix.
 */
template<class DataType>
void out_matrix(uint64_t m, uint64_t n, DataType *ptr, uint64_t ldim, const std::string &pre = "") {
    std::cout << pre;
    for (uint64_t i = 0; i < m; ++i) {
        for (uint64_t j = 0; j < n; ++j) {
            std::cout << ptr[i + ldim * j] << " ";
        }
        std::cout << std::endl;
    }
}


/**
 * @brief Solve a system of linear equations with multiple right-hand sides with a matrix in
 * mosaic-skeleton format using GMRES method.
 *
 * @tparam DataType Matrix element type.
 * @param[in] msk_matrix Matrix in mosaic-skeleton format.
 * @param[in] n_rhs Number of system's right-hand sides.
 * @param[in, out] rhs Right-hand sides array on input, residuals array on output.
 * @param[in] ld_rhs Leading dimension of RHS array.
 * @param[in, out] sol Initial guess array on input, solutions array on output.
 * @param[in] ld_sol Leading dimension of solutions array.
 * @param[in] params Solving parameters.
 * @param[in] log Logger.
 */
template<class DataType, Trans trans = Trans::gen>
void run_gmres(const MosaicSkeleton<DataType, trans> &msk_matrix,
               uint64_t n_rhs, DataType *rhs, uint64_t ld_rhs, DataType *sol, uint64_t ld_sol,
               solvers::linsys::Parameters<decltype(std::real(DataType()))> params = {},
               std::shared_ptr<log::Logger> log = nullptr) {
    for (uint64_t j = 0; j < n_rhs; ++j) {
        msk_matrix.permute_row(rhs + j * ld_rhs);
        msk_matrix.permute_col(sol + j * ld_sol);
    }
    auto t1 = std::chrono::steady_clock::now();
    solvers::linsys::GMRES(msk_matrix, n_rhs, rhs, ld_rhs, sol, ld_sol, params);
    auto t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> delta = t2 - t1;
    double full_time = delta.count();
    if (log) log->info("GMRES time = {} on {} rhs", full_time, n_rhs);
    for (uint64_t j = 0; j < n_rhs; ++j) {
        msk_matrix.unpermute_row(rhs + j * ld_rhs);
        msk_matrix.unpermute_col(sol + j * ld_sol);
    }
}


/**
 * @brief Log compression rate of a matrix in mosaic-skeleton format.
 *
 * @details In case of multiprocessing mode also checks memory imbalance between processes.
 *
 * @tparam DataType Matrix element type.
 * @param[in] msk Matrix in mosaic-skeleton format.
 * @param[in] log Logger.
 */
template<class DataType, Trans trans = Trans::gen>
double check_compression(const MosaicSkeleton<DataType, trans> &msk, std::shared_ptr<log::Logger> log = nullptr) {
    uint64_t n_rows, n_cols;
    std::tie(n_rows, n_cols) = msk.size();
    uint64_t n_elements = n_rows * n_cols;

    int proc_id{0};
#ifdef MSK_USE_MPI
    MPI_Comm comm = msk.comm();
    MPI_Comm_rank(comm, &proc_id);
    int n_proc{1};
    MPI_Comm_size(comm, &n_proc);
#endif

    uint64_t global_mem{msk.memory()};
#ifdef MSK_USE_MPI
    uint64_t local_mem{msk.memory(false)};
    if (log) log->info("Local memory = {:3.4}%", (100.0 * local_mem) / global_mem);
    double imbalance{100.0 * std::abs((1.0 * local_mem * n_proc) / global_mem - 1.0)};
    if (proc_id == 0) {
        MPI_Reduce(MPI_IN_PLACE, &imbalance, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (log) log->info("Max memory imbalance = {:3.4}%", imbalance);
    } else {
        MPI_Reduce(&imbalance, &imbalance, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    }
#endif
    double compression = (1.0 * global_mem) / n_elements;
    if (proc_id == 0) {
        if (log) log->info("Compression is {:3.4}%", 100.0 * compression);
    }
    return compression;
}


/**
 * @brief Log approximation relative error and matrix norm of a mosaic-skeleton matrix.
 *
 * @tparam DataType Matrix element type.
 * @tparam Matrix Functional matrix type of MosaicSkeleton.
 * @param[in] msk Matrix in mosaic-skeleton format.
 * @param[in] log Logger.
 */
template<class DataType, class Matrix, Trans trans = Trans::gen>
double check_error(MosaicSkeleton<DataType, trans> &msk, Matrix &matrix, std::shared_ptr<log::Logger> log = nullptr) {
    double f_norm{msk.frob_norm()};
    auto t_start = std::chrono::steady_clock::now();
    double f_error{msk.check_error(matrix)};
    auto t_end = std::chrono::steady_clock::now();
    double t_diff = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    if (log) log->debug("MSk error check time = {:.4f}ms", t_diff);
    int proc_id{0};
#ifdef MSK_USE_MPI
    MPI_Comm comm = msk.comm();
    MPI_Comm_rank(comm, &proc_id);
#endif
    if (proc_id == 0) {
        if (log) log->info("F_error = {:.4e}", f_error);
        if (log) log->info("F_norm = {:.4e}", f_norm);
    }
    return f_error;
}


/**
 * @brief Log mosaic-skeleton matrix-vector product performance time and error optionally.
 *
 * @details If to compare to dense matrix then it is acllocated and calculated.
 * In case of multiprocessing mode it is done on every process, so think of your memory limits.
 *
 * @tparam DataType Matrix element type.
 * @tparam Matrix Functional matrix type of MosaicSkeleton.
 * @param[in] msk Matrix in mosaic-skeleton format.
 * @param[in] log Logger.
 * @param[in] compare_to_dense If to calculate error in comperison to a dense matrix.
 * @param[in] n_checks How many runs to average.
 */
template<class DataType, class Matrix, Trans trans = Trans::gen>
void check_matvec(MosaicSkeleton<DataType, trans> &msk, Matrix &matrix, std::shared_ptr<log::Logger> log,
                  bool compare_to_dense = true, int n_checks = 3) {
    if (n_checks < 1) {
        return;
    }

    uint64_t n_rows, n_cols;
    std::tie(n_rows, n_cols) = msk.size();

    int proc_id{0};
#ifdef MSK_USE_MPI
    MPI_Comm comm = msk.comm();
    MPI_Comm_rank(comm, &proc_id);
#endif

    // Calculate true dense matrix if required
    block::Block<DataType> *full_block{nullptr};
    block::Dense<DataType> *data_block{nullptr};
    comp::Dense<DataType> *dense{nullptr};
    if (compare_to_dense) {
        index::Index<uint64_t> I(0, n_rows), J(0, n_cols);
        data_block = new block::Dense<DataType>(n_rows, n_cols);
        dense = new comp::Dense<DataType>{data_block};
        full_block = new block::Block<DataType>{I, J, data_block, dense};
        full_block->compute(matrix);
    }

    uint64_t x_local_size{msk.local_col_size()}, y_local_size{msk.local_row_size()};
    uint64_t x_local_off{msk.local_col_off()}, y_local_off{msk.local_row_off()};

    SmartArray<DataType> x(n_cols), z(n_cols), y(n_rows), w(n_rows);
    DataType *x_local{x.data() + x_local_off}, *z_local{z.data() + x_local_off};
    DataType *y_local{y.data() + y_local_off}, *w_local{w.data() + y_local_off};

    DataType zero(0);
    std::uniform_real_distribution<double> unif(-1, 1);
    std::default_random_engine eng(703);

    double msk_matvec_time{0}, dense_matvec_time{0}, max_err{0};
    // Direct matvec check:
    // (true) w = A x
    // (appr) y = A x
    // check ||w - y|| / ||w||
    for (int s = 0; s < n_checks; ++s) {
        for (uint64_t i = 0; i < n_cols; ++i) {
            x[i] = unif(eng);
        }
        BLAS::rscal(n_cols, 1 / BLAS::nrm2(n_cols, x.data(), 1), x.data(), 1);

        BLAS::copy(y_local_size, &zero, 0, y_local, 1);
        auto t_start = std::chrono::steady_clock::now();
        msk.matvec(x_local, y_local, false);
        auto t_end = std::chrono::steady_clock::now();
        double t_diff = std::chrono::duration <double, std::milli> (t_end - t_start).count();
        msk_matvec_time += t_diff;

        if (compare_to_dense) {
            BLAS::copy(n_rows, &zero, 0, w.data(), 1);
            t_start = std::chrono::steady_clock::now();
            full_block->matvec(x.data(), w.data());
            t_end = std::chrono::steady_clock::now();
            t_diff = std::chrono::duration<double, std::milli>(t_end - t_start).count();
            dense_matvec_time += t_diff;

            BLAS::axpy(y_local_size, DataType{-1.0}, w_local, 1, y_local, 1);
            double y_local_err{BLAS::nrm2(y_local_size, y_local, 1)};
#ifdef MSK_USE_MPI
            y_local_err = y_local_err * y_local_err;
            if (proc_id == 0) {
                MPI_Reduce(MPI_IN_PLACE, &y_local_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            } else {
                MPI_Reduce(&y_local_err, &y_local_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
            y_local_err = std::sqrt(y_local_err);
#endif
            if (proc_id == 0) {
                double Ax_err = y_local_err / BLAS::nrm2(n_rows, w.data(), 1);
                max_err = std::max(max_err, Ax_err);
                if (log) log->debug("Matvec #{} err = {:.3e}", s, Ax_err);
            }
        }
    }
    msk_matvec_time /= n_checks;
    dense_matvec_time /= n_checks;
    if (proc_id == 0) {
        if (dense_matvec_time > 0) {
            if (log) log->info("Mean matvec time: msk = {:.4f}ms, dense = {:.4f}ms", msk_matvec_time, dense_matvec_time);
        } else {
            if (log) log->info("Mean matvec time: msk = {:.4f}ms", msk_matvec_time);
        }
        if (compare_to_dense) {
            if (log) log->info("Max matvec err = {:.3e}", max_err);
        }
    }
    
    
    msk_matvec_time = dense_matvec_time = max_err = 0;
    // Conjugated matvec check:
    // (true) z = A^* y
    // (appr) x = A^* y
    // check ||z - x|| / ||z||
    for (int s = 0; s < n_checks; ++s) {
        for (uint64_t i = 0; i < n_rows; ++i) {
            y[i] = unif(eng);
        }
        BLAS::rscal(n_rows, 1 / BLAS::nrm2(n_rows, y.data(), 1), y.data(), 1);

        BLAS::copy(x_local_size, &zero, 0, x_local, 1);
        auto t_start = std::chrono::steady_clock::now();
        msk.conj_matvec(y_local, x_local, false);
        auto t_end = std::chrono::steady_clock::now();
        double t_diff = std::chrono::duration <double, std::milli> (t_end - t_start).count();
        msk_matvec_time += t_diff;

        if (compare_to_dense) {
            BLAS::copy(n_cols, &zero, 0, z.data(), 1);
            t_start = std::chrono::steady_clock::now();
            full_block->conj_matvec(y.data(), z.data());
            t_end = std::chrono::steady_clock::now();
            t_diff = std::chrono::duration<double, std::milli>(t_end - t_start).count();
            dense_matvec_time += t_diff;

            BLAS::axpy(x_local_size, DataType{-1.0}, z_local, 1, x_local, 1);
            double x_local_err = BLAS::nrm2(x_local_size, x_local, 1);
#ifdef MSK_USE_MPI
            x_local_err = x_local_err * x_local_err;
            if (proc_id == 0) {
                MPI_Reduce(MPI_IN_PLACE, &x_local_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            } else {
                MPI_Reduce(&x_local_err, &x_local_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
            x_local_err = std::sqrt(x_local_err);
#endif

            if (proc_id == 0) {
                double ATy_err = x_local_err / BLAS::nrm2(n_cols, z.data(), 1);
                max_err = std::max(max_err, ATy_err);
                if (log) log->debug("Conjugated matvec #{} err = {:.3e}", s, ATy_err);
            }
        }
    }
    msk_matvec_time /= n_checks;
    dense_matvec_time /= n_checks;
    if (proc_id == 0) {
        if (dense_matvec_time > 0) {
            if (log) log->info("Mean conj_matvec time: msk = {:.4f}ms, dense = {:.4f}ms", msk_matvec_time, dense_matvec_time);
        } else {
            if (log) log->info("Mean conj_matvec time: msk = {:.4f}ms", msk_matvec_time);
        }
        if (compare_to_dense) {
            if (log) log->info("Max conj_matvec err = {:.3e}", max_err);
        }
    } 
    
    
    msk_matvec_time = dense_matvec_time = max_err = 0;
    // Transposed matvec check:
    // (true) z = A^T y
    // (appr) x = A^T y
    // check ||z - x|| / ||z||
    for (int s = 0; s < n_checks; ++s) {
        for (uint64_t i = 0; i < n_rows; ++i) {
            y[i] = unif(eng);
        }
        BLAS::rscal(n_rows, 1 / BLAS::nrm2(n_rows, y.data(), 1), y.data(), 1);

        BLAS::copy(x_local_size, &zero, 0, x_local, 1);
        auto t_start = std::chrono::steady_clock::now();
        msk.transp_matvec(y_local, x_local, false);
        auto t_end = std::chrono::steady_clock::now();
        double t_diff = std::chrono::duration <double, std::milli> (t_end - t_start).count();
        msk_matvec_time += t_diff;

        if (compare_to_dense) {
            BLAS::copy(n_cols, &zero, 0, z.data(), 1);
            t_start = std::chrono::steady_clock::now();
            full_block->transp_matvec(y.data(), z.data());
            t_end = std::chrono::steady_clock::now();
            t_diff = std::chrono::duration<double, std::milli>(t_end - t_start).count();
            dense_matvec_time += t_diff;

            BLAS::axpy(x_local_size, DataType{-1.0}, z_local, 1, x_local, 1);
            double x_local_err = BLAS::nrm2(x_local_size, x_local, 1);
#ifdef MSK_USE_MPI
            x_local_err = x_local_err * x_local_err;
            if (proc_id == 0) {
                MPI_Reduce(MPI_IN_PLACE, &x_local_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            } else {
                MPI_Reduce(&x_local_err, &x_local_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
            x_local_err = std::sqrt(x_local_err);
#endif

            if (proc_id == 0) {
                double ATy_err = x_local_err / BLAS::nrm2(n_cols, z.data(), 1);
                max_err = std::max(max_err, ATy_err);
                if (log) log->debug("Transposed matvec #{} err = {:.3e}", s, ATy_err);
            }
        }
    }
    msk_matvec_time /= n_checks;
    dense_matvec_time /= n_checks;
    if (proc_id == 0) {
        if (dense_matvec_time > 0) {
            if (log) log->info("Mean transp_matvec time: msk = {:.4f}ms, dense = {:.4f}ms", msk_matvec_time, dense_matvec_time);
        } else {
            if (log) log->info("Mean transp_matvec time: msk = {:.4f}ms", msk_matvec_time);
        }
        if (compare_to_dense) {
            if (log) log->info("Max transp_matvec err = {:.3e}", max_err);
        }
    }

    if (compare_to_dense) {
        delete full_block;
    }
}


/**
 * @brief Log mosaic-skeleton convolution performance time and error optionally.
 *
 * @details If to compare to dense matrix then it is acllocated and calculated.
 *
 * @tparam DataType Matrix element type.
 * @tparam Matrix Functional matrix type of MosaicSkeleton.
 * @param[in] msk Matrix in mosaic-skeleton format.
 * @param[in] log Logger.
 * @param[in] compare_to_dense If to calculate error in comperison to a dense matrix.
 * @param[in] n_checks How many runs to average.
 * @param[in] size Convolution size.
 */
template<class DataType, class Matrix, Trans trans = Trans::gen>
void check_convolution(MosaicSkeleton<DataType, trans> &msk, Matrix &matrix, std::shared_ptr<log::Logger> log,
                       bool compare_to_dense = true, int n_checks = 3, uint64_t size = 0) {
    if (n_checks < 1) {
        return;
    }

    uint64_t n_rows, n_cols;
    std::tie(n_rows, n_cols) = msk.size();

    // Calculate true dense matrix if required
    block::Block<DataType> *full_block{nullptr};
    block::Dense<DataType> *data_block{nullptr};
    comp::Dense<DataType> *dense{nullptr};
    if (compare_to_dense) {
        index::Index<uint64_t> I(0, n_rows), J(0, n_cols);
        data_block = new block::Dense<DataType>(n_rows, n_cols);
        dense = new comp::Dense<DataType>{data_block};
        full_block = new block::Block<DataType>{I, J, data_block, dense};
        full_block->compute(matrix);
    }

    if (size == 0 || size > n_rows + n_cols - 1) {
        size = n_rows + n_cols - 1;
    }
    uint64_t conv_size{size};
    SmartArray<DataType> y(n_rows), x(n_cols), z(conv_size), w(conv_size);

    constexpr DataType zero(0);
    std::uniform_real_distribution<double> unif(-1, 1);
    std::default_random_engine eng(703);

    double msk_conv_time{0}, dense_conv_time{0}, max_err{0};
    // Convolution check:
    // (true) w = conv(y, A, x)
    // (appr) z = conv(y, A_msk, x)
    // check ||w - z|| / ||w||
    for (int s = -1; s < n_checks; ++s) {
        for (uint64_t i = 0; i < n_rows; ++i) {
            y[i] = unif(eng);
        }
        BLAS::rscal(n_rows, 1 / BLAS::nrm2(n_rows, y.data(), 1), y.data(), 1);
        for (uint64_t i = 0; i < n_cols; ++i) {
            x[i] = unif(eng);
        }
        BLAS::rscal(n_cols, 1 / BLAS::nrm2(n_cols, x.data(), 1), x.data(), 1);

        auto t_start = std::chrono::steady_clock::now();
        msk.convolve(y.data(), x.data(), z.data(), conv_size);
        auto t_end = std::chrono::steady_clock::now();
        double t_diff = std::chrono::duration <double, std::milli> (t_end - t_start).count();
        if (s >= 0) {
            msk_conv_time += t_diff;
        }
        if (compare_to_dense) {
            BLAS::copy(conv_size, &zero, 0, w.data(), 1);
            t_start = std::chrono::steady_clock::now();
            full_block->convolve(y.data(), x.data(), w.data(), conv_size);
            t_end = std::chrono::steady_clock::now();
            t_diff = std::chrono::duration<double, std::milli>(t_end - t_start).count();
            if (s >= 0) {
                dense_conv_time += t_diff;
            }

            BLAS::axpy(conv_size, DataType{-1.0}, w.data(), 1, z.data(), 1);
            double conv_err = BLAS::nrm2(conv_size, z.data(), 1) / BLAS::nrm2(n_rows, w.data(), 1);
            max_err = std::max(max_err, conv_err);
            if (log) log->debug("Convolution #{} err = {:.3e}", s, conv_err);
        }
    }
    msk_conv_time /= n_checks;
    dense_conv_time /= n_checks;
    if (dense_conv_time > 0) {
        if (log) log->info("Mean convolution time: msk = {:.4f}ms, dense = {:.4f}ms", msk_conv_time, dense_conv_time);
    } else {
        if (log) log->info("Mean convolution time: msk = {:.4f}ms", msk_conv_time);
    }
    if (compare_to_dense) {
        if (log) log->info("Max convolution err = {:.3e}", max_err);
    }
    if (compare_to_dense) {
        delete full_block;
    }
}


/**
 * @brief Calculate singular value decomposition of a low-rank block column matrix.
 *
 * @details The procedure calculates the single orthogonal basis for rows.
 * Then it computes SVD for a small matrix (around sum of ranks) and recomputes factors.
 *
 * @tparam DataType Matrix element type.
 * @param[in] n_blocks Number of blocks in column.
 * @param[in] blocks List of low-rank blocks.
 * @param[out] U Left singular vectors.
 * @param[out] V Right singular vectors.
 * @param[out] S Singular values.
 * @param[in] tol SVD tolerance.
 * @param[in] frob If to perform SVD truncation by Frobenius norm or by spectral norm.
 * @return Resulting matrix rank.
 */
template<class DataType>
uint64_t block_svd(uint64_t n_blocks, const std::vector<const block::LowRank<DataType> *> &blocks,
                   SmartArray<DataType> &U, SmartArray<decltype(std::real(DataType{}))> &S, SmartArray<DataType> &V,
                   decltype(std::real(DataType{})) tol = 1e-12, bool frob = false) {
    using RealType = decltype(std::real(DataType()));

    const uint64_t n_common{blocks[0]->size().second};
    uint64_t rank_total{0}, col_basis_elements{0}, m_total{0};
    for (uint64_t i = 0; i < n_blocks; ++i) {
        uint64_t m, n, r;
        std::tie(m, n) = blocks[i]->size();
        r = blocks[i]->rank();
        rank_total += r;
        m_total += m;
        col_basis_elements += m * r;
        if (n != n_common) {
            return -1;
        }
    }
    // Extract row factors (V)
    SmartArray<DataType> U_total(col_basis_elements), tau_u(rank_total);
    SmartArray<DataType> V_total(n_common * rank_total);
    DataType *U_ptr{U_total.data()}, *v_ptr{V_total.data()};
    for (uint64_t i = 0; i < n_blocks; ++i) {
        uint64_t m, n, r;
        std::tie(m, n) = blocks[i]->size();
        r = blocks[i]->rank();
        blocks[i]->extract_col_factor(U_ptr, m);
        U_ptr += m * r;
        blocks[i]->extract_row_factor(v_ptr, n_common);
        v_ptr += n_common * r;
    }

    // Get orthogonal basis of rows
    uint64_t lwork = -1;
    int info{0};
    SmartArray<DataType> tau_v(rank_total), work(1);
    LAPACK::geqrf(n_common, rank_total, V_total.data(), n_common, tau_v.data(), work.data(), lwork, info);
    lwork = (uint64_t)std::real(work[0]) + 1;
    work.resize(lwork);
    LAPACK::geqrf(n_common, rank_total, V_total.data(), n_common, tau_v.data(), work.data(), lwork, info);

    // Extract R from V = Q R
    uint64_t min_nrt{std::min(n_common, rank_total)};
    SmartArray<DataType> R(min_nrt * rank_total);
    constexpr DataType zero{0};
    BLAS::copy(min_nrt * rank_total, &zero, 0, R.data(), 1);
    for (uint64_t j = 0; j < rank_total; ++j) {
        BLAS::copy(std::min(j + 1, n_common), V_total.data() + j * n_common, 1, R.data() + j * min_nrt, 1);
    }

    // Perform QR for every U_k if required
    U_ptr = U_total.data();
    DataType *r_ptr = R.data();
    DataType *tau_u_ptr{tau_u.data()};
    for (uint64_t i = 0; i < n_blocks; ++i) {
        uint64_t m, n, r;
        std::tie(m, n) = blocks[i]->size();
        r = blocks[i]->rank();
        if (blocks[i]->orth_cols()) { // Only scale by norms
            for (uint64_t j = 0; j < r; ++j) { // r <= m is guaranteed
                RealType norm = BLAS::nrm2(m, U_ptr + j * m, 1);
                if (std::abs(norm - 1.0) > 5 * std::numeric_limits<RealType>::epsilon()) {
                    BLAS::rscal(min_nrt, norm, r_ptr + j * min_nrt, 1);
                    BLAS::rscal(m, RealType{1.0} / norm, U_ptr + j * m, 1);
                }
            }
        } else { // True QR of U_k
            LAPACK::geqrf(m, r, U_ptr, m, tau_u_ptr, work.data(), lwork, info);
            BLAS::trmm('R', 'U', 'C', 'N', min_nrt, r, DataType{1.0}, U_ptr, m, r_ptr, min_nrt);
        }
        U_ptr += m * r;
        r_ptr += min_nrt * r;
        tau_u_ptr += r;
    }

    // Perform SVD over R = v s u^*
    SmartArray<DataType> V_tmp, U_tmp;
    auto result = perform_svd(min_nrt, rank_total, R, V_tmp, S, U_tmp, {0, tol, 0}, frob);
    uint64_t new_rank{result.memory / (min_nrt + rank_total)};

    // Compute final Vfinal = Q v
    V.resize(n_common * new_rank);
    for (uint64_t j = 0; j < new_rank; ++j) {
        BLAS::copy(min_nrt, V_tmp.data() + j * min_nrt, 1, V.data() + j * n_common, 1);
        if (n_common > min_nrt) {
            BLAS::copy(n_common - min_nrt, &zero, 0, V.data() + j * n_common + min_nrt, 1);
        }
    }
    LAPACK::unmqr('L', 'N', n_common, new_rank, std::min(n_common, rank_total), V_total.data(), n_common,
                  tau_v.data(), V.data(), n_common);

    // Compute final Ufinal = U u
    U.resize(m_total * new_rank);
    DataType *U_final_ptr{U.data()}, *U_tmp_ptr{U_tmp.data()};
    U_ptr = U_total.data();
    tau_u_ptr = tau_u.data();
    for (uint64_t i = 0; i < n_blocks; ++i) {
        uint64_t m, n, r;
        std::tie(m, n) = blocks[i]->size();
        r = blocks[i]->rank();
        if (blocks[i]->orth_cols()) {
            BLAS::gemm('N', 'N', m, new_rank, r, DataType{1.0}, U_ptr, m,
                       U_tmp_ptr, rank_total, DataType{0.0}, U_final_ptr, m_total);
        } else { // Applu reflections back
            for (uint64_t j = 0; j < new_rank; ++j) {
                BLAS::copy(r, U_tmp_ptr + j * rank_total, 1, U_final_ptr + j * m_total, 1);
                if (m > r) {
                    BLAS::copy(m - r, &zero, 0, U_final_ptr + j * m_total + r, 1);
                }
            }
            LAPACK::unmqr('L', 'N', m, new_rank, r, U_ptr, m, tau_u_ptr, U_final_ptr, m_total);
        }
        U_final_ptr += m;
        U_tmp_ptr += r;
        U_ptr += m * r;
        tau_u_ptr += r;
    }

    return new_rank;
}


extern template void run_gmres(const MosaicSkeleton<float> &msk_matrix,
               uint64_t n_rhs, float *rhs, uint64_t ld_rhs, float *sol, uint64_t ld_sol,
               solvers::linsys::Parameters<float> params, std::shared_ptr<log::Logger> log);

extern template void run_gmres(const MosaicSkeleton<double> &msk_matrix,
               uint64_t n_rhs, double *rhs, uint64_t ld_rhs, double *sol, uint64_t ld_sol,
               solvers::linsys::Parameters<double> params, std::shared_ptr<log::Logger> log);

extern template void run_gmres(const MosaicSkeleton<std::complex<float>> &msk_matrix,
               uint64_t n_rhs, std::complex<float> *rhs, uint64_t ld_rhs, std::complex<float> *sol, uint64_t ld_sol,
               solvers::linsys::Parameters<float> params, std::shared_ptr<log::Logger> log);

extern template void run_gmres(const MosaicSkeleton<std::complex<double>> &msk_matrix,
               uint64_t n_rhs, std::complex<double> *rhs, uint64_t ld_rhs, std::complex<double> *sol, uint64_t ld_sol,
               solvers::linsys::Parameters<double> params, std::shared_ptr<log::Logger> log);


extern template uint64_t block_svd(uint64_t n_blocks, const std::vector<const block::LowRank<float> *> &blocks,
                   SmartArray<float> &U, SmartArray<float> &S, SmartArray<float> &V,
                   float tol, bool frob);

extern template uint64_t block_svd(uint64_t n_blocks, const std::vector<const block::LowRank<double> *> &blocks,
                   SmartArray<double> &U, SmartArray<double> &S, SmartArray<double> &V,
                   double tol, bool frob);

extern template uint64_t block_svd(uint64_t n_blocks, const std::vector<const block::LowRank<std::complex<float>> *> &blocks,
                   SmartArray<std::complex<float>> &U, SmartArray<float> &S, SmartArray<std::complex<float>> &V,
                   float tol, bool frob);

extern template uint64_t block_svd(uint64_t n_blocks, const std::vector<const block::LowRank<std::complex<double>> *> &blocks,
                   SmartArray<std::complex<double>> &U, SmartArray<double> &S, SmartArray<std::complex<double>> &V,
                   double tol, bool frob);


}

#endif //MOSAIC_SKELETON_UTILITY_TESTING_HPP
