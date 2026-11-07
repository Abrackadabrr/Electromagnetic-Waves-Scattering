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
 * @brief Defines C++ language side of C/C++ interface to class MosaicSkeleton.
 */


#ifndef MOSAIC_SKELETON_INTERFACE_C_MSK_HPP
#define MOSAIC_SKELETON_INTERFACE_C_MSK_HPP


#include "msk_config.h"

#ifdef MSK_USE_MPI
#include "cpp_mpi/MPI.hpp"
#endif

#include <cmath>
#include <cstdint>
#include <string>
#include <fstream>
#include <functional>
#include <exception>

#include "msk/matrix/basic.hpp"
#include "msk/geometry/point.hpp"
#include "msk/oracle/distant.hpp"
#include "msk/mosaic_skeleton.hpp"
#include "msk/utility/logging.hpp"
#include "msk/utility/testing.hpp"


#define MSK_C_RESULT_OK 0
#define MSK_C_RESULT_ERR 1
#define MSK_C_RESULT_ERR_RUNTIME 2
#define MSK_C_RESULT_ERR_MEM 3

namespace MSk {

namespace C {

enum Result : int {
    OK = MSK_C_RESULT_OK,
    Err = MSK_C_RESULT_ERR,
    ErrRT = MSK_C_RESULT_ERR_RUNTIME,
    ErrMem = MSK_C_RESULT_ERR_MEM,
};

/**
 * @brief Basic functional matrix with arbitrary geometry.
 *
 * @tparam DataType Matrix element type.
 * @tparam PointType Points coordinates type.
 * @tparam dim Points spatial dimension.
 */
template<class DataType, class PointType = decltype(std::abs(DataType{})), uint64_t dim = 3>
class RawCenter : public matrix::Matrix<DataType> {
public:
    using obj_type = geom::Point<PointType, dim>;

private:
    std::function<DataType(uint64_t, uint64_t)> _function; ///< Element evaluation function
    geom::PointGrid<PointType, dim> rows, cols; ///< Rows and columns related geometry

public:
    /**
     * @brief Initialize element evaluation rule and geometry.
     *
     * @param[in] element_function Function to evaluate matrix elements.
     * @param[in, out] row_pts Rows related points, moved out.
     * @param[in, out] col_pts Columns related points, moved out.
     */
    RawCenter(std::function<DataType(uint64_t, uint64_t)> element_function,
              geom::PointGrid<PointType, dim> &&row_pts, geom::PointGrid<PointType, dim> &&col_pts) :
              matrix::Matrix<DataType>{row_pts.size(), col_pts.size()},
              _function{element_function}, rows{row_pts}, cols{col_pts} {}

    DataType element(uint64_t i, uint64_t j) const override {
        return _function(i, j);
    }

    /**
     * @brief Get point corresponding to a row.
     *
     * @param[in] i Row index.
     * @return Point.
     */
    const obj_type &row_point(uint64_t i) const {
        return rows[i];
    }

    /**
     * @brief Get point corresponding to a column.
     *
     * @param[in] j Column index.
     * @return Point.
     */
    const obj_type &col_point(uint64_t j) const {
        return cols[j];
    }
};


template<class DataType, class PointType = decltype(std::abs(DataType{})), uint64_t dim = 3>
Result build_basic(MosaicSkeleton<DataType> **msk_ptr, 
                   uint64_t m, uint64_t n, std::function<DataType(uint64_t, uint64_t)> f,
                   const PointType *row_pts, const PointType *col_pts,
                   oracle::Parameters params, uint64_t n_jobs = omp_get_max_threads()
#ifdef MSK_USE_MPI
                   , MPI_Comm comm = MPI_COMM_WORLD
#endif
)
{
    // Build logger
    auto main_logger = params.log->clone("MSk");

    // Build geometry
    using RealType = decltype(std::abs(DataType(1.0)));
    geom::PointGrid<RealType, dim> rows, cols;
    try {
        rows = geom::PointGrid<RealType, dim>(m, reinterpret_cast<const geom::Point<PointType, dim> *>(row_pts));
        cols = geom::PointGrid<RealType, dim>(n, reinterpret_cast<const geom::Point<PointType, dim> *>(col_pts));
    } catch (const std::exception &ex) {
        std::cerr << "Error while processing geometry: " << ex.what() << std::endl;
        return Result::ErrMem;
    }

    // Build matrix
    using MatrixType = RawCenter<DataType, PointType, dim>;
    MatrixType matrix{f, std::move(rows), std::move(cols)};

    // Build oracle
    using OracleType = oracle::Distant<DataType, MatrixType>;
    OracleType oracle{matrix, params};

    // Do approximation
    MosaicSkeleton<DataType> *msk{nullptr};
    try {
        msk = new MosaicSkeleton<DataType>{};
    } catch (const std::exception &ex) {
        std::cerr << "Error while allocating MSk matrx: " << ex.what() << std::endl;
        return Result::ErrMem;
    }
    auto t_start = std::chrono::steady_clock::now();
    try {
#ifdef MSK_USE_MPI
        *msk = std::move(MosaicSkeleton<DataType>::approximate(matrix, oracle, n_jobs, comm, main_logger));
#else
        *msk = std::move(MosaicSkeleton<DataType>::approximate(matrix, oracle, n_jobs, main_logger));
#endif
    } catch (const std::exception &ex) {
        std::cerr << "Error while building MSk matrx: " << ex.what() << std::endl;
        return Result::ErrRT;
    }
    auto t_end = std::chrono::steady_clock::now();
    double t_approx = std::chrono::duration<double>(t_end - t_start).count();

    int proc_id{0};
#ifdef MSK_USE_MPI
    MPI_Comm_rank(comm, &proc_id);
#endif
    if (main_logger && proc_id == 0) main_logger->debug("Approximation time = {:.3e}s", t_approx);

    *msk_ptr = msk;
    return Result::OK;
}

template<class DataType, class PointType = decltype(std::abs(DataType{})), uint64_t dim = 3>
Result build_block(MosaicSkeleton<DataType> **msk_ptr, 
                   uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                   uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                   std::function<DataType(uint64_t, uint64_t)> f,
                   const PointType *row_pts, const PointType *col_pts,
                   oracle::Parameters params, uint64_t n_jobs = omp_get_max_threads()
#ifdef MSK_USE_MPI
                   , MPI_Comm comm = MPI_COMM_WORLD
#endif
) {
    // Basic building in case of no real block structure
    if (m_blocks <= 1 && n_blocks <= 1) {
#ifdef MSK_USE_MPI
        return build_basic(msk_ptr, m, n, f, row_pts, col_pts, params, n_jobs, comm);
#else
        return build_basic(msk_ptr, m, n, f, row_pts, col_pts, params, n_jobs);
#endif
    }

    // Build logger
    auto main_logger = params.log->clone("MSk");

    // Build geometry
    using RealType = decltype(std::abs(DataType(1.0)));
    geom::PointGrid<RealType, dim> rows, cols;
    try {
        rows = geom::PointGrid<RealType, dim>(m, reinterpret_cast<const geom::Point<PointType, dim> *>(row_pts));
        cols = geom::PointGrid<RealType, dim>(n, reinterpret_cast<const geom::Point<PointType, dim> *>(col_pts));
    } catch (const std::exception &ex) {
        std::cerr << "Error while processing geometry: " << ex.what() << std::endl;
        return Result::ErrMem;
    }

    // Build matrix
    using MatrixType = RawCenter<DataType, PointType, dim>;
    MatrixType matrix{f, std::move(rows), std::move(cols)};

    // Build oracle
    using OracleType = oracle::Distant<DataType, MatrixType>;
    OracleType oracle{matrix, params};
    using TBlock = block::Block<DataType>;
    using TIdx = index::Index<uint64_t>;
    TIdx idx_rows(0, m), idx_cols(0, n);
    TBlock base(idx_rows, idx_cols);

    decltype(oracle.get_action(&base)) answer{};
    answer.n_child = m_blocks * n_blocks;
    answer.child.resize(answer.n_child);
    idx_rows = {0, m};
    for (uint64_t i{0}; i < m_blocks; ++i) {
        TIdx idx_rows_left, idx_rows_right;
        std::tie(idx_rows_left, idx_rows_right) = idx_rows.split(idx_rows.begin() + m_block_sizes[i]);
        idx_rows = idx_rows_right;
        idx_cols = {0, n};
        for (uint64_t j{0}; j < n_blocks; ++j) {
            TIdx idx_cols_up, idx_cols_down;
            std::tie(idx_cols_up, idx_cols_down) = idx_cols.split(idx_cols.begin() + n_block_sizes[j]);
            idx_cols = idx_cols_down;
            answer.child[i * n_blocks + j] = std::move(TBlock(idx_rows_left, idx_cols_up));
        }
    }
    oracle.add_manual(std::move(base), std::move(answer));

    // Do approximation
    MosaicSkeleton<DataType> *msk{nullptr};
    try {
        msk = new MosaicSkeleton<DataType>{};
    } catch (const std::exception &ex) {
        std::cerr << "Error while allocating MSk matrx: " << ex.what() << std::endl;
        return Result::ErrMem;
    }
    auto t_start = std::chrono::steady_clock::now();
    try {
#ifdef MSK_USE_MPI
        *msk = std::move(MosaicSkeleton<DataType>::approximate(matrix, oracle, n_jobs, comm, main_logger));
#else
        *msk = std::move(MosaicSkeleton<DataType>::approximate(matrix, oracle, n_jobs, main_logger));
#endif
    } catch (const std::exception &ex) {
        std::cerr << "Error while building block MSk matrx: " << ex.what() << std::endl;
        return Result::ErrRT;
    }
    auto t_end = std::chrono::steady_clock::now();
    double t_approx = std::chrono::duration<double>(t_end - t_start).count();

    int proc_id{0};
#ifdef MSK_USE_MPI
    MPI_Comm_rank(comm, &proc_id);
#endif
    if (main_logger && proc_id == 0) main_logger->debug("Approximation time = {:.3e}s", t_approx);

    *msk_ptr = msk;
    return Result::OK;
}

template<class DataType>
Result recompress(void *msk_matrix, double rel_tol) {
    using MSkType = MosaicSkeleton<DataType>;

    auto msk_ptr = reinterpret_cast<MSkType *>(msk_matrix);

    // Do recompression
    auto t_start = std::chrono::steady_clock::now();
    try {
        msk_ptr->recompress(rel_tol);
    } catch (const std::exception &ex) {
        std::cerr << "Error while recompressing MSk matrx: " << ex.what() << std::endl;
        return Result::ErrRT;
    }
    auto t_end = std::chrono::steady_clock::now();
    double t_recomp = std::chrono::duration<double>(t_end - t_start).count();
    int proc_id{0};
#ifdef MSK_USE_MPI
    MPI_Comm_rank(msk_ptr->comm(), &proc_id);
#endif
    if (proc_id == 0) {
        msk_ptr->logger()->info("Recompression time = {:.3e}s", t_recomp);
    }
    return Result::OK;
}

template<class DataType>
Result destroy(void *msk_matrix) {
    using MSkType = MosaicSkeleton<DataType>;

    auto msk_ptr = reinterpret_cast<MSkType *>(msk_matrix);

    int proc_id{0};
#ifdef MSK_USE_MPI
    MPI_Comm_rank(msk_ptr->comm(), &proc_id);
#endif

    auto main_logger = msk_ptr->logger();

    try {
        delete msk_ptr;
    } catch (const std::exception &ex) {
        std::cerr << "Error while deallocating MSk matrix: " << ex.what() << std::endl;
    }

    if (proc_id == 0) {
        main_logger->info("Destroyed MSk");
    }
    return Result::OK;
}

template<class DataType>
Result matvec(void *msk_matrix, const char trans, const DataType *x, DataType *y) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->matvec(trans, x, y);
    } catch (const std::exception &ex) {
        std::cerr << "Error while performing MSk matvec: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
#ifdef MSK_USE_MPI
Result to_dense(void *msk_matrix, DataType *ptr, uint64_t ldim, int root = 0)
#else
Result to_dense(void *msk_matrix, DataType *ptr, uint64_t ldim)
#endif
{
    using MSkType = MosaicSkeleton<DataType>;

    auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
    try {
#ifdef MSK_USE_MPI
        msk_ptr->to_dense(ptr, ldim, root);
#else
        msk_ptr->to_dense(ptr, ldim);
#endif
    } catch (const std::exception &ex) {
        std::cerr << "Error while densifying MSk matrx: " << ex.what() << std::endl;
    }
    return Result::OK;
}

#ifdef MSK_USE_MPI
template<class DataType>
Result permute_col(void *msk_matrix, const DataType *x_in, DataType *x_out, bool global = true) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->permute_col(x_in, x_out, global);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-permuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result permute_col(void *msk_matrix, DataType *x, bool global = true) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->permute_col(x, global);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-permuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result unpermute_col(void *msk_matrix, const DataType *x_in, DataType *x_out, bool global = true) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->unpermute_col(x_in, x_out, global);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-unpermuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result unpermute_col(void *msk_matrix, DataType *x, bool global = true) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->unpermute_col(x, global);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-unpermuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result permute_row(void *msk_matrix, const DataType *y_in, DataType *y_out, bool global = true) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->permute_row(y_in, y_out, global);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-permuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result permute_row(void *msk_matrix, DataType *y, bool global = true) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->permute_row(y, global);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-permuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result unpermute_row(void *msk_matrix, const DataType *y_in, DataType *y_out, bool global = true) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->unpermute_row(y_in, y_out, global);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-unpermuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result unpermute_row(void *msk_matrix, DataType *y, bool global = true) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->unpermute_row(y, global);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-unpermuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}
#else
template<class DataType>
Result permute_col(void *msk_matrix, const DataType *x_in, DataType *x_out) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        const auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->permute_col(x_in, x_out);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-permuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result permute_col(void *msk_matrix, DataType *x) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        const auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->permute_col(x);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-permuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result unpermute_col(void *msk_matrix, const DataType *x_in, DataType *x_out) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        const auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->unpermute_col(x_in, x_out);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-unpermuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result unpermute_col(void *msk_matrix, DataType *x) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        const auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->unpermute_col(x);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-unpermuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result permute_row(void *msk_matrix, const DataType *y_in, DataType *y_out) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        const auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->permute_row(y_in, y_out);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-permuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result permute_row(void *msk_matrix, DataType *y) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        const auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->permute_row(y);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-permuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result unpermute_row(void *msk_matrix, const DataType *y_in, DataType *y_out) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        const auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->unpermute_row(y_in, y_out);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-unpermuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}

template<class DataType>
Result unpermute_row(void *msk_matrix, DataType *y) {
    using MSkType = MosaicSkeleton<DataType>;

    try {
        const auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
        msk_ptr->unpermute_row(y);
    } catch (const std::exception &ex) {
        std::cerr << "Error while MSk-unpermuting vector: " << ex.what() << std::endl;
    }
    return Result::OK;
}
#endif

template<class DataType>
#ifdef MSK_USE_MPI
int check_error(void *msk_matrix, std::function<DataType(uint64_t, uint64_t)> f, bool relative, bool global, double *error) {
#else
int check_error(void *msk_matrix, std::function<DataType(uint64_t, uint64_t)> f, bool relative, double *error) {
#endif
    using MSkType = MosaicSkeleton<DataType>;
    using RefMatrixType = matrix::Elementary<DataType>;

    auto msk_ptr = reinterpret_cast<const MSkType *>(msk_matrix);
    uint64_t m{}, n{};
    std::tie(m, n) = msk_ptr->size();
    auto ref_matrix = RefMatrixType(m, n, f);
    try {
#ifdef MSK_USE_MPI
        *error = msk_ptr->check_error(ref_matrix, relative, global);
#else
        *error = msk_ptr->check_error(ref_matrix, relative);
#endif
    } catch (const std::exception &ex) {
        std::cerr << "Error while cimpiting MSk approximation error: " << ex.what() << std::endl;
        return Result::ErrRT;
    }
    return Result::OK;
}


template<class DataType, class PointType = decltype(std::abs(DataType{})), uint64_t dim = 3>
Result approx_block_solve(
        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
        std::function<DataType(uint64_t, uint64_t)> f,
        const PointType *row_pts, const PointType *col_pts,
        oracle::Parameters params,
        uint64_t n_rhs, DataType *B, DataType *X,
        solvers::linsys::Parameters<decltype(std::abs(DataType{}))> solve_params,
        uint64_t n_jobs = omp_get_max_threads()
#ifdef MSK_USE_MPI
        , MPI_Comm comm = MPI_COMM_WORLD
#endif
) {
    MosaicSkeleton<DataType> *msk_ptr;
#ifdef MSK_USE_MPI
    Result result = build_block(&msk_ptr, m, m_blocks, m_block_sizes, n, n_blocks, n_block_sizes,
                                f, row_pts, col_pts, params, n_jobs, comm);
#else
    Result result = build_block(&msk_ptr, m, m_blocks, m_block_sizes, n, n_blocks, n_block_sizes,
                                f, row_pts, col_pts, params, n_jobs);
#endif
    if (result != Result::OK) {
        return result;
    }
    try {
        msk_ptr->recompress(params.block_comp.rel_tol);
    } catch (const std::exception &ex) {
        std::cerr << "Error while recompressing MSk matrx: " << ex.what() << std::endl;
        return Result::ErrRT;
    }
    try {
        run_gmres(*msk_ptr, n_rhs,
                  B, msk_ptr->local_row_size(),
                  X, msk_ptr->local_col_size(),
                  solve_params, msk_ptr->logger());
    } catch (const std::exception &ex) {
        std::cerr << "Error while solving a system with MSk matrx using GMRES: " << ex.what() << std::endl;
        return Result::ErrRT;
    }
    try {
        delete msk_ptr;
    } catch (const std::exception &ex) {
        std::cerr << "Error while deallocating MSk matrix: " << ex.what() << std::endl;
        return Result::ErrMem;
    }
    return Result::OK;
}

} // namespace C
} // namespace MSk

#endif // MOSAIC_SKELETON_INTERFACE_C_MSK_HPP
