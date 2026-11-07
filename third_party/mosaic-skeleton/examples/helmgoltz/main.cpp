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
// Created by bulat_v on 04.02.2022.
//


#include <cstdint>
#include <sys/types.h>
#ifdef MSK_EXAMPLES_USE_MPI
#include <mpi.h>
#endif

#include <random>
#include <complex>
#include <functional>

#include "toml.hpp"
#include "common.hpp"

#include "msk/matrix/basic.hpp"
#include "msk/geometry/point_grid.hpp"
#include "msk/oracle/distant.hpp"
#include "msk/mosaic_skeleton.hpp"
#include "msk/utility/testing.hpp"
#include "msk/utility/logging.hpp"


int main(int argc, char *argv[]) {
    constexpr uint64_t PROBLEM_DIMENSION = 2;

#ifdef MSK_EXAMPLES_USE_MPI
    int mpi_thread_req{MPI_THREAD_MULTIPLE}, mpi_thread_prov{0};
    MPI_Init_thread(&argc, &argv, mpi_thread_req, &mpi_thread_prov);
    MPI_Comm comm{MPI_COMM_WORLD};
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);

    assert(mpi_thread_prov >= mpi_thread_req);
#endif

    int proc_id{0};
#ifdef MSK_EXAMPLES_USE_MPI
    MPI_Comm_rank(comm, &proc_id);
#endif

    using DataType = std::complex<double>;
    using RealType = decltype(std::abs(DataType{}));

    // Parse problem parameters
    std::string problem_file;
    if (argc > 1) {
        problem_file = argv[1];
    } else {
        problem_file = std::string(MSK_EXAMPLES_DIRECTORY) + "/helmgoltz/parameters.toml";
    }
    auto params = toml::parse(problem_file);

    // Build geometry
    uint64_t n_pts{0};
    MSk::geom::Point<RealType, 2> *pts{nullptr};

    uint64_t N = toml::find_or(params, "geometry", "side", 100);
    std::string problem = toml::find_or(params, "geometry", "type", std::string{"line"});
    if (problem == "line") {
        n_pts = N;
        pts = new MSk::geom::Point<RealType, PROBLEM_DIMENSION>[n_pts];
        MSk::geom::generate_line<RealType>(pts, n_pts);
    } else if (problem == "circle") {
        n_pts = N;
        pts = new MSk::geom::Point<RealType, PROBLEM_DIMENSION>[n_pts];
        MSk::geom::generate_circle<RealType>(pts, n_pts, 1.0);
    } else if (problem == "square") {
        n_pts = N * N;
        pts = new MSk::geom::Point<RealType, PROBLEM_DIMENSION>[n_pts];
        MSk::geom::generate_cube<RealType>(pts, N, PROBLEM_DIMENSION, {0.0, 0.0}, {1.0, 1.0});
    }
    MSk::geom::PointGrid<RealType, PROBLEM_DIMENSION> geometry(n_pts, pts);

    // Build functional matrix
    using MatrixType = MSk::matrix::BiPoint<DataType, RealType, 2>;
    MatrixType matr(std::bind(MSk::helmgoltz_fundamental<RealType, 2>,
                    std::placeholders::_1, std::placeholders::_2,
                    DataType{100.0}), geometry, geometry);

    auto &msk_params = toml::find(params, "msk");
    auto main_logger = MSk::build_logger(msk_params["log"], "MSk", "stdout", "info");

    // Build oracle
    auto &oracle_params = toml::find(params, "oracle");
    auto &approx_params = toml::find(params, "approximation");
    using OracleType = MSk::oracle::Distant<DataType, MatrixType>;
    OracleType oracle{matr, MSk::oracle::read_parameters(oracle_params, approx_params)};
    oracle.symmetry(true);

    constexpr MSk::Trans TRANS = MSk::Trans::symm;
    // Do approximation
    uint64_t n_jobs = toml::find<toml::integer>(params, "msk", "n_jobs");
    auto t_start = std::chrono::steady_clock::now();
#ifdef MSK_EXAMPLES_USE_MPI
    auto msk = MSk::MosaicSkeleton<DataType, TRANS>::approximate(matr, oracle, n_jobs, comm, main_logger);
#else
    auto msk = MSk::MosaicSkeleton<DataType, TRANS>::approximate(matr, oracle, n_jobs, main_logger);
#endif
    auto t_end = std::chrono::steady_clock::now();
    double t_diff = std::chrono::duration <double> (t_end - t_start).count();
    if (proc_id == 0) {
        main_logger->info("Approximation time = {:.3e}s", t_diff);
    }

    // Do recompression
    double rel_tol = toml::find_or<double>(approx_params, "rel_tol", 1e-3);
    t_start = std::chrono::steady_clock::now();
    msk.recompress(rel_tol);
    t_end = std::chrono::steady_clock::now();
    t_diff = std::chrono::duration<double>(t_end - t_start).count();
    if (proc_id == 0) {
        main_logger->info("Recompression time = {:.3e}s", t_diff);
    }

    // Do checks
    auto &check_params = toml::find(params, "check");
    if (check_params.contains("compression") && check_params["compression"].as_boolean()) {
        MSk::check_compression<DataType, TRANS>(msk, main_logger);
    }
    if (check_params.contains("error") && check_params["error"].as_boolean()) {
        MSk::check_error<DataType, MatrixType, TRANS>(msk, matr, main_logger);
    }
    if (check_params.contains("matvec") && check_params["matvec"].as_boolean()) {
        MSk::check_matvec<DataType, MatrixType, TRANS>(msk, matr, main_logger, true, 1);
    }

    constexpr uint64_t block_size = 100;
    constexpr bool pre_permuted{false};
    uint64_t local_row_size{msk.local_row_size()}, local_col_size{msk.local_col_size()};
    SmartArray<DataType> X_block{local_col_size * block_size}, Y_block{local_row_size * block_size}, Y_block_true{local_row_size * block_size};

    std::uniform_real_distribution<double> unif(-1, 1);
    std::default_random_engine eng(703);
    for (uint64_t i = 0; i < X_block.size(); ++i) {
        X_block[i] = unif(eng);
    }

    t_start = std::chrono::steady_clock::now();
    for (uint64_t k = 0; k < block_size; ++k) {
        msk.matvec(X_block.data() + k * local_col_size, Y_block_true.data() + k * local_row_size, pre_permuted);
    }
    t_end = std::chrono::steady_clock::now();
    t_diff = std::chrono::duration<double>(t_end - t_start).count();
    if (proc_id == 0) {
        main_logger->info("Separate matvec time = {:.3e}s", t_diff);
    }

    t_start = std::chrono::steady_clock::now();
    msk.matvec(block_size, X_block.data(), Y_block.data(), pre_permuted);
    t_end = std::chrono::steady_clock::now();
    t_diff = std::chrono::duration<double>(t_end - t_start).count();
    if (proc_id == 0) {
        main_logger->info("Block matvec time = {:.3e}s", t_diff);
    }

    BLAS::axpy(local_row_size * block_size, DataType{-1}, Y_block_true.data(), 1, Y_block.data(), 1);
    
    RealType nrm = BLAS::nrm2(local_row_size * block_size, Y_block_true.data(), 1);
    RealType err = BLAS::nrm2(local_row_size * block_size, Y_block.data(), 1);
    
    main_logger->debug("{} :: Block matvec err = {:.3e}", proc_id, err / nrm);

#ifdef MSK_USE_MPI
    err = err * err;
    nrm = nrm * nrm;
    MPI::Allreduce_inplace_sum(&err, 1, comm);
    MPI::Allreduce_inplace_sum(&nrm, 1, comm);
    err = std::sqrt(err);
    nrm = std::sqrt(nrm);
#endif
    if (proc_id == 0) {
        main_logger->info("Block matvec err = {:.3e}", err / nrm);
    }

    delete[] pts;
    

#ifdef MSK_EXAMPLES_USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
