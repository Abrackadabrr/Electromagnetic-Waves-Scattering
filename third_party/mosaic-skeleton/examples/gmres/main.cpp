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
// Created by b.valiakhmetov on 27.04.23.
//


#ifdef MSK_EXAMPLES_USE_MPI
#include <mpi.h>
#endif

#include <random>

#include "toml.hpp"
#include "common.hpp"
#include "smart_array/smart_array.hpp"

#include "msk/mosaic_skeleton.hpp"
#include "msk/geometry/point_grid.hpp"
#include "msk/matrix/basic.hpp"
#include "msk/oracle/distant.hpp"
#include "msk/utility/testing.hpp"
#include "msk/solvers/linsys/GMRES.hpp"
#include "msk/utility/logging.hpp"


constexpr uint64_t PROBLEM_DIMENSION = 2;


int main(int argc, char *argv[]) {
    int proc_id{0};
#ifdef MSK_EXAMPLES_USE_MPI
    const MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_thread_req{MPI_THREAD_MULTIPLE}, mpi_thread_prov{0};
    MPI_Init_thread(&argc, &argv, mpi_thread_req, &mpi_thread_prov);
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);

    assert(mpi_thread_prov >= mpi_thread_req);
    MPI_Comm_rank(comm, &proc_id);
#endif

    using DataType = std::complex<double>;
    using RealType = decltype(std::abs(DataType(1.0)));

    // Parse problem parameters
    std::string problem_file;
    if (argc > 1) {
        problem_file = argv[1];
    } else {
        problem_file = std::string(MSK_EXAMPLES_DIRECTORY) + "/gmres/parameters.toml";
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
        MSk::geom::generate_line(pts, n_pts);
    } else if (problem == "circle") {
        n_pts = N;
        pts = new MSk::geom::Point<RealType, PROBLEM_DIMENSION>[n_pts];
        MSk::geom::generate_circle(pts, n_pts, 1.0);
    } else if (problem == "square") {
        n_pts = N * N;
        pts = new MSk::geom::Point<RealType, PROBLEM_DIMENSION>[n_pts];
        MSk::geom::generate_cube(pts, N, PROBLEM_DIMENSION, {0.0, 0.0}, {1.0, 1.0});
    }
    MSk::geom::PointGrid<RealType, PROBLEM_DIMENSION> geometry(n_pts, pts);

    // Build matrix
    using MatrixType = MSk::matrix::BiPoint<DataType, RealType, 2>;
    MatrixType matr(std::bind(MSk::helmgoltz_fundamental<RealType, 2>, std::placeholders::_1, std::placeholders::_2,
                              DataType(1.0)), geometry, geometry);

    // Logger
    auto &msk_params = toml::find(params, "msk");
    auto main_logger = MSk::build_logger(msk_params["log"], "MSk", "stdout", "info");

    // Oracle
    auto &oracle_params = toml::find(params, "oracle");
    auto &approx_params = toml::find(params, "approximation");
    using OracleType = MSk::oracle::Distant<DataType, MatrixType>;
    OracleType oracle{matr, MSk::oracle::read_parameters(oracle_params, approx_params)};

    double rel_tol = toml::find_or<double>(approx_params, "rel_tol", 1e-3);

    // Do approximation and recompression
#ifdef MSK_EXAMPLES_USE_MPI
    auto msk = MSk::MosaicSkeleton<DataType>::approximate(matr, oracle, omp_get_max_threads(), comm, main_logger);
#else
    auto msk = MSk::MosaicSkeleton<DataType>::approximate(matr, oracle, omp_get_max_threads(), main_logger);
#endif
    msk.recompress(rel_tol);

    // Check approximation
    auto &check_params = toml::find(params, "check");
    if (check_params.contains("compression") && check_params["compression"].as_boolean()) {
        MSk::check_compression<DataType>(msk, main_logger);
    }
    if (check_params.contains("error") && check_params["error"].as_boolean()) {
        MSk::check_error<DataType, MatrixType>(msk, matr, main_logger);
    }
    if (check_params.contains("matvec") && check_params["matvec"].as_boolean()) {
        MSk::check_matvec<DataType, MatrixType>(msk, matr, main_logger);
    }

    // Generate random solution X and corresponding rhs B = AX
    uint64_t n_rhs{std::min(n_pts, 10ul)};

    uint64_t local_row_size{msk.local_row_size()}, local_col_size{msk.local_col_size()};
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    std::default_random_engine eng(42 * (proc_id + 1));
    SmartArray<DataType> B(local_row_size * n_rhs), B_copy(local_row_size * n_rhs);
    SmartArray<DataType> X(local_col_size * n_rhs), X_true(local_col_size * n_rhs);
    constexpr DataType zero{0};
    BLAS::copy(local_row_size * n_rhs, &zero, 0, B.data(), 1);
    for (uint64_t j = 0; j < n_rhs; ++j) {
        for (uint64_t i = 0; i < local_col_size; ++i) {
            X[i + local_col_size * j] = unif(eng) + DataType(0.0, -1.0) * unif(eng);
        }
        RealType nrm{1.0 / BLAS::nrm2(local_col_size, &X[local_col_size * j], 1)};
        BLAS::rscal(local_col_size, nrm, &X[local_col_size * j], 1);
        msk.matvec(X.data() + local_col_size * j, B.data() + local_row_size * j, false);
    }
    BLAS::copy(local_col_size * n_rhs, X.data(), 1, X_true.data(), 1);
    BLAS::copy(local_row_size * n_rhs, B.data(), 1, B_copy.data(), 1);
    BLAS::copy(local_col_size * n_rhs, &zero, 0, X.data(), 1);

    // Get GMRES parameters
    MSk::solvers::linsys::Parameters<RealType> solve_params;
    MSk::set_gmres_parameters<RealType>(params["solve"], solve_params);

    // Solve linear system AX = B using GMRES
    MSk::run_gmres(msk, n_rhs, B.data(), local_row_size, X.data(), local_col_size, solve_params);

    // Check maximal error and residual
    SmartArray<DataType> Ax(local_row_size), x(local_col_size);
    double max_res{0}, max_err{0};
    for (uint64_t j = 0; j < n_rhs; ++j) {
        BLAS::copy(local_col_size, X.data() + local_col_size * j, 1, x.data(), 1);
        BLAS::axpy(local_col_size, DataType(-1.0), X_true.data() + local_col_size * j, 1, x.data(), 1);
        double err = BLAS::nrm2(local_col_size, x.data(), 1);
        double nrm = BLAS::nrm2(local_col_size, X_true.data() + local_col_size * j, 1);
#ifdef MSK_EXAMPLES_USE_MPI
        err = err * err;
        nrm = nrm * nrm;
        MPI::Allreduce_inplace_sum(&err, 1, comm);
        MPI::Allreduce_inplace_sum(&nrm, 1, comm);
        err = std::sqrt(err);
        nrm = std::sqrt(nrm);
#endif
        err /= nrm;
        max_err = std::max(err, max_err);

        BLAS::copy(local_row_size, &zero, 0, Ax.data(), 1);
        msk.matvec(X.data() + local_col_size * j, Ax.data(), false);
        BLAS::axpy(local_row_size, DataType(-1.0), B_copy.data() + local_row_size * j, 1, Ax.data(), 1);
        double res = BLAS::nrm2(local_row_size, Ax.data(), 1);
        nrm = BLAS::nrm2(local_row_size, B_copy.data() + local_row_size * j, 1);

#ifdef MSK_EXAMPLES_USE_MPI
        res = res * res;
        nrm = nrm * nrm;
        MPI::Allreduce_inplace_sum(&res, 1, comm);
        MPI::Allreduce_inplace_sum(&nrm, 1, comm);
        res = std::sqrt(res);
        nrm = std::sqrt(nrm);
#endif
        res /= nrm;
        max_res = std::max(res, max_res);
    }
    if (proc_id == 0) {
        main_logger->info("Max error: {:.4e}", max_err);
        main_logger->info("Max residual: {:.4e}", max_res);
    }
    delete[] pts;

#ifdef MSK_EXAMPLES_USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
