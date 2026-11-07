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


#ifdef MSK_EXAMPLES_USE_MPI
#include <mpi.h>
#endif

#include <memory>
#include <random>
#include <chrono>

#include "toml.hpp"
#include "common.hpp"

#include "msk/matrix/basic.hpp"
#include "msk/mosaic_skeleton.hpp"
#include "msk/oracle/distant.hpp"
#include "msk/utility/testing.hpp"
#include "msk/utility/logging.hpp"


constexpr uint64_t PROBLEM_DIMENSION = 2;


int main(int argc, char *argv[]) {
#ifdef MSK_EXAMPLES_USE_MPI
    int mpi_thread_req{MPI_THREAD_MULTIPLE}, mpi_thread_prov{0};
    MPI_Init_thread(&argc, &argv, mpi_thread_req, &mpi_thread_prov);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    assert(mpi_thread_prov >= mpi_thread_req);
#endif

    using DataType = double;

    // Parse problem parameters
    std::string problem_file;
    if (argc > 1) {
        problem_file = argv[1];
    } else {
        problem_file = std::string(MSK_EXAMPLES_DIRECTORY) + "/dense/parameters.toml";
    }
    auto params = toml::parse(problem_file);

    // Build logger
    auto &msk_params = toml::find(params, "msk");
    auto main_logger = MSk::build_logger(msk_params["log"], "MSk", "stdout", "info");

    // Parse matrix parameters
    using MatrixType = MSk::matrix::Dense<DataType, PROBLEM_DIMENSION>;
    auto &matrix_params = toml::find(params, "matrix");
    std::string matrix_fname = toml::find<std::string>(matrix_params, "file");
    std::string matrix_type = toml::find_or<toml::string>(matrix_params, "type", "rect");
    std::string matrix_storage = toml::find_or<toml::string>(matrix_params, "storage", "columns");
    auto mode = MatrixType::StoreMode::columns;
    if (matrix_storage == "rows") {
        mode = MatrixType::StoreMode::rows;
    }
    auto type = MatrixType::Type::rect;
    if (matrix_type == "sym_up") {
        type = MatrixType::Type::sym_up;
    } if (matrix_type == "sym_low") {
        type = MatrixType::Type::sym_low;
    } else if (matrix_type == "tri_up") {
        type = MatrixType::Type::tri_up;
    } else if (matrix_type == "tri_low") {
        type = MatrixType::Type::tri_low;
    }

    // Geometry
    std::string p_row_fname = toml::find<std::string>(params, "geometry", "rows_file");
    std::string p_col_fname = toml::find<std::string>(params, "geometry", "columns_file");

    // Build matrix
    MatrixType matr(matrix_fname, type, mode, p_row_fname, p_col_fname);

    uint64_t n_rows, n_cols;
    std::tie(n_rows, n_cols) = matr.size();
    main_logger->debug("Matrix size: {} x {}", n_rows, n_cols);

    // Build oracle
    auto &oracle_params = toml::find(params, "oracle");
    auto &approx_params = toml::find(params, "approximation");
    using OracleType = MSk::oracle::Distant<DataType, MatrixType>;
    OracleType oracle{matr, MSk::oracle::read_parameters(oracle_params, approx_params)};

    // Do approximation
    auto t_start = std::chrono::steady_clock::now();
#ifdef MSK_EXAMPLES_USE_MPI
    auto msk = MSk::MosaicSkeleton<DataType>::approximate(matr, oracle, omp_get_max_threads(), MPI_COMM_WORLD, main_logger);
#else
    auto msk = MSk::MosaicSkeleton<DataType>::approximate(matr, oracle, omp_get_max_threads(), main_logger);
#endif
    auto t_end = std::chrono::steady_clock::now();
    double t_diff = std::chrono::duration <double, std::milli> (t_end - t_start).count();
    main_logger->info("MSk building time = {}ms", t_diff);

    // Recompress
    double rel_tol = toml::find_or<double>(approx_params, "rel_tol", 1e-3);
    t_start = std::chrono::steady_clock::now();
    msk.recompress(rel_tol);
    t_end = std::chrono::steady_clock::now();
    t_diff = std::chrono::duration <double, std::milli> (t_end - t_start).count();
    main_logger->info("MSk recompression time = {}ms", t_diff);

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

#ifdef MSK_EXAMPLES_USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
