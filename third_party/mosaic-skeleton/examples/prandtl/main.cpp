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

#include <random>
#include <complex>
#include <functional>

#include "toml.hpp"
#include "common.hpp"
#include "kernel.hpp"

#include "msk/oracle/distant.hpp"
#include "msk/mosaic_skeleton.hpp"
#include "msk/utility/testing.hpp"
#include "msk/utility/logging.hpp"




template<class DataType>
void read_vector(const std::string &fname, uint64_t n, DataType *x) {
    std::ifstream fin(fname, std::ios_base::in);
    for (uint64_t k = 0; k < n; ++k) {
        fin >> x[k];
    }
    // std::ifstream fin(fname, std::ios_base::in | std::ios_base::binary);
    // fin.read(reinterpret_cast<char *>(x), sizeof(DataType) * n);
}

template<class DataType>
void write_vector(const std::string &fname, uint64_t n, const DataType *x) {
    std::ofstream fout(fname, std::ios_base::out);
    for (uint64_t k = 0; k < n; ++k) {
        fout << x[k] << '\n';
    }
}


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

    using DataType = double;
    using RealType = decltype(std::abs(DataType{}));

    // Parse problem parameters
    std::string problem_file;
    if (argc > 1) {
        problem_file = argv[1];
    } else {
        problem_file = std::string(MSK_EXAMPLES_DIRECTORY) + "/prandtl/parameters.toml";
    }
    auto params = toml::parse(problem_file);
    
    auto &msk_params = toml::find(params, "msk");
    auto main_logger = MSk::build_logger(msk_params["log"], "MSk", "stdout", "info");

    for (const uint64_t N : {50, 100, 200, 300, 400, 500, 600}) {
        for (const auto approx_tol : {1e-4, 1e-5, 1e-6, 1e-7}) {

            const uint64_t n_pts = N * N;
            const std::string geom_type = toml::find_or(params, "geometry", "type", "uniform");

            // Build functional matrix
            using MatrixType = Prandtl<DataType>;
            MatrixType matr(N, N, geom_type);


            // Build oracle
            auto &oracle_params = toml::find(params, "oracle");
            auto &approx_params = toml::find(params, "approximation");
            using OracleType = MSk::oracle::Distant<DataType, MatrixType>;
            OracleType oracle{matr, MSk::oracle::read_parameters(oracle_params, approx_params)};
            oracle.min_rows(2 * N);
            oracle.min_cols(2 * N);
            oracle.block_tolerance(approx_tol);
            main_logger->info("Approximate for N = {}, tol = {:.2e}", N, approx_tol);

            // Do approximation
            uint64_t n_jobs = toml::find<toml::integer>(params, "msk", "n_jobs");
            auto t_start = std::chrono::steady_clock::now();
#ifdef MSK_EXAMPLES_USE_MPI
            auto msk = MSk::MosaicSkeleton<DataType>::approximate(matr, oracle, n_jobs, comm, main_logger);
#else
            auto msk = MSk::MosaicSkeleton<DataType>::approximate(matr, oracle, n_jobs, main_logger);
#endif
            auto t_end = std::chrono::steady_clock::now();
            double t_diff = std::chrono::duration <double> (t_end - t_start).count();
            if (proc_id == 0) {
                main_logger->info("Approximation time = {:.3e}s", t_diff);
            }
            double approx_time = t_diff;

            // Do recompression
            // auto &approx_params = toml::find(params, "approximation");
            double rel_tol = toml::find_or<double>(approx_params, "rel_tol", 1e-3);
            t_start = std::chrono::steady_clock::now();
            msk.recompress(approx_tol);
            t_end = std::chrono::steady_clock::now();
            t_diff = std::chrono::duration<double>(t_end - t_start).count();
            if (proc_id == 0) {
                main_logger->info("Recompression time = {:.3e}s", t_diff);
            }
            approx_time += t_diff;

            msk.dump_blocks("./logs/blocks.txt", true);

            // Do checks
            auto &check_params = toml::find(params, "check");
            if (check_params.contains("compression") && check_params["compression"].as_boolean()) {
                MSk::check_compression<DataType>(msk, main_logger);
            }
            if (check_params.contains("error") && check_params["error"].as_boolean()) {
                MSk::check_error<DataType, MatrixType>(msk, matr, main_logger);
            }
            if (check_params.contains("matvec") && check_params["matvec"].as_boolean()) {
                MSk::check_matvec<DataType, MatrixType>(msk, matr, main_logger, false, 10);
            }
            double matvec_time{0};
            
            const std::string pre = "./examples/prandtl/data/", post = "_" + std::to_string(N) + "_chebyshev.txt";
            /*std::string solution_file = toml::find_or<std::string>(params, "solve", "solution", std::string{""});
            if (solution_file.empty()) {
                main_logger->warn("Can't read solution (X) file: {}", solution_file);
            }*/
            const std::string solution_file = pre + "exact_sol" + post;
            /*std::string eff_solution_file = toml::find_or<std::string>(params, "solve", "eff_solution", std::string{""});
            if (eff_solution_file.empty()) {
                main_logger->warn("Can't read effective solution (X_eff) file: {}", eff_solution_file);
            }*/
            std::string eff_solution_file = pre + "sol" + post;
            /*std::string rhs_file = toml::find_or<std::string>(params, "solve", "rhs", std::string{""});
            if (rhs_file.empty()) {
                main_logger->warn("Can't read rhs (B) file: {}", rhs_file);
            }*/
            std::string rhs_file = pre + "rhs" + post;
            /*std::string msk_solution_file = toml::find_or<std::string>(params, "solve", "msk_solution", std::string{""});
            if (eff_solution_file.empty()) {
                main_logger->warn("Can't find msk solution file name: {}", eff_solution_file);
            }*/

            // Generate random solution X and corresponding rhs B = AX
            uint64_t local_row_size{msk.local_row_size()}, local_col_size{msk.local_col_size()};
            SmartArray<DataType> B(local_row_size), B_true(local_row_size);
            SmartArray<DataType> X(local_col_size), X_true(local_col_size), X_eff(local_col_size);
            read_vector(solution_file, n_pts, X_true.data());
            read_vector(rhs_file, n_pts, B_true.data());
            read_vector(eff_solution_file, n_pts, X_eff.data());

            for (auto gmres_tol : {1e-2, 1e-3, 1e-4}) {
            main_logger->info("Solve for tol = {:.2e}", gmres_tol);
            constexpr DataType zero{0};
            BLAS::copy(local_row_size, B_true.data(), 1, B.data(), 1);
            BLAS::copy(local_col_size, &zero, 0, X.data(), 1);

            // Get GMRES parameters
            MSk::solvers::linsys::Parameters<RealType> solve_params;
            MSk::set_gmres_parameters<RealType>(params["solve"], solve_params);
            solve_params.set_relative_tolerance(gmres_tol);

            // Solve linear system AX = B using GMRES
            uint64_t n_iter{0};

            t_start = std::chrono::steady_clock::now();
            run_gmres(msk, 1, B.data(), local_row_size, X.data(), local_col_size, solve_params);
            t_end = std::chrono::steady_clock::now();
            t_diff = std::chrono::duration<double>(t_end - t_start).count();
            double solve_time = t_diff;
            std::string msk_solution_file = pre + "msk_sol_" + std::to_string(gmres_tol) + post;
            write_vector(msk_solution_file, local_col_size, X.data());

            // Check maximal error and residual
            SmartArray<DataType> Ax(local_row_size), x(local_col_size);
            double max_res{0}, max_err{0}, max_eff_err{0};
            {
                // Error
                BLAS::copy(local_col_size, X.data(), 1, x.data(), 1);
                BLAS::axpy(local_col_size, DataType(-1.0), X_true.data(), 1, x.data(), 1);
                double err = BLAS::nrm2(local_col_size, x.data(), 1);
                double nrm = BLAS::nrm2(local_col_size, X_true.data(), 1);
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

                // Another method diff
                BLAS::copy(local_col_size, X.data(), 1, x.data(), 1);
                BLAS::axpy(local_col_size, DataType(-1.0), X_eff.data(), 1, x.data(), 1);
                err = BLAS::nrm2(local_col_size, x.data(), 1);
                nrm = BLAS::nrm2(local_col_size, X_eff.data(), 1);
#ifdef MSK_EXAMPLES_USE_MPI
                err = err * err;
                nrm = nrm * nrm;
                MPI::Allreduce_inplace_sum(&err, 1, comm);
                MPI::Allreduce_inplace_sum(&nrm, 1, comm);
                err = std::sqrt(err);
                nrm = std::sqrt(nrm);
#endif
                err /= nrm;
                max_eff_err = std::max(err, max_eff_err);
                
                // Residual
                BLAS::copy(local_row_size, &zero, 0, Ax.data(), 1);
                msk.matvec(X.data(), Ax.data(), false);
                BLAS::axpy(local_row_size, DataType(-1.0), B_true.data(), 1, Ax.data(), 1);
                double res = BLAS::nrm2(local_row_size, Ax.data(), 1);
                nrm = BLAS::nrm2(local_row_size, B_true.data(), 1);

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
                main_logger->info("Error: {:.4e}", max_err);
                main_logger->info("Effective error: {:.4e}", max_eff_err);
                main_logger->info("Residual: {:.4e}", max_res);
            }

            main_logger->info("TOTAL: {} {:.2e} {:.2e} {} {} {} {} {} {} {} ",
                    N, approx_tol, gmres_tol, 
                    approx_time, matvec_time, solve_time,
                    n_iter, max_err, max_eff_err, max_res);

            }
        }
    }

#ifdef MSK_EXAMPLES_USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
