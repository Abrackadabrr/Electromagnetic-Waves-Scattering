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
#include <string>
#include <tuple>
#include <vector>
#ifdef MSK_EXAMPLES_USE_MPI
#include <mpi.h>
#endif

#include <cmath>
#include <fstream>
#include <random>
#include <chrono>
#include <tuple>
#include <functional>

#include "toml.hpp"
#include "common.hpp"

#include "cpp_blas/BLAS.hpp"
#include "smart_array/smart_array.hpp"

#include "magnetic.hpp"
#include "msk/oracle/distant.hpp"
#include "msk/mosaic_skeleton.hpp"
#include "msk/utility/testing.hpp"
#include "msk/utility/logging.hpp"
#include "msk/solvers/linsys/CG.hpp"
#include "msk/solvers/linsys/precondition.hpp"


template<typename DataType>
class DenseMatrix {
private:
    const MSk::block::Block<DataType> &data;
    uint64_t m{}, n{};

public:
    DenseMatrix(const MSk::block::Block<DataType> &block) : data{block} {
        std::tie(m, n) = block.size();
    }

    void matvec(const DataType *x, DataType *y) const {
        constexpr DataType zero{0};
        BLAS::copy(m, &zero, 0, y, 1);
        data.matvec(x, y);
    }

    void conj_matvec(const DataType *x, DataType *y) const {
        constexpr DataType zero{0};
        BLAS::copy(n, &zero, 0, y, 1);
        data.conj_matvec(x, y);
    }

    uint64_t get_row_number() const {
        return m;
    }

    uint64_t get_column_number() const {
        return n;
    }

    std::pair<uint64_t, uint64_t> size() const {
        return {m, n};
    }

    decltype(std::abs(DataType{})) norm() const {
        return data.frob_norm();
    }
};


template<class DataType>
int read_vector(const std::string &fname, uint64_t n, DataType *x) {
    std::ifstream fin(fname, std::ios_base::in | std::ios_base::binary);
    if (!fin) {
        return -1;
    }
    fin.read(reinterpret_cast<char *>(x), sizeof(DataType) * n);
    return 0;
}

template<class DataType>
int write_vector(const std::string &fname, uint64_t n, const DataType *x) {
    std::ofstream fout(fname, std::ios_base::out | std::ios_base::binary);
    if (!fout) {
        return -1;
    }
    fout.write(reinterpret_cast<const char *>(x), sizeof(DataType) * n);
    return 0;
}



template<typename DataType, class Matrix>
void make_matrix_svd(const Matrix &matrix, const std::string &prefix = "", bool write_dense = false) {
    using RealType = decltype(std::abs(DataType()));
    uint64_t n_rows{0}, n_cols{0};
    std::tie(n_rows, n_cols) = matrix.size();

    SmartArray<DataType> raw_dense(n_rows * n_cols);
    matrix.to_dense(raw_dense.data(), n_rows);
    if (write_dense) {
        write_vector(prefix + "matrix.bin", n_rows * n_cols, raw_dense.data());
    }

    SmartArray<DataType> U_true, V_true;
    SmartArray<RealType> S_true;
    auto true_svd_result = MSk::perform_svd(n_rows, n_cols, raw_dense, U_true, S_true, V_true, {0.0, 1e-12, 0});
    uint64_t true_svd_rank = true_svd_result.memory / (n_rows + n_cols);
    write_vector("U.bin", n_rows * true_svd_rank, U_true.data());
    write_vector("V.bin", n_cols * true_svd_rank, V_true.data());
    write_vector("S.bin", true_svd_rank, S_true.data());
}



template<typename DataType, class Matrix>
void generate_const(const Matrix &matrix, DataType *X, DataType *B) {
    uint64_t n_rows, n_cols;
    std::tie(n_rows, n_cols) = matrix.size();
    for (uint64_t i = 0; i < n_cols; ++i) {
        X[i] = DataType{1} / std::sqrt(n_cols);
    }
    matrix.matvec(X, B);
}

template<typename DataType, class Matrix>
void generate_from_file(const Matrix &matrix, DataType *X, DataType *B, const std::string &filename) {
    uint64_t n_rows, n_cols;
    std::tie(n_rows, n_cols) = matrix.size();
    read_vector(filename, n_cols, X);
    matrix.matvec(X, B);
}

template<typename DataType, class Matrix>
void generate_conj(const Matrix &matrix, DataType *X, DataType *B) {
    std::normal_distribution<double> normal(0, 1);
    std::default_random_engine eng(42);
    uint64_t n_rows, n_cols;
    std::tie(n_rows, n_cols) = matrix.size();
    SmartArray<DataType> coefs{n_rows};
    for (uint64_t i = 0; i < n_rows; ++i) {
        coefs[i] = normal(eng);
    }
    matrix.conj_matvec(coefs.data(), X);
    using RealType = decltype(std::abs(DataType()));
    RealType x_nrm = BLAS::nrm2(n_cols, X, 1);
    BLAS::rscal(n_cols, RealType{1.0} / x_nrm, X, 1);
    matrix.matvec(X, B);
}

template<class DataType, class Matrix>
std::tuple<decltype(std::abs(DataType())), decltype(std::abs(DataType()))> 
calculate_error_residual(const DataType *Y_guess, const DataType *Y_base, const Matrix &ref_matrix) {
    using RealType = decltype(std::abs(DataType()));
    constexpr DataType zero{0};
    uint64_t n_rows{0}, n_cols{0};
    std::tie(n_rows, n_cols) = ref_matrix.size();
    SmartArray<DataType> Y{n_cols};
    BLAS::copy(n_cols, Y_guess, 1, Y.data(), 1);
    BLAS::axpy(n_cols, DataType{-1.0}, Y_base, 1, Y.data(), 1);
    // error
    RealType err = BLAS::nrm2(n_cols, Y.data(), 1) / BLAS::nrm2(n_cols, Y_base, 1);

    SmartArray<DataType> Z_base{n_rows}, Z_guess{n_rows};
    BLAS::copy(n_rows, &zero, 0, Z_base.data(), 1);
    BLAS::copy(n_rows, &zero, 0, Z_guess.data(), 1);
    ref_matrix.matvec(Y_base, Z_base.data());
    ref_matrix.matvec(Y_guess, Z_guess.data());
    
    // residual
    BLAS::axpy(n_rows, DataType{-1.0}, Z_base.data(), 1, Z_guess.data(), 1);
    RealType res = BLAS::nrm2(n_rows, Z_guess.data(), 1) / BLAS::nrm2(n_rows, Z_base.data(), 1);

    return std::tuple<RealType, RealType>(err, res);
}

template<class DataType, class Matrix, class Prec, class RefMatrix>
void make_experiment_cgnr(
        const Matrix &matrix, const std::string &matrix_name, 
        const Prec &prec, const std::string &prec_name,
        const DataType *rhs,
        const decltype(std::abs(DataType())) lambda,
        const DataType noise_level,
        int seed,
        MSk::solvers::linsys::Parameters<decltype(std::abs(DataType()))> params,
        DataType *sol_num,
        const RefMatrix &matrix_ref,
        const DataType *sol_true, const DataType *sol_ref,
        std::shared_ptr<MSk::log::Logger> log) {
    
    using RealType = decltype(std::abs(DataType()));

    constexpr DataType zero{0};
    uint64_t n_rows{0}, n_cols{0};
    std::tie(n_rows, n_cols) = matrix.size();
    SmartArray<DataType> X(n_cols), B(n_rows);
    BLAS::copy(n_cols, &zero, 0, X.data(), 1);
    
    BLAS::copy(n_rows, rhs, 1, B.data(), 1);

    std::normal_distribution<double> normal(0, 1);
    std::default_random_engine eng(seed);
    SmartArray<DataType> Noise(n_rows);
    for (uint64_t i = 0; i < n_rows; ++i) {
        Noise[i] = normal(eng);
    }
    const DataType noise_scale = 1 / BLAS::nrm2(n_rows, Noise.data(), 1) * noise_level;
    BLAS::axpy(n_rows, noise_scale, Noise.data(), 1, B.data(), 1);

    auto t_start = std::chrono::steady_clock::now();
    int64_t iterations = MSk::solvers::linsys::rCGNR(matrix, prec, B.data(), X.data(), params, lambda * lambda);
    auto t_end = std::chrono::steady_clock::now();
    double t_diff = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    const double sol_norm = BLAS::nrm2(n_cols, X.data(), 1);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, time, {:.3f}, ms", matrix_name, prec_name, lambda, noise_level, t_diff);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, iterations, {}", matrix_name, prec_name, lambda, noise_level, iterations);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, solution norm, {:.3e}", matrix_name, prec_name, lambda, noise_level, sol_norm);

    BLAS::copy(n_cols, X.data(), 1, sol_num, 1);
    
    RealType error{0}, residual{0};
    std::tie(error, residual) = calculate_error_residual(sol_num, sol_true, matrix_ref);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, error vs true, {:.3e}", matrix_name, prec_name, lambda, noise_level, error);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, residual vs true, {:.3e}", matrix_name, prec_name, lambda, noise_level, residual);

    std::tie(error, residual) = calculate_error_residual(sol_num, sol_ref, matrix_ref);



    log->info("solve, {}, {}, {:.3e}, {:.3e}, error vs ref, {:.3e}", matrix_name, prec_name, lambda, noise_level, error);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, residual vs ref, {:.3e}", matrix_name, prec_name, lambda, noise_level, residual);

}


template<class DataType, class RefMatrix>
void make_experiment_svd(
        const uint64_t svd_rank,
        const DataType *U, const decltype(std::abs(DataType())) *S,
        const DataType *V,
        const DataType *rhs,
        const bool cut_sing_tail,
        const decltype(std::abs(DataType())) lambda,
        const DataType noise_level,
        int seed,
        MSk::solvers::linsys::Parameters<decltype(std::abs(DataType()))> params,
        DataType *sol_num,
        const RefMatrix &matrix_ref,
        const DataType *sol_true, const DataType *sol_ref,
        std::shared_ptr<MSk::log::Logger> log) {
    
    using RealType = decltype(std::abs(DataType()));
    const std::string matrix_name{"svd"}, prec_name{"svd"};

    constexpr DataType zero{0};
    uint64_t n_rows{0}, n_cols{0};
    std::tie(n_rows, n_cols) = matrix_ref.size();
    SmartArray<DataType> X(n_cols), B(n_rows);
    BLAS::copy(n_cols, &zero, 0, X.data(), 1);
    
    BLAS::copy(n_rows, rhs, 1, B.data(), 1);

    std::normal_distribution<double> normal(0, 1);
    std::default_random_engine eng(seed);
    SmartArray<DataType> Noise(n_rows);
    for (uint64_t i = 0; i < n_rows; ++i) {
        Noise[i] = normal(eng);
    }
    const DataType noise_scale = 1 / BLAS::nrm2(n_rows, Noise.data(), 1) * noise_level;
    BLAS::axpy(n_rows, noise_scale, Noise.data(), 1, B.data(), 1);

    SmartArray<DataType> SinvUcB{svd_rank};

    auto t_start = std::chrono::steady_clock::now();
    uint64_t rank{0};
    if (cut_sing_tail) {
        rank = 0;
        for (; rank < svd_rank && S[rank] > 2 * lambda; ++rank) {}
        BLAS::gemv('C', n_rows, rank, DataType{1.0}, U, n_rows, B.data(), 1, DataType{0}, SinvUcB.data(), 1);
        for (uint64_t i = 0; i < rank; ++i) {
            SinvUcB[i] /= S[i];
        }
    } else {
        rank = svd_rank;
        BLAS::gemv('C', n_rows, rank, DataType{1.0}, U, n_rows, B.data(), 1, DataType{0}, SinvUcB.data(), 1);
        for (uint64_t i = 0; i < rank; ++i) {
            SinvUcB[i] = S[i] / (S[i] * S[i] + lambda * lambda);
        }
    }
    BLAS::gemv('N', n_cols, rank, DataType{1}, V, n_cols, SinvUcB.data(), 1,
            DataType{0}, X.data(), 1);
    auto t_end = std::chrono::steady_clock::now();
    double t_diff = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    const double sol_norm = BLAS::nrm2(n_cols, X.data(), 1);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, time, {:.3f}, ms", matrix_name, prec_name, lambda, noise_level, t_diff);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, solution norm, {:.3e}", matrix_name, prec_name, lambda, noise_level, sol_norm);

    BLAS::copy(n_cols, X.data(), 1, sol_num, 1);
    
    RealType error{0}, residual{0};
    std::tie(error, residual) = calculate_error_residual(sol_num, sol_true, matrix_ref);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, error vs true, {:.3e}", matrix_name, prec_name, lambda, noise_level, error);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, residual vs true, {:.3e}", matrix_name, prec_name, lambda, noise_level, residual);

    std::tie(error, residual) = calculate_error_residual(sol_num, sol_ref, matrix_ref);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, error vs ref, {:.3e}", matrix_name, prec_name, lambda, noise_level, error);
    log->info("solve, {}, {}, {:.3e}, {:.3e}, residual vs ref, {:.3e}", matrix_name, prec_name, lambda, noise_level, residual);
}




int main(int argc, char *argv[]) {
#ifdef MSK_EXAMPLES_USE_MPI
    int mpi_thread_req{MPI_THREAD_MULTIPLE}, mpi_thread_prov{0};
    MPI_Init_thread(&argc, &argv, mpi_thread_req, &mpi_thread_prov);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    assert(mpi_thread_prov >= mpi_thread_req);
#endif

    using DataType = double;
    using RealType = decltype(std::abs(DataType{}));

    // Parse problem parameters
    std::string problem_file;
    if (argc > 1) {
        problem_file = argv[1];
    } else {
        problem_file = std::string(MSK_EXAMPLES_DIRECTORY) + "/magnetic/parameters.toml";
    }
    auto params = toml::parse(problem_file);

    // Build logger
    auto &msk_params = toml::find(params, "msk");
    auto main_logger = MSk::build_logger(msk_params["log"], "MSk", "stdout", "info");
    auto stat_log = MSk::build_logger(std::string("stat"), "statistics.txt", "debug");

    // Build matrix
    const auto &matrix_params = toml::find(params, "matrix");
    uint64_t matrix_fields_type{0};
    std::string m_type = toml::find_or<toml::string>(matrix_params, "type", "field");

    if (m_type == "grad") {
        matrix_fields_type = 1;
    } else if (m_type == "field_grad" || m_type == "grad_field") {
        matrix_fields_type = 2;
    } else {
        matrix_fields_type = 0;
    }
    main_logger->debug("Matrix type: {}", m_type);

    
    // Geometry for matrix
    std::string rows_file = toml::find<toml::string>(params, "geometry", "rows_file");
    std::string columns_file = toml::find<toml::string>(params, "geometry", "columns_file");

    using MatrixType = MSk::Magnetic<DataType>;
    MatrixType matr(matrix_fields_type, rows_file, columns_file);

    uint64_t n_rows, n_cols;
    std::tie(n_rows, n_cols) = matr.size();
    main_logger->debug("Matrix size: {} x {}", n_rows, n_cols);
    stat_log->info("build, msk, size, {}, {}", n_rows, n_cols);

    // Build oracle
    auto &oracle_params_table = toml::find(params, "oracle");
    auto &approx_params = toml::find(params, "approximation");
    auto oracle_params = MSk::oracle::read_parameters(oracle_params_table, approx_params);
    using OracleType = MSk::oracle::Distant<DataType, MatrixType>;
    OracleType oracle{matr, oracle_params};

    uint64_t min_cols{oracle_params.min_cols};
    if (n_cols == min_cols) {
        oracle.force_lrank(true);
    }

    uint64_t n_jobs = toml::find<toml::integer>(params, "msk", "n_jobs");
    auto t_start = std::chrono::steady_clock::now();
#ifdef MSK_EXAMPLES_USE_MPI
    auto msk = MSk::MosaicSkeleton<DataType>::approximate(matr, oracle, n_jobs, MPI_COMM_WORLD, main_logger);
#else
    auto msk = MSk::MosaicSkeleton<DataType>::approximate(matr, oracle, n_jobs, main_logger);
#endif
    auto t_end = std::chrono::steady_clock::now();
    double t_diff = std::chrono::duration <double, std::milli> (t_end - t_start).count();
    const double t_msk_approx{t_diff};
    main_logger->info("MSk building time = {:.4f}ms", t_diff);

    // Recompress
    const double rel_tol = toml::find_or<double>(approx_params, "rel_tol", 1e-3);
    stat_log->info("build, msk, tol, {:.4e}", rel_tol);
    const double sol_tol = toml::find_or<double>(params, "solve", "rel_tol", 1e-3);
    const double rel_svd = toml::find_or<double>(params, "approximation", "rel_svd", 1e-3);
    stat_log->info("build, svd, tol, {:.4e}", rel_svd);
    uint64_t iter_limit = toml::find_or<int>(params, "solve", "iter_limit", 0);
    t_start = std::chrono::steady_clock::now();
    msk.recompress(rel_tol);
    t_end = std::chrono::steady_clock::now();
    t_diff = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    main_logger->info("MSk recompression time = {:.4f}ms", t_diff);
    stat_log->info("build, msk, time, {:.4f}, ms", t_msk_approx + t_diff);


    // Approximation checks
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

    // Block-SVD solution test
    std::string solution_file = toml::find_or<std::string>(params, "solve", "solution", "file", std::string{""});
    if (solution_file.empty()) {
        main_logger->warn("Can't read solution (X) file: {}", solution_file);
    }
    stat_log->info("solve, solution, filename, {}", solution_file);


    MSk::index::Index<uint64_t> I(0, n_rows), J(0, n_cols);
    auto data_block = new MSk::block::Dense<DataType>{n_rows, n_cols};
    auto dense = new MSk::comp::Dense<DataType>{data_block, false};
    MSk::block::Block<DataType> full_block{I, J, data_block, dense};
    t_start = std::chrono::steady_clock::now();
    full_block.compute(matr);
    t_end = std::chrono::steady_clock::now();
    t_diff = std::chrono::duration <double, std::milli> (t_end - t_start).count();
    main_logger->info("FULL MATRIX building time = {:.4f}ms", t_diff);
    stat_log->info("build, full, time, {:.4f}, ms", t_diff);
    DenseMatrix<DataType> full_matrix{full_block};

    DataType zero{0.0};

    if (min_cols == n_cols) {
        // perform block SVD factorization
        std::vector<const MSk::block::LowRank<DataType> *> lrank_blocks;
        msk.extract_blocks(lrank_blocks);

        // Block_SVD solver
        SmartArray<DataType> U, V;
        SmartArray<RealType> S;
        t_start = std::chrono::steady_clock::now();
        uint64_t svd_rank = MSk::block_svd(lrank_blocks.size(), lrank_blocks, U, S, V, rel_svd);
        t_end = std::chrono::steady_clock::now();
        t_diff = std::chrono::duration<double, std::milli> (t_end - t_start).count();
        main_logger->info("Block-SVD time = {:.4f}ms", t_diff);
        main_logger->info("Block-SVD rank = {}", svd_rank);
        stat_log->info("build, svd, time, {:.4f}, ms", t_diff);
        stat_log->info("build, svd, rank, {}", svd_rank);
        for (uint64_t j = 0; j < svd_rank; ++j) {
            msk.unpermute_row(U.data() + n_rows * j);
        }
        for (uint64_t j = 0; j < svd_rank; ++j) {
            msk.unpermute_col(V.data() + n_cols * j);
        }


        SmartArray<DataType> X_true{n_cols}; // solution
        SmartArray<DataType> B_true{n_rows}; // system rhs: A^* B

        if (!solution_file.empty()) {
            generate_from_file(full_matrix, X_true.data(), B_true.data(), solution_file);
        } else {
            generate_conj(full_matrix, X_true.data(), B_true.data());
            // generate_const(full_matrix, X_true.data(), B_true.data());
        }
        
        SmartArray<DataType> X_ref;
        const std::string ref_solution_file = toml::find_or<std::string>(params, "solve", "solution_ref", "file", std::string{""});
        if (!ref_solution_file.empty()) {
            X_ref.resize(n_cols);
            read_vector(ref_solution_file, n_cols, X_ref.data());
        }
        
        SmartArray<DataType> B{n_rows};
        BLAS::copy(n_rows, B_true.data(), 1, B.data(), 1);

        SmartArray<DataType> VcX{svd_rank};
        t_start = std::chrono::steady_clock::now();
        BLAS::gemv('C', n_cols, svd_rank, DataType{1.0}, V.data(), n_cols,
                   X_true.data(), 1, DataType{0}, VcX.data(), 1);
        t_end = std::chrono::steady_clock::now();
        double VcX_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        RealType nrm = BLAS::nrm2(n_cols, X_true.data(), 1);
        RealType inner = BLAS::nrm2(svd_rank, VcX.data() , 1) / nrm;
        main_logger->info("Norm of sol: inner = {:.3e}, outer = {:.3e}", inner, std::sqrt(1.0 - inner * inner));
        stat_log->info("svd, solution, in norm, {:.4f}, ms", inner);
        stat_log->info("svd, solution, out norm, {:.4f}, ms", std::sqrt(1.0 - inner * inner));
        

        // Solve via SVD: x = V S^{-1} U^* b
        main_logger->info("Solve via SVD. Singular range: {:.3e} -- {:.3e}", S[0], S[svd_rank - 1]);
        SmartArray<DataType> UcB{svd_rank};
        BLAS::gemv('C', n_rows, svd_rank, DataType{1.0}, U.data(), n_rows, B.data(), 1,
                    DataType{0}, UcB.data(), 1);
        nrm = BLAS::nrm2(n_rows, B.data(), 1);
        inner = BLAS::nrm2(svd_rank, UcB.data(), 1) / nrm;
        main_logger->info("Norm of rhs: inner = {:.3e}, outer = {:.3e}", inner, (1.0 - inner * inner));

        auto reg_coefs = toml::get<std::vector<double>>(toml::find(params, "solve", "regular"));
        auto noise_levels = toml::get<std::vector<double>>(toml::find(params, "solve", "noise_levels"));

        SmartArray<DataType> X_svd{n_cols}, X_msk_prec{n_cols}, X_msk_noprec{n_cols}, X_full_prec{n_cols}, X_full_noprec{n_cols};
        
        for (auto noise_level: noise_levels) {
            for (auto lam_rel : reg_coefs) {
                RealType lambda = lam_rel * S[0];

                MSk::solvers::linsys::precond::None<DataType> none_prec{n_cols};
                uint64_t rank{0};
                for (; rank < svd_rank && S[rank] > 2 * lambda; ++rank) {}
                MSk::solvers::linsys::precond::RegNormal<DataType> SVD_prec{n_cols, rank, V.data(), S.data(), lambda};
                MSk::solvers::linsys::Parameters<RealType> params{sol_tol, 0.0, 0};
                int seed = 42;
                make_experiment_cgnr(full_matrix, "full", none_prec, "none", B_true.data(), lambda, noise_level, seed, params, X_full_noprec.data(), full_matrix, X_true.data(), X_ref.data(), stat_log);
                write_vector("sol_full_noprec_" + std::to_string(noise_level) + "_" + std::to_string(lambda) + ".bin", n_cols, X_full_noprec.data());
                make_experiment_cgnr(full_matrix, "full", SVD_prec, "svd", B_true.data(), lambda, noise_level, seed, params, X_full_prec.data(), full_matrix, X_true.data(), X_full_noprec.data(), stat_log);
                make_experiment_svd(svd_rank, U.data(), S.data(), V.data(), B.data(), true, lambda, noise_level, seed, params, X_svd.data(), full_matrix, X_true.data(), X_full_noprec.data(), stat_log);
                make_experiment_cgnr(msk, "msk", none_prec, "none", B_true.data(), lambda, noise_level, seed, params, X_msk_noprec.data(), full_matrix, X_true.data(), X_full_noprec.data(), stat_log);
                make_experiment_cgnr(msk, "msk", SVD_prec, "svd", B_true.data(), lambda, noise_level, seed, params, X_msk_prec.data(), full_matrix, X_true.data(), X_full_noprec.data(), stat_log);
            }
        }
        
        
    }

#ifdef MSK_EXAMPLES_USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
