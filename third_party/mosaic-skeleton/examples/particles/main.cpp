/*
 * Copyright 2024-2026 Bulat Valiakhmetov
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
// Created by bulat_v on 29.01.2024.
//


#ifdef MSK_EXAMPLES_USE_MPI
#include <mpi.h>
#endif

#include <random>
#include <complex>
#include <functional>

#include "toml.hpp"
#include "common.hpp"
#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"

#include "particles.hpp"
#include "msk/matrix/basic.hpp"
#include "msk/oracle/elementary.hpp"
#include "msk/mosaic_skeleton.hpp"
#include "msk/solvers/ode/runge_kutta.hpp"
#include "msk/utility/testing.hpp"
#include "msk/utility/logging.hpp"



template<typename DataType, class Matrix>
class Smoluch {
    const uint64_t _size;
    const Matrix &_matrix;
    const DataType _src{0};
public:
    Smoluch(const uint64_t size, const Matrix &matrix, DataType src = 0) : _size{size}, _matrix{matrix}, _src{src} {}

    void apply(const DataType *x, DataType *y) const {
        SmartArray<DataType> tmp(_size);
        tmp[0] = 0;
        _matrix.convolve(x, x, tmp.data() + 1, _size - 1);
        _matrix.matvec(x, y);
        for (uint64_t i = 0; i < _size; ++i) {
            y[i] = 0.5 * tmp[i] - x[i] * y[i];
        }
        y[0] += _src;
    }
};


int main(int argc, char *argv[]) {

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
    using RealType = decltype(std::abs(DataType(1.0)));

    // Parse problem parameters
    std::string problem_file;
    if (argc > 1) {
        problem_file = argv[1];
    } else {
        problem_file = std::string(MSK_EXAMPLES_DIRECTORY) + "/particles/parameters.toml";
    }
    auto params = toml::parse(problem_file);

    auto &msk_params = toml::find(params, "msk");
    auto main_logger = MSk::build_logger(msk_params["log"], "MSk", "stdout", "info");

    // uint64_t n_pts = toml::find_or(params, "problem", "size", 1024);

    std::function<DataType(uint64_t, uint64_t)> kernel;
    std::string type = toml::find_or(params, "problem", "type", "atmos");
    if (type == "atmosphere" || type == "atmos") {
        main_logger->info("Kernel: Atmospheric");
        kernel = MSk::ker_atmos<RealType>;
    } else if (type == "stream") {
        main_logger->info("Kernel: Stream");
        kernel = MSk::ker_stream<RealType>;
    } else if (type == "fusion") {
        RealType alpha = toml::find_or(params, "problem", "alpha", 1.0);
        main_logger->info("Kernel: Fusion, alpha = {}", alpha);
        kernel = std::bind(MSk::ker_fusion<RealType>, std::placeholders::_1, std::placeholders::_2, alpha);
    } else if (type == "modif" || type == "modified") {
        RealType c = toml::find_or(params, "problem", "c", 1.0);
        main_logger->info("Kernel: Modified, C = {}", c);
        kernel = std::bind(MSk::ker_modified<RealType>, std::placeholders::_1, std::placeholders::_2, c);
    } else if (type == "emul" || type == "emulsion") {
        RealType c1 = toml::find_or(params, "problem", "c1", 1.0);
        RealType c2 = toml::find_or(params, "problem", "c2", 1.0);
        main_logger->info("Kernel: Emulsion, C1 = {}, C2 = {}", c1, c2);
        kernel = std::bind(MSk::ker_emulsion<RealType>, std::placeholders::_1, std::placeholders::_2, c1, c2);
    } else if (type == "const" || type == "constant") {
        main_logger->info("Kernel: Constant, C = 2");
        kernel = [](uint64_t, uint64_t) { return RealType{2}; };
    } else if (type == "powered" || type == "pow") {
        RealType p = toml::find_or(params, "problem", "p", 1.0 / 3.0);
        main_logger->info("Kernel: Powered module, p = {}", p);
        kernel = std::bind(MSk::ker_module_powered<RealType>, std::placeholders::_1, std::placeholders::_2, p);
    } else if (type == "module" || type == "diff") {
        main_logger->info("Kernel: Module");
        kernel = MSk::ker_module<RealType>;
    } else {
        main_logger->info("Kernel: Zero");
        kernel = [](uint64_t, uint64_t) { return RealType{0}; };
    }


    // Build oracle
    auto &oracle_params = toml::find(params, "oracle");
    auto &approx_params = toml::find(params, "approximation");
    using OracleType = MSk::oracle::Elementary<DataType>;
    OracleType oracle{MSk::oracle::read_parameters(oracle_params, approx_params)};

    // uint64_t n_pts = toml::find_or(params, "problem", "size", 1024);
    uint64_t min_pts = toml::find_or(params, "problem", "min_size", 1024);
    uint64_t max_pts = toml::find_or(params, "problem", "max_size", 1024);
    for (uint64_t n_pts = min_pts; n_pts <= max_pts; n_pts *= 2) {
        // Build functional matrix
        using MatrixType = MSk::matrix::Elementary<DataType>;
        MatrixType matr(n_pts, n_pts, kernel);

        // blocks
        const uint64_t row_blocks{4}, col_blocks{4};
        const uint64_t n_blocks{row_blocks * col_blocks};
        assert(n_pts % row_blocks == 0 && "Invalid row blocks number");
        assert(n_pts % col_blocks == 0 && "Invalid column blocks number");
        const uint64_t row_block_size{n_pts / row_blocks}, col_block_size{n_pts / col_blocks};
        std::vector<MSk::block::Block<DataType>> children(n_blocks);
        for (uint64_t r = 0; r < row_blocks; ++r) {
            for (uint64_t c = 0; c < col_blocks; ++c) {
                MSk::index::Index<uint64_t> b_idx(row_block_size * r, row_block_size * (r + 1));
                MSk::index::Index<uint64_t> b_jdx(col_block_size * c, col_block_size * (c + 1));
                children[r * col_blocks + c] = MSk::block::Block<DataType>(b_idx, b_jdx);
            }
        }
        MSk::block::Block<DataType> base_block(MSk::index::Index<uint64_t>(0, n_pts), MSk::index::Index<uint64_t>(0, n_pts));
        MSk::oracle::ActionAnswer<DataType> base_answer(n_blocks, std::move(children));
        oracle.add_manual(std::move(base_block), std::move(base_answer));

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
            MSk::check_compression<DataType>(msk, main_logger);
        }
        if (check_params.contains("error") && check_params["error"].as_boolean()) {
            MSk::check_error<DataType, MatrixType, TRANS>(msk, matr, main_logger);
        }
        if (check_params.contains("matvec") && check_params["matvec"].as_boolean()) {
            MSk::check_matvec<DataType, MatrixType, TRANS>(msk, matr, main_logger, true, 10);
        }
        if (check_params.contains("convolution") && check_params["convolution"].as_boolean()) {
            MSk::check_convolution<DataType, MatrixType, TRANS>(msk, matr, main_logger, true, 10, n_pts);
        }

        uint64_t n_rows, n_cols;
        std::tie(n_rows, n_cols) = msk.size();
        assert(n_rows == n_cols && "Nonsquare matrix ODE solving");
        SmartArray<DataType> y_distrib(n_rows);
        constexpr DataType zero{0};
        BLAS::copy(n_rows, &zero, 0, y_distrib.data(), 1);
        RealType initial = toml::find_or<double>(params, "problem", "initial", 1.0);
        //y_distrib[0] = initial;
        for (uint64_t i = 0; i < n_rows; ++i) {
            y_distrib[i] = RealType{1} / RealType(i + 2);
        }

        RealType source = toml::find_or<double>(params, "problem", "source", 0.0);
        Smoluch<DataType, MSk::MosaicSkeleton<DataType, TRANS>> smoluch(n_rows, msk, source);
        RealType T = toml::find_or(params, "solve", "time", 1.0);
        uint64_t n_steps = toml::find_or(params, "solve", "n_steps", 5 * n_rows);
        RealType tau = T / n_steps;
        main_logger->info("Start solving for T = {}, n = {} (tau = {}), init = {}, src = {}", T, n_steps, tau, initial, source);
        t_start = std::chrono::steady_clock::now();
        MSk::solvers::ode::runge_kutta4_adaptive(n_rows, smoluch, y_distrib.data(), T, tau);
        t_end = std::chrono::steady_clock::now();
        t_diff = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        main_logger->info("Smoluchowsky solve time = {:.4f}ms", t_diff);

        std::string result_filename = toml::find_or(params, "solve", "file", "stdout");
        if (result_filename == "stdout") {
            MSk::out_matrix(1, n_rows, y_distrib.data(), 1);
        } else {
            std::ofstream f(std::string(type) + "_" + std::to_string(n_pts) + "_" + std::to_string(T) + "_" + std::to_string(n_steps) + + "_" + result_filename);
            for (uint64_t i = 0; i < n_rows; ++i) {
                f << y_distrib[i] << '\n';
            }
        }
    }

#ifdef MSK_EXAMPLES_USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
