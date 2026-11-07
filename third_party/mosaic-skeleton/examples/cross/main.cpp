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
/*
 * Bare usage of MSk::comp::Dense for MSk::block::Dense evaluation and
 * MSk::comp::Cross for MSk::block::LowRank evaluation using
 * matrix cross approximation method.
 */


#include <memory>
#include <random>
#include <complex>

#include "common.hpp"
#include "toml.hpp"

#include "cpp_blas/BLAS.hpp"

#include "msk/index/index.hpp"
#include "msk/block/dense.hpp"
#include "msk/block/lowrank.hpp"
#include "msk/computation/lowrank.hpp"
#include "msk/computation/dense.hpp"
#include "msk/matrix/basic.hpp"
#include "msk/utility/logging.hpp"



template<typename DataType>
void log_statistics(uint64_t N, const SmartArray<DataType> &a_true, const MSk::block::LowRank<DataType> &cross_block,
                      MSk::comp::Result approx_res, std::shared_ptr<MSk::log::Logger> log) {
    using RealType = decltype(std::abs(DataType(1.0)));
    constexpr DataType zero{0};

    uint64_t n_elements{N * N};
    SmartArray<DataType> a_cross(n_elements);
    BLAS::copy(n_elements, &zero, 0, a_cross.data(), 1);
    cross_block.to_dense(a_cross.data(), N);

    RealType true_f_norm = BLAS::nrm2(n_elements, a_true.data(), 1);
    BLAS::axpy(n_elements, DataType(-1.0), a_true.data(), 1, a_cross.data(), 1);
    RealType err_norm = BLAS::nrm2(n_elements, a_cross.data(), 1);

    log->info("Approximation: rank = {}", approx_res.memory / (2 * N));
    log->info("Approx result: err = {:.4e}, nrm = {:.4e}", approx_res.err_norm / approx_res.matr_norm, approx_res.matr_norm);
    log->info("True result: err = {:.3e}, nrm = {:.4e}", err_norm / true_f_norm, true_f_norm);
}

int main(int argc, char *argv[]) {
    using DataType = std::complex<double>;
    using RealType = decltype(std::abs(DataType(1.0)));

    // Parse problem's parameters
    std::string problem_file;
    if (argc > 1) {
        problem_file = argv[1];
    } else {
        problem_file = std::string(MSK_EXAMPLES_DIRECTORY) + "/cross/parameters.toml";
    }
    auto problem = toml::parse(problem_file);

    auto main_logger = MSk::build_logger(problem["log"], "Cross", "stdout", "info");

    // Generate matrix of required size and type
    uint64_t N = toml::find_or(problem, "matrix", "size", 1000);
    uint64_t n_elements{N * N};
    std::string matrix_type = toml::find_or(problem, "matrix", "type", std::string{"bilinear"});
    MSk::matrix::Matrix<DataType> *matrix_ptr;
    if (matrix_type == "hilbert") {
        matrix_ptr = new MSk::matrix::Elementary<DataType>{N, N, MSk::hilbert<DataType>};
    } else {
        matrix_ptr = new MSk::matrix::BiLinear<DataType>{N, N};
    }

    // Get approximation parameters
    auto &approx_params = toml::find(problem, "approximation");
    RealType abs_tol = toml::find_or<toml::floating>(approx_params, "abs_tol", 0.0);
    RealType rel_tol = toml::find_or<toml::floating>(approx_params, "rel_tol", 1e-3);
    uint64_t max_rank = toml::find_or<toml::integer>(approx_params, "max_rank", 0);
    MSk::comp::Parameters compress_params{abs_tol, rel_tol, max_rank};

    // Get recompression parameters
    auto &recompression_params = toml::find(problem, "recompression");
    abs_tol = toml::find_or<toml::floating>(recompression_params, "abs_tol", 0.0);
    rel_tol = toml::find_or<toml::floating>(recompression_params, "rel_tol", 1e-3);
    max_rank = toml::find_or<toml::integer>(recompression_params, "max_rank", 0);
    MSk::comp::Parameters recompress_params{abs_tol, rel_tol, max_rank};

    // Compute full matrix for future checks
    MSk::index::Index<uint64_t> I(N), J(N);
    auto dense_block = new MSk::block::Dense<DataType>{N, N};
    auto denser = new MSk::comp::Dense<DataType>{dense_block};
    denser->compute(*matrix_ptr, I, J); // actual evaluation
    SmartArray<DataType> a_true{n_elements};
    dense_block->to_dense(a_true.data(), N);

    // Approximated matrix data
    auto cross_block = new MSk::block::LowRank<DataType>{N, N};
    auto cross = new MSk::comp::Cross<DataType>{cross_block, compress_params};

    // Data for elements evaluation
    MSk::eval::Request<DataType> request;
    MSk::comp::Result res;
    std::vector<uint64_t> idx{N}, jdx{N};

    main_logger->info("Start Cross on size {} x {}", N, N);

    uint64_t rank{0}, iter{0};
    while (!res.finished) {
        res = cross->step_enqueue(request);
        uint64_t new_rank = res.memory / (2 * N);
        if (rank != new_rank) {
            main_logger->debug("Iter = {}: rank = {}, err = {:.3e}, nrm = {:.4e}", iter, new_rank, res.err_norm, res.matr_norm);
            rank = new_rank;
        }
        // compute elements
        DataType *ptr{request.result_ptr()};
        uint64_t n_elems{request.size()};
        request.extract_list(idx, jdx);
        for (uint64_t i = 0; i < n_elems; ++i) {
            ptr[i] = matrix_ptr->element(idx[i], jdx[i]);
        }
        request.free();

        ++iter;
    }
    log_statistics(N, a_true, *cross_block, res, main_logger);
    main_logger->info("Perform recompression");
    res = cross_block->recompress(recompress_params);
    log_statistics(N, a_true, *cross_block, res, main_logger);

    delete matrix_ptr;
    delete dense_block;
    delete denser;
    delete cross_block;
    delete cross;

    return 0;
}
