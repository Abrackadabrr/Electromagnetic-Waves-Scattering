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


/**
 * @file
 *
 * @brief Defines mosaic-skeleton method of a functional matrix approximation.
 */


#ifndef MOSAIC_SKELETON_MOSAIC_SKELETON_HPP
#define MOSAIC_SKELETON_MOSAIC_SKELETON_HPP


#include "msk_config.h"

#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <sys/types.h>
#include <utility>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <tuple>

#ifdef MSK_USE_MPI
#include "cpp_mpi/MPI.hpp"
#endif

#include "msk/index/index.hpp"
#include "msk/block/block.hpp"
#include "msk/block/lowrank.hpp"
#include "msk/evaluation/queue.hpp"
#include "msk/block/manager.hpp"
#include "msk/thread_pool/thread_pool.hpp"
#include "msk/utility/logging.hpp"
#include "msk/utility/concurrent.hpp"
#ifdef MSK_USE_FFT
#include "msk/utility/fourier.hpp"
#endif
// For specific templates
#include "msk/matrix/basic.hpp"
#include "msk/oracle/distant.hpp"
#include "msk/oracle/elementary.hpp"


namespace MSk {


enum class Trans {
    gen,
    symm,
    herm,
};


/**
 * @brief Mosaic-skeleton approximation of an arbitrary matrix.
 *
 * @details The class manages the approximation process and provides the interface to matrix
 * in mosaic-skeleton format.
 *
 * @tparam DataType Matrix element type.
 * @tparam Trans Matrix type of storage.
 */
template<typename DataType, Trans msk_trans = Trans::gen>
class MosaicSkeleton {
private:
    using RealType = decltype(std::abs(DataType()));
    using MatrixIdx = index::Index<uint64_t>; ///< Interface type for indexing over matrix rows and columns
    using MatrixBlock = block::Block<DataType, MatrixIdx>; ///< Block type

    uint64_t n_rows, n_cols; ///< Matrix sizes

    std::unique_ptr<index::Permutation<uint64_t>> _p_row_perm, _p_col_perm; ///< Rows and columns indices permutation

    uint64_t n_blocks{0}; ///< Total number of matrix subblocks
    mutable block::Manager<DataType> _manager; ///< Blocks manager

    RealType f_norm_sq{0}; ///< Frobenius norm of approximated matrix
    uint64_t _memory{0}; ///< Number of parameters for storing approximated matrix

    mutable std::unique_ptr<exec::ThreadPool> p_pool; ///< Pool of threads for asynchronous approximation

#ifdef MSK_USE_MPI
    MPI_Comm msk_comm; ///< MPI communicator
#endif
    SmartArray<int> local_row_offs, local_col_offs, local_row_sizes, local_col_sizes; ///< Vector slicing data for matvec
    int proc_id{0}, n_procs{1};

    std::shared_ptr<log::Logger> msk_log; ///< Logger

    struct BlockShape {
    public:
        std::vector<std::pair<uint64_t, uint64_t>> row_sizes_offs{};
        uint64_t col_size{0};
        uint64_t col_ldim{0};
    public:
        BlockShape() {}
        BlockShape(uint64_t m, uint64_t i0, uint64_t n = 1, uint64_t ld = 0) :
                   row_sizes_offs{{m, i0}}, col_size{n}, col_ldim{ld} {}

        uint64_t size() const {
            uint64_t rows_total{0};
            for (auto rs : row_sizes_offs) {
                rows_total += rs.first;
            }
            return rows_total * col_size;
        }
    };

    enum class MatvecType {
        direct,
        conj,
        transp,
    };


private:
    template<class Matrix>
    void approximate(Matrix &matrix) {
        // Prepare Queue for elements evaluation
        const std::string eq_log_fname = "./logs/log_eq_" + std::to_string(proc_id) + ".txt";
        auto eq_logger = log::build_logger("EQ_" + std::to_string(proc_id), eq_log_fname);
        eq_logger->set_level(log::Level::info);
        uint64_t max_block_len{0};
        for (auto block_it = _manager.begin(false); block_it != _manager.end(false); ++block_it) {
            uint64_t m, n;
            std::tie(m, n) = block_it->size();
            max_block_len = std::max(max_block_len, std::max(m, n));
        }
        uint64_t batch_size = std::min(std::max(max_block_len, 10000ul), 50000ul);
        if (msk_log) msk_log->debug("batch_size = {}", batch_size);
        // The Queue
        eval::Queue<DataType> elements_queue{&_manager, batch_size, eq_logger};

        // Start Queue evaluation
        elements_queue.start(matrix, p_pool.get());

        // Synchronize approximation procedure
#ifdef MSK_USE_MPI
        MPI_Barrier(msk_comm);
#endif

        // Start approximation of the blocks, locks until it's done
        _manager.start_approximation(p_pool.get(), &elements_queue);

        // Stop elements evaluation
        elements_queue.stop(true);
        p_pool->batch_size(1);

        // The approximation has finished
        if (msk_log && proc_id == 0) msk_log->debug("Finished approximation");

        // Postprocessing
        // Rebalance blocks between processes
        _manager.balance_blocks();
        if (msk_log && proc_id == 0) msk_log->debug("Finished balancing");
#ifdef MSK_USE_MPI
        MPI_Barrier(msk_comm);
#endif
        // Calculate statistics
        refresh_stat();

#ifdef MSK_USE_FFT
        uint64_t n_lrank{0}, max_rank{0}, max_size{0};
        for (auto block_it = _manager.begin(false); block_it != _manager.end(false); ++block_it) {
            auto ptr = dynamic_cast<block::LowRank<DataType> *>(block_it->data());
            if (ptr != nullptr) {
                ++n_lrank;
                uint64_t m, n;
                std::tie(m, n) = block_it->size();
                uint64_t r = ptr->rank();
                max_size = std::max(max_size, m + n - 1);
                max_rank = std::max(max_rank, r);
            }
        }
        auto& daemon = FFT::Daemon::get();
        daemon.init_memory<DataType>(p_pool->max_threads() + 1, max_size, max_rank);
#endif
    }

    /**
     * @brief Generate mosaic partitioning for the matrix.
     *
     * @details Recursively processes blocks.
     * Asks oracle about the block and gathers evaluatable blocks.
     * If the block should be divided then recursively processes subblocks.
     *
     * @tparam Oracle Oracle entity type.
     *
     * @param[in] oracle Oracle that has the knowledge for matrix partitioning.
     */
    template<class Oracle>
    void build_mosaic(Oracle &oracle) {
        // Blocks for future processing
        std::deque<MatrixBlock> remaining_blocks;

        // Resulting blocks to be computed at the approximation stage
        std::vector<MatrixBlock> computable_blocks;

        auto process_next_block = [this, &oracle, &remaining_blocks, &computable_blocks]() {
            static uint64_t next_id{0};
            MatrixBlock &block = remaining_blocks.front();

            // Ask Oracle
            auto action_ans = oracle.get_action(&block);

            if (action_ans.data != nullptr && action_ans.computer != nullptr) { // If block is to compute
                block.set_id(n_blocks);
                ++n_blocks;
                computable_blocks.push_back(std::move(block));
            }
            remaining_blocks.pop_front();

            // Set children for processing if there are any
            for (uint64_t i = 0; i < action_ans.n_child; ++i) {
                action_ans.child[i].set_id(++next_id);
                remaining_blocks.push_back(std::move(action_ans.child[i]));
            }
        };

        // Initialize the first block as the matrix itself
        MatrixBlock starting(MatrixIdx(0, n_rows, _p_row_perm.get()), MatrixIdx(0, n_cols, _p_col_perm.get()));
        starting.set_id(0);
        remaining_blocks.emplace_back(std::move(starting));

        _p_col_perm->lock(true);
        _p_row_perm->lock(true);
        while (!remaining_blocks.empty()) {
            process_next_block();
        }
        // Unlock permutations since they're immutable
        _p_col_perm->lock(false);
        _p_row_perm->lock(false);

        // Transfer computable blocks to Manager
        computable_blocks.shrink_to_fit();
#ifdef MSK_USE_MPI
        MPI_Barrier(msk_comm);
#endif
        _manager.init_blocks(std::move(computable_blocks));
        if (msk_log && proc_id == 0) msk_log->debug("Mosaic is built");
    }

#ifdef MSK_USE_MPI
    void _gather_vector(uint64_t size, DataType *vec, const int *local_sizes, const int *local_offs) const {
        MPI::Allgatherv_inplace(vec, size, local_sizes, local_offs, msk_comm);
    }

    void _gather_block_vector(uint64_t block_size, uint64_t global_size, const DataType *local_block, DataType *global_block, const int *local_sizes, const int *local_offs) const {
        if (block_size == 0) {
            return;
        }
        if (n_procs <= 1) {
            BLAS::copy(global_size * block_size, local_block, 1, global_block, 1);
            return;
        }
        if (block_size == 1) {
            MPI_Allgatherv(local_block, local_sizes[proc_id], MPI::data_type<DataType>(), 
                           global_block, local_sizes, local_offs, MPI::data_type<DataType>(), msk_comm);
            return;
        }

        SmartArray<MPI_Datatype> block_types(n_procs);
        SmartArray<MPI_Request> requests(n_procs);
        for (int root = 0; root < n_procs; ++root) {
            const int local_size = local_sizes[root], local_off = local_offs[root];
            MPI_Type_vector(block_size, local_size, global_size, MPI::data_type<DataType>(), &block_types[root]);
            MPI_Type_commit(&block_types[root]);
            requests[root] = MPI_REQUEST_NULL;
            if (root == proc_id) {
                MPI_Ibcast(const_cast<DataType *>(local_block), local_size * block_size, MPI::data_type<DataType>(), root, msk_comm, &requests[root]);
            } else {
                MPI_Ibcast(global_block + local_off, 1, block_types[root], root, msk_comm, &requests[root]);
            }
        }

        const int local_size = local_sizes[proc_id], local_off = local_offs[proc_id];
        for (uint64_t k = 0; k < block_size; ++k) {
            BLAS::copy(local_size, local_block + k * local_size, 1, global_block + local_off + k * global_size, 1);
        }

        SmartArray<MPI_Status> statuses(n_procs);
        MPI_Waitall(n_procs, requests.data(), statuses.data());
        for (int root = 0; root < n_procs; ++root) {
            MPI_Type_free(&block_types[root]);
        }
        MPI_Barrier(msk_comm);
    }

    void _distreduce_block_vector(uint64_t block_size, uint64_t global_size, DataType *global_block, DataType *local_block, const int *local_sizes, const int *local_offs) const {
        if (n_procs <= 1 || block_size == 0) {
            return;
        }
        if (block_size == 1) {
            const int local_size = local_sizes[proc_id], local_off = local_offs[proc_id];            
            SmartArray<MPI_Request> result_req(n_procs);
            for (int root = 0; root < n_procs; ++root) {
                result_req[root] = MPI_REQUEST_NULL;
                MPI_Ireduce(global_block + local_offs[root], local_block, local_sizes[root], 
                    MPI::data_type<DataType>(), MPI_SUM, root, msk_comm, &result_req[root]);
            }
            SmartArray<MPI_Status> result_stat(n_procs);
            MPI_Waitall(n_procs, result_req.data(), result_stat.data());
            return;
        }

        SmartArray<MPI_Request> requests(block_size);
        for (uint64_t k = 0; k < block_size; ++k) {
            requests[k] = MPI_REQUEST_NULL;
            MPI_Ireduce_scatter(global_block + k * global_size, local_block + k * local_sizes[proc_id], local_sizes, MPI::data_type<DataType>(), MPI_SUM, msk_comm, &requests[k]);
        }

        SmartArray<MPI_Status> statuses(block_size);
        MPI_Waitall(block_size, requests.data(), statuses.data());
    }


#endif
    
    void permute_block_vector(uint64_t size, uint64_t block_size, const DataType *V_in, DataType *V_out, const std::unique_ptr<index::Permutation<uint64_t>> &perm, bool direct) const {
        if (block_size == 1) {
            if (direct) {
                perm->permute(V_in, V_out);
            } else {
                perm->unpermute(V_in, V_out);
            }
            return;
        }

        std::mutex m_perm;
        std::condition_variable cv_perm;
        uint64_t count_perm{block_size};

        auto permute_cols = [&m_perm, &cv_perm, &count_perm, &perm, direct](const DataType *input, DataType *output) {
            if (direct) {
                perm->permute(input, output);
            } else {
                perm->unpermute(input, output);
            }
            std::unique_lock<std::mutex> lock{m_perm};
            --count_perm;
            if (count_perm == 0) {
                lock.unlock();
                cv_perm.notify_one();
            }
        };
        for (uint64_t k = 0; k < block_size; ++k) {
            p_pool->execute(permute_cols, V_in + k * size, V_out + k * size);
        }
        std::unique_lock<std::mutex> lock{m_perm};
        cv_perm.wait(lock, [&count_perm]() {return count_perm == 0;});
    }


    void apply_batched_block_linear(const std::function<BlockShape(MatrixBlock &, DataType *)> &block_operation,
            DataType *result = nullptr, const std::function<bool(const MatrixBlock &)> &condition = [](const MatrixBlock &) {return true;}) const {

        using TBlockIt = decltype(_manager.begin(true));
        struct BlockBatch {
            TBlockIt begin{}, end{};
            typename exec::MemoryBank<DataType>::Cell result{};
        };

        // Synchronization
        class SyncQueue {
        private:
            mutable std::mutex _m{};
            mutable std::condition_variable _cv{};
            std::deque<BlockBatch> _q{}; // batches for master to process
        public:
            void put(BlockBatch &&batch) {
                std::unique_lock<std::mutex> lock{_m};
                _q.push_back(std::move(batch));
                lock.unlock();
                _cv.notify_one();
            }

            void extract(std::deque<BlockBatch> &to) {
                std::unique_lock<std::mutex> lock{_m};
                _cv.wait(lock, [this]{return !_q.empty();});
                to.insert(to.end(), _q.begin(), _q.end());
                _q.clear();
            }
        };
        SyncQueue queue;
        // Block-local memory
        exec::MemoryBank<DataType> bank;

        auto batch_op = [&block_operation, &condition, &bank]
                        (BlockBatch batch, bool output, SyncQueue *queue) {
            DataType *result_ptr{nullptr};
            if (output) { // Prepare local memory for a batch, if required
                uint64_t total_size{0};
                // Accumulate result total size
                for (auto it = batch.begin; it != batch.end; ++it) {
                    if (condition(*it)) {
                        auto block_shape = block_operation(*it, nullptr);
                        total_size += block_shape.size();
                    }
                }
                if (total_size > 0) {
                    batch.result = std::move(bank.get(total_size));
                    result_ptr = batch.result.data();
                    constexpr DataType zero{0.0};
                    BLAS::copy(total_size, &zero, 0, result_ptr, 1);
                }
            }

            // Apply operation blockwise
            for (auto it = batch.begin; it != batch.end; ++it) {
                if (condition(*it)) {
                    auto block_shape = block_operation(*it, result_ptr);
                    result_ptr += block_shape.size();
                }
            }
            // Submit completed blocks
            queue->put(std::move(batch));
        };

        // Form and schedule batches
        constexpr uint64_t max_batch_work_size{100000};
        uint64_t batch_work_size{0};
        auto it_begin = _manager.begin(true);
        for (auto it = _manager.begin(true); it <= _manager.end(true); ++it) {
            if (batch_work_size >= max_batch_work_size || it == _manager.end(true)) {
                BlockBatch batch;
                batch.begin = it_begin;
                batch.end = it;
                p_pool->execute(batch_op, std::move(batch), result != nullptr, &queue);

                it_begin = it;
                batch_work_size = 0;
            }
            if (it != _manager.end(true)) {
                if (condition(*it)) {
                    batch_work_size += it->memory();
                }
            }
        }

        // Accumulate bathes' results
        uint64_t blocks_left{_manager.size(true)};

        while (blocks_left > 0) {
            std::deque<BlockBatch> local_batches{};
            queue.extract(local_batches);
            // Extract ready batches

            // Process batches without lock
            if (result != nullptr) {
                for (auto &batch : local_batches) {
                    uint64_t result_off{0};
                    for (auto it = batch.begin; it != batch.end; ++it) {
                        if (condition(*it)) {
                            auto block_shape = block_operation(*it, nullptr);
                            const uint64_t col_size{block_shape.col_size}, col_ldim{block_shape.col_ldim};
                            for (auto size_off : block_shape.row_sizes_offs) {
                                for (uint64_t k = 0; k < col_size; ++k) {
                                    BLAS::axpy(size_off.first, DataType{1.0}, 
                                               batch.result.data() + result_off + k * size_off.first, 1,
                                               result + size_off.second + k * col_ldim, 1);
                                }
                                result_off += size_off.first * col_size;
                            }
                        }
                        --blocks_left;
                    }
                    if (batch.result.size() > 0) {
                        bank.ret(std::move(batch.result));
                    }
                }
            } else {
                for (auto &batch : local_batches) {
                    blocks_left -= std::distance(batch.begin, batch.end);
                }
            }
        }
    }


    void _general_matvec(MatvecType type, const uint64_t block_size, const DataType *X, DataType *Y, bool pre_perm) const {
        if (block_size == 0) {
            return;
        }

        uint64_t X_size{}, Y_size{};
        using TPerm = index::Permutation<uint64_t>;
        const std::unique_ptr<TPerm> *X_perm_ptr{nullptr};
        const std::unique_ptr<TPerm> *Y_perm_ptr{nullptr};
        const int *X_local_sizes{nullptr}, *X_local_offs{nullptr};
        const int *Y_local_sizes{nullptr}, *Y_local_offs{nullptr};
        switch (type) {
        case MatvecType::direct:
            X_size = n_cols;
            Y_size = n_rows;
            X_perm_ptr = &(this->_p_col_perm);
            Y_perm_ptr = &(this->_p_row_perm);

            X_local_sizes = local_col_sizes.data();
            Y_local_sizes = local_row_sizes.data();
            
            X_local_offs = local_col_offs.data();
            Y_local_offs = local_row_offs.data();
            
            break;
        case MatvecType::transp:
        case MatvecType::conj:
            X_size = n_rows;
            Y_size = n_cols;
            X_perm_ptr = &(this->_p_row_perm);
            Y_perm_ptr = &(this->_p_col_perm);

            X_local_sizes = local_row_sizes.data();
            Y_local_sizes = local_col_sizes.data();
            
            X_local_offs = local_row_offs.data();
            Y_local_offs = local_col_offs.data();
            break;
        }

        const std::unique_ptr<TPerm> &X_permutation{*X_perm_ptr};
        const std::unique_ptr<TPerm> &Y_permutation{*Y_perm_ptr};

        SmartArray<DataType> X_full;
#ifdef MSK_USE_MPI
        if (n_procs > 1) {
            X_full.resize(X_size * block_size);
            _gather_block_vector(block_size, X_size, X, X_full.data(), X_local_sizes, X_local_offs);
            X = X_full.data();
        }
#endif

        const DataType *X_in = X;
        SmartArray<DataType> X_perm;
        if (!pre_perm) {
            X_perm.resize(X_size * block_size);
            permute_block_vector(X_size, block_size, X_in, X_perm.data(), X_permutation, true);
            X_in = X_perm.data();
            X_full.free();
        }
        
        DataType *Y_out = Y;
        SmartArray<DataType> Y_full_perm;

        if (!pre_perm || n_procs > 1) { // For non-MPI case n_procs == 1
            Y_full_perm.resize(Y_size * block_size);
            Y_out = Y_full_perm.data();
        }

        constexpr DataType zero{0};
        BLAS::copy(Y_size * block_size, &zero, 0, Y_out, 1);

        std::function<BlockShape(const MatrixBlock &block, DataType *Y_local)> local_matvec;
        std::function<bool(const MatrixBlock &)> block_condition = [](const MatrixBlock &) { return true; };

        switch (msk_trans) {
        case Trans::gen:
            switch (type) {
            case MatvecType::direct:
                if (block_size == 1) {
                    local_matvec = [X_in](const MatrixBlock &block, DataType *y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        if (y_local != nullptr) {
                            block.matvec(X_in + j0, y_local);
                        }
                        return BlockShape{m, i0};
                    };
                } else {
                    local_matvec = [X_in, block_size, X_size, Y_size](const MatrixBlock &block, DataType *Y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        if (Y_local != nullptr) {
                            block.matvec(block_size, X_in + j0, X_size, Y_local, m);
                        }
                        return BlockShape{m, i0, block_size, Y_size};
                    };
                }
                break;
            case MatvecType::transp:
                if (block_size == 1) {
                    local_matvec = [X_in](const MatrixBlock &block, DataType *y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        if (y_local != nullptr) {
                            block.transp_matvec(X_in + i0, y_local);
                        }
                        return BlockShape{n, j0};
                    };
                } else {
                    local_matvec = [X_in, block_size, X_size, Y_size](const MatrixBlock &block, DataType *Y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        if (Y_local != nullptr) {
                            block.transp_matvec(block_size, X_in + i0, X_size, Y_local, n);
                        }
                        return BlockShape{n, j0, block_size, Y_size};
                    };
                }
                break;
            case MatvecType::conj:
                if (block_size == 1) {
                    local_matvec = [X_in](const MatrixBlock &block, DataType *y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        if (y_local != nullptr) {
                            block.conj_matvec(X_in + i0, y_local);
                        }
                        return BlockShape{n, j0};
                    };
                } else {
                    local_matvec = [X_in, block_size, X_size, Y_size](const MatrixBlock &block, DataType *Y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        if (Y_local != nullptr) {
                            block.conj_matvec(block_size, X_in + i0, X_size, Y_local, n);
                        }
                        return BlockShape{n, j0, block_size, Y_size};
                    };
                }
                break;
            }
            break;
        case Trans::herm:
            block_condition = [](const MatrixBlock &block) {
                uint64_t i0{}, j0{};
                std::tie(i0, j0) = block.corner();
                return i0 <= j0;
            };
            switch (type) {
            case MatvecType::direct:
            case MatvecType::conj:
                if (block_size == 1) {
                    local_matvec = [X_in](const MatrixBlock &block, DataType *y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        BlockShape result;
                        result.col_size = 1;
                        result.col_ldim = 0;
                        const bool is_diag{i0 == j0};
                        if (is_diag) { // block is diagonal
                            result.row_sizes_offs.resize(1);
                            result.row_sizes_offs[0] = {m, i0};
                        } else {
                            result.row_sizes_offs.resize(2);
                            result.row_sizes_offs[0] = {m, i0};
                            result.row_sizes_offs[1] = {n, j0};
                        }

                        if (y_local != nullptr) {
                            if (is_diag) {
                                block.matvec(X_in + j0, y_local);
                            } else {
                                block.matvec(X_in + j0, y_local);
                                block.conj_matvec(X_in + i0, y_local + m);
                            }
                        }
                        return result;
                    };
                } else {
                    local_matvec = [X_in, block_size, X_size, Y_size](const MatrixBlock &block, DataType *Y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        BlockShape result;
                        result.col_size = block_size;
                        result.col_ldim = Y_size;
                        const bool is_diag{i0 == j0};
                        if (is_diag) { // block is diagonal
                            result.row_sizes_offs.resize(1);
                            result.row_sizes_offs[0] = {m, i0};
                        } else {
                            result.row_sizes_offs.resize(2);
                            result.row_sizes_offs[0] = {m, i0};
                            result.row_sizes_offs[1] = {n, j0};
                        }
                        if (Y_local != nullptr) {
                            if (is_diag) {
                                block.matvec(block_size, X_in + j0, X_size, Y_local, m);
                            } else {
                                block.matvec(block_size, X_in + j0, X_size, Y_local, m);
                                block.conj_matvec(block_size, X_in + i0, X_size, Y_local + m * block_size, n);
                            }
                        }
                        return result;
                    };
                }
                break;
            case MatvecType::transp:
                if (block_size == 1) {
                    local_matvec = [X_in](const MatrixBlock &block, DataType *y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        BlockShape result;
                        result.col_size = 1;
                        result.col_ldim = 0;
                        const bool is_diag{i0 == j0};
                        if (is_diag) { // block is diagonal
                            result.row_sizes_offs.resize(1);
                            result.row_sizes_offs[0] = {m, i0};
                        } else {
                            result.row_sizes_offs.resize(2);
                            result.row_sizes_offs[0] = {m, i0};
                            result.row_sizes_offs[1] = {n, j0};
                        }

                        if (y_local != nullptr) {
                            if (is_diag) {
                                block.transp_matvec(X_in + j0, y_local);
                            } else {
                                block.compl_conj_matvec(X_in + j0, y_local);
                                block.transp_matvec(X_in + i0, y_local + m);
                            }
                        }
                        return result;
                    };
                } else {
                    local_matvec = [X_in, block_size, X_size, Y_size](const MatrixBlock &block, DataType *Y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        BlockShape result;
                        result.col_size = block_size;
                        result.col_ldim = Y_size;
                        const bool is_diag{i0 == j0};
                        if (is_diag) { // block is diagonal
                            result.row_sizes_offs.resize(1);
                            result.row_sizes_offs[0] = {m, i0};
                        } else {
                            result.row_sizes_offs.resize(2);
                            result.row_sizes_offs[0] = {m, i0};
                            result.row_sizes_offs[1] = {n, j0};
                        }
                        if (Y_local != nullptr) {
                            if (is_diag) {
                                block.transp_matvec(block_size, X_in + j0, X_size, Y_local, m);
                            } else {
                                block.compl_conj_matvec(block_size, X_in + j0, X_size, Y_local, m);
                                block.transp_matvec(block_size, X_in + i0, X_size, Y_local + m * block_size, n);
                            }
                        }
                        return result;
                    };
                }
                break;
            }
            break;
        case Trans::symm:
            block_condition = [](const MatrixBlock &block) {
                uint64_t i0{}, j0{};
                std::tie(i0, j0) = block.corner();
                return i0 <= j0;
            };
            switch (type) {
            case MatvecType::direct:
            case MatvecType::transp:
                if (block_size == 1) {
                    local_matvec = [X_in](const MatrixBlock &block, DataType *y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        BlockShape result;
                        result.col_size = 1;
                        result.col_ldim = 0;
                        const bool is_diag{i0 == j0};
                        if (is_diag) { // block is diagonal
                            result.row_sizes_offs.resize(1);
                            result.row_sizes_offs[0] = {m, i0};
                        } else {
                            result.row_sizes_offs.resize(2);
                            result.row_sizes_offs[0] = {m, i0};
                            result.row_sizes_offs[1] = {n, j0};
                        }
                        if (y_local != nullptr) {
                            if (is_diag) {
                                block.matvec(X_in + j0, y_local);
                            } else {
                                block.matvec(X_in + j0, y_local);
                                block.transp_matvec(X_in + i0, y_local + m);
                            }
                        }
                        return result;
                    };
                } else {
                    local_matvec = [X_in, block_size, X_size, Y_size](const MatrixBlock &block, DataType *Y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        BlockShape result;
                        result.col_size = block_size;
                        result.col_ldim = Y_size;
                        const bool is_diag{i0 == j0};
                        if (is_diag) { // block is diagonal
                            result.row_sizes_offs.resize(1);
                            result.row_sizes_offs[0] = {m, i0};
                        } else {
                            result.row_sizes_offs.resize(2);
                            result.row_sizes_offs[0] = {m, i0};
                            result.row_sizes_offs[1] = {n, j0};
                        }
                        if (Y_local != nullptr) {
                            if (is_diag) {
                                block.matvec(block_size, X_in + j0, X_size, Y_local, m);
                            } else {
                                block.matvec(block_size, X_in + j0, X_size, Y_local, m);
                                block.transp_matvec(block_size, X_in + i0, X_size, Y_local + m * block_size, n);
                            }
                        }
                        return result;
                    };
                }
                break;
            case MatvecType::conj:
                if (block_size == 1) {
                    local_matvec = [X_in](const MatrixBlock &block, DataType *y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        BlockShape result;
                        result.col_size = 1;
                        result.col_ldim = 0;
                        const bool is_diag{i0 == j0};
                        if (is_diag) { // block is diagonal
                            result.row_sizes_offs.resize(1);
                            result.row_sizes_offs[0] = {m, i0};
                        } else {
                            result.row_sizes_offs.resize(2);
                            result.row_sizes_offs[0] = {m, i0};
                            result.row_sizes_offs[1] = {n, j0};
                        }
                        if (y_local != nullptr) {
                            if (is_diag) {
                                block.conj_matvec(X_in + j0, y_local);
                            } else {
                                block.compl_conj_matvec(X_in + j0, y_local);
                                block.conj_matvec(X_in + i0, y_local + m);
                            }
                        }
                        return result;
                    };
                } else {
                    local_matvec = [X_in, block_size, X_size, Y_size](const MatrixBlock &block, DataType *Y_local) {
                        uint64_t m{}, n{}, i0{}, j0{};
                        std::tie(m, n) = block.size();
                        std::tie(i0, j0) = block.corner();
                        BlockShape result;
                        result.col_size = block_size;
                        result.col_ldim = Y_size;
                        const bool is_diag{i0 == j0};
                        if (is_diag) { // block is diagonal
                            result.row_sizes_offs.resize(1);
                            result.row_sizes_offs[0] = {m, i0};
                        } else {
                            result.row_sizes_offs.resize(2);
                            result.row_sizes_offs[0] = {m, i0};
                            result.row_sizes_offs[1] = {n, j0};
                        }
                        if (Y_local != nullptr) {
                            if (is_diag) {
                                block.conj_matvec(block_size, X_in + j0, X_size, Y_local, m);
                            } else {
                                block.compl_conj_matvec(block_size, X_in + j0, X_size, Y_local, m);
                                block.conj_matvec(block_size, X_in + i0, X_size, Y_local + m * block_size, n);
                            }
                        }
                        return result;
                    };
                }
            }
        
        }

        // Perform matvec
        apply_batched_block_linear(local_matvec, Y_out, block_condition);
        

        // Permute back
        if (n_procs == 1) { // For non-MPI case n_procs == 1
            if (!pre_perm) {
                permute_block_vector(Y_size, block_size, Y_out, Y, Y_permutation, false);
            }
            return;
        }

#ifdef MSK_USE_MPI
        // Y_out maps to a tmp array Y_full_perm
        SmartArray<DataType> Y_unperm;
        if (!pre_perm) {
            Y_unperm.resize(Y_size * block_size);
            permute_block_vector(Y_size, block_size, Y_out, Y_unperm.data(), Y_permutation, false);
            Y_out = Y_unperm.data();
            Y_full_perm.free();
        }
        MPI_Barrier(msk_comm);

        // Gather only needed part of y at each process
        _distreduce_block_vector(block_size, Y_size, Y_out, Y, Y_local_sizes, Y_local_offs);
#endif

    }

    static void set_default_split_rule(int n_procs, int size, SmartArray<int> &offs, SmartArray<int> &sizes) {
        offs.resize(n_procs + 1);
        sizes.resize(n_procs);
        const int p_common = size / n_procs;
        const int n_extra = size % n_procs;
        for (int p = 0; p < n_procs; ++p) {
            sizes[p] = p_common + (p < n_extra ? 1 : 0);
        }
        offs[0] = 0;
        for (uint64_t p = 0; p < n_procs; ++p) {
            offs[p + 1] = offs[p] + sizes[p];
        }        
    }

    static void set_split_rule(const int n_procs, const int global_size,
                               const int *offs, const int *sizes,
                               SmartArray<int> &local_offs, SmartArray<int> &local_sizes) {
        local_offs.resize(n_procs + 1);
        local_sizes.resize(n_procs);
        if (offs != nullptr) {
            for (uint64_t p = 0; p < n_procs; ++p) {
                local_offs[p] = offs[p];
            }
            local_offs[n_procs] = global_size;
            for (uint64_t p = 0; p < n_procs; ++p) {
                local_sizes[p] = local_offs[p + 1] - local_offs[p];
            }
        } else if (sizes != nullptr) {
            for (uint64_t p = 0; p < n_procs; ++p) {
                local_sizes[p] = sizes[p];
            }
            local_offs[0] = 0;
            for (uint64_t p = 0; p < n_procs; ++p) {
                local_offs[p + 1] = local_offs[p] + sizes[p];
            }            
        } else {
            set_default_split_rule(global_size, n_procs, local_offs, local_sizes);
        }
    }

    void refresh_stat() {
        // Calculate Frobenius norm of approximation
        f_norm_sq = 0;
        _memory = 0;
        for (auto block_it = _manager.begin(true); block_it != _manager.end(true); ++block_it) {
            _memory += block_it->memory();
            auto nrm = block_it->frob_norm();
            f_norm_sq += nrm * nrm;
        }

    }

    static uint64_t limit_workers(const uint64_t n_jobs, std::shared_ptr<log::Logger> log = nullptr) {
        if (n_jobs == 0) {
            if (log) {
                log->warn("No threads required, defaulting to single");
            } else {
                std::cerr << "No threads required, defaulting to single" << std::endl;
            }
            return 1;
        }
        const uint64_t system_jobs = std::thread::hardware_concurrency();
        const uint64_t max_jobs = system_jobs > 0 ? system_jobs : 1024;
        if (n_jobs > max_jobs) {
            if (log) {
                log->warn("Strange number of threads required: {}", n_jobs);
            } else {
                std::cerr << "Strange number of threads required: " << n_jobs << std::endl;
            }
            return max_jobs;
        }
        return n_jobs;
    }

    MosaicSkeleton(const uint64_t m, const uint64_t n, uint64_t n_jobs = 4,
#ifdef MSK_USE_MPI
                   MPI_Comm comm = MPI_COMM_WORLD,
#endif
                   std::shared_ptr<log::Logger> logger_ptr = nullptr) :
            n_rows{m}, n_cols{n},
            _p_row_perm{new index::Permutation<uint64_t>{n_rows}},
            _p_col_perm{new index::Permutation<uint64_t>{n_cols}},
            _manager{
#ifdef MSK_USE_MPI
                comm,
            #endif
                logger_ptr
            },
            p_pool{ new exec::ThreadPool{"MSk", n_jobs, n_jobs, 10}},
            #ifdef MSK_USE_MPI
            msk_comm{comm},
            #endif
            msk_log{std::move(logger_ptr)} {
        if (msk_log && proc_id == 0) msk_log->info("MSk: init on size {} x {}", n_rows, n_cols);

        // Set up multiprocessing data (MPI)
#ifdef MSK_USE_MPI
        MPI_Comm_rank(msk_comm, &proc_id);
        MPI_Comm_size(msk_comm, &n_procs);
#endif
        // Resolve local sizes for matvec
        set_default_split_rule(n_procs, n_rows, local_row_offs, local_row_sizes);
        set_default_split_rule(n_procs, n_cols, local_col_offs, local_col_sizes);

        // Start pool of threads
        p_pool->start();

#ifdef MSK_USE_MPI
        MPI_Barrier(msk_comm);
#endif
        if (msk_log && proc_id == 0) msk_log->debug("MSk: initialized");
    }



/********* Public methods *********/



public:
    explicit MosaicSkeleton() :
            n_rows{0}, n_cols{0}, _p_row_perm{nullptr}, _p_col_perm{nullptr},
            _manager{
#ifdef MSK_USE_MPI
                MPI_COMM_NULL,
#endif
                nullptr
            },
            p_pool{nullptr},
#ifdef MSK_USE_MPI
            msk_comm{MPI_COMM_NULL},
#endif
            msk_log{nullptr} {}

    /**
     * @brief Peform mosaic-skeleton approximation procedure.
     *
     * @tparam Matrix Functional matrix type, should provide element evaluation, for examples see msk/matrix/matrix.hpp.
     * @tparam Oracle Oracle type, should provide blocks division information, for examples see msk/oracle/oracle.hpp.
     *
     * @param[in] matrix Functional matrix to approximate.
     * @param[in] oracle Oracle entity to define mosaic (block hierarchical) structure.
     * @param[in] n_jobs Number of threads to use.
     * @if MSK_USE_MPI
     *   @param[in] comm MPI communicator to operate on.
     * @endif
     * @param[in] logger_ptr Logger.
     */
    template<class Matrix, class Oracle>
    static MosaicSkeleton approximate(Matrix &matrix, Oracle &oracle, uint64_t n_jobs = 4,
    #ifdef MSK_USE_MPI
                                      MPI_Comm comm = MPI_COMM_WORLD,
    #endif
                                      std::shared_ptr<log::Logger> logger_ptr = nullptr) {
        uint64_t m{}, n{};
        std::tie(m, n) = matrix.size();
        uint64_t n_workers{limit_workers(n_jobs, logger_ptr)};
#ifdef MSK_USE_MPI
        MosaicSkeleton msk(m, n, n_workers, comm, std::move(logger_ptr));
#else
        MosaicSkeleton msk(m, n, n_workers, std::move(logger_ptr));
#endif

        // Start approximation procedure
        // 1. Build mosaic blocks tree
        msk.build_mosaic(oracle);

        // 2. Perform blocks approximation
        msk.approximate(matrix);

        if (msk.msk_log && msk.proc_id == 0) msk.msk_log->info("MSk: approximated");
        return msk;
    }

    MosaicSkeleton(const MosaicSkeleton &) = delete;

    MosaicSkeleton &operator=(const MosaicSkeleton &) = delete;

    MosaicSkeleton(MosaicSkeleton &&other) noexcept :
            n_rows{other.n_rows}, n_cols{other.n_cols},
            _p_row_perm{std::move(other._p_row_perm)}, _p_col_perm{std::move(other._p_col_perm)},
            n_blocks{other.n_blocks}, _manager{std::move(other._manager)},
            f_norm_sq{other.f_norm_sq}, _memory{other._memory},
            p_pool{std::move(other.p_pool)},
#ifdef MSK_USE_MPI
            msk_comm{std::move(other.msk_comm)},
#endif
            local_row_offs{std::move(other.local_row_offs)}, local_col_offs{std::move(other.local_col_offs)},
            local_row_sizes{std::move(other.local_row_sizes)}, local_col_sizes{std::move(other.local_col_sizes)},
            proc_id{other.proc_id}, n_procs{other.n_procs},
            msk_log{std::move(other.msk_log)} {}

    MosaicSkeleton &operator=(MosaicSkeleton &&other) noexcept {
        if (this == &other) {
            return *this;
        }

        n_rows = other.n_rows, other.n_rows = 0;
        n_cols = other.n_cols, other.n_cols = 0;
        _p_row_perm = std::move(other._p_row_perm);
        _p_col_perm = std::move(other._p_col_perm);

        n_blocks = other.n_blocks, other.n_blocks = 0;
        _manager = std::move(other._manager);

        f_norm_sq = other.f_norm_sq, other.f_norm_sq = 0;
        _memory = other._memory, other._memory = 0;

        p_pool = std::move(other.p_pool);

#ifdef MSK_USE_MPI
        msk_comm = std::move(other.msk_comm), other.msk_comm = MPI_COMM_NULL;
#endif
        local_row_offs = std::move(other.local_row_offs);
        local_col_offs = std::move(other.local_col_offs);
        local_row_sizes = std::move(other.local_row_sizes);
        local_col_sizes = std::move(other.local_col_sizes);
        proc_id = other.proc_id, other.proc_id = 0;
        n_procs = other.n_procs, other.n_procs = 0;

        msk_log = std::move(other.msk_log);

        return *this;
    }

    ~MosaicSkeleton() {
        if (p_pool) {
            p_pool->stop(true);
        }
    }

    /**
     * @brief Perform matrix-vector product: y = A x.
     *
     * @details Mosaic-skeleton format permutes matrix's rows and columns.
     * For the better performance it is recommended to work with permuted vectors
     * and unpermute them at the end.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     * @param[in] pre_perm If the input vector is already permuted, output vector is treated the same way.
     */
    void matvec(const DataType *x, DataType *y, bool pre_perm = true) const {
        _general_matvec(MatvecType::direct, 1, x, y, pre_perm);
    }


     /**
     * @brief Perform block matrix-block-vector product: Y = A X.
     *
     * @details Mosaic-skeleton format permutes matrix's rows and columns.
     * For the better performance it is recommended to work with permuted vectors
     * and unpermute them at the end.
     *
     * @param[in] block_size Number of block columns.
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     * @param[in] pre_perm If the input vector is already permuted, output vector is treated the same way.
     */
    void matvec(const uint64_t block_size, const DataType *X, DataType *Y, bool pre_perm = true) const {
        _general_matvec(MatvecType::direct, block_size, X, Y, pre_perm);
    }

    /**
     * @brief Perform transposed matrix-vector product: y = A^T x.
     *
     * @details Mosaic-skeleton format permutes matrix's rows and columns.
     * For the better performance it is recommended to work with permuted vectors
     * and unpermute them at the end.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     * @param[in] pre_perm If the input vector is already permuted, output vector is treated the same way.
     */
    void transp_matvec(const DataType *x, DataType *y, bool pre_perm = true) const {
        _general_matvec(MatvecType::transp, 1, x, y, pre_perm);
    }

    /**
     * @brief Perform transposed matrix-block-vector product: Y = A^T X.
     *
     * @details Mosaic-skeleton format permutes matrix's rows and columns.
     * For the better performance it is recommended to work with permuted vectors
     * and unpermute them at the end.
     *
     * @param[in] block_size Number of block columns.
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     * @param[in] pre_perm If the input vector is already permuted, output vector is treated the same way.
     */
     void transp_matvec(const uint64_t block_size, const DataType *X, DataType *Y, bool pre_perm = true) const {
        _general_matvec(MatvecType::transp, block_size, X, Y, pre_perm);
    }

    /**
     * @brief Perform conjugated matrix-vector product: y = A^* x.
     *
     * @details Mosaic-skeleton format permutes matrix's rows and columns.
     * For the better performance it is recommended to work with permuted vectors
     * and unpermute them at the end.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     * @param[in] pre_perm If the input vector is already permuted, output vector is treated the same way.
     */
     void conj_matvec(const DataType *x, DataType *y, bool pre_perm = true) const {
        _general_matvec(MatvecType::conj, 1, x, y, pre_perm);
    }

    /**
     * @brief Perform conjugated matrix-block-vector product: Y = A^* X.
     *
     * @details Mosaic-skeleton format permutes matrix's rows and columns.
     * For the better performance it is recommended to work with permuted vectors
     * and unpermute them at the end.
     *
     * @param[in] block_size Number of block columns.
     * @param[in] X Input vector.
     * @param[out] Y Output vector.
     * @param[in] pre_perm If the input vector is already permuted, output vector is treated the same way.
     */
     void conj_matvec(const uint64_t block_size, const DataType *X, DataType *Y, bool pre_perm = true) const {
        _general_matvec(MatvecType::conj, block_size, X, Y, pre_perm);
    }

    /**
     * @brief Perform BLAS-style general matrix-vector product: y = op(A) x.
     *
     * @details Mosaic-skeleton format permutes matrix's rows and columns.
     * For the better performance it is recommended to work with permuted vectors
     * and unpermute them at the end.
     *
     * @param[in] trans Operation on the matrix:
     * if 'N' then op(A) = A, if 'T' then op(A) = A^T, if 'C' then op(A) = A^*.
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     * @param[in] pre_perm If the input vector is already permuted, output vector is treated the same way.
     */
    void matvec(const char trans, const DataType *x, DataType *y, bool pre_perm = true) const {
        if (trans == 'n' || trans == 'N') {
            matvec(x, y, pre_perm);
        } else if (trans == 'c' || trans == 'C') {
            conj_matvec(x, y, pre_perm);
        } else if (trans == 't' || trans == 'T') {
            transp_matvec(x, y, pre_perm);
        }
    }

    /**
     * @brief Perform BLAS-style general matrix-vector product: Y = op(A) X.
     *
     * @details Mosaic-skeleton format permutes matrix's rows and columns.
     * For the better performance it is recommended to work with permuted vectors
     * and unpermute them at the end.
     *
     * @param[in] trans Operation on the matrix:
     * if 'N' then op(A) = A, if 'T' then op(A) = A^T, if 'C' then op(A) = A^*.
     * @param[in] block_size Number of input/ourput vectors.
     * @param[in] X Input vectors.
     * @param[out] Y Output vectors.
     * @param[in] pre_perm If the input vector is already permuted, output vector is treated the same way.
     */
    void matvec(const char trans, const uint64_t block_size, const DataType *X, DataType *Y, bool pre_perm = true) const {
        if (trans == 'n' || trans == 'N') {
            matvec(block_size, X, Y, pre_perm);
        } else if (trans == 'c' || trans == 'C') {
            conj_matvec(block_size, X, Y, pre_perm);
        } else if (trans == 't' || trans == 'T') {
            transp_matvec(block_size, X, Y, pre_perm);
        }
    }

    /**
     * @brief Convolve two vectors using the matrix as the kernel.
     *
     * @details Convolution result contains first size entries of 'full' convolution.
     * So z_k = SUM[over j] y_j A_{j, k - j} x_{k - j} is computed for k < size (indexing starts from zero).
     * Parameter size=0 is also treated as 'full' of size M + N - 1 for A of size (M, N).
     * Mosaic partitioning is assumed to have identical rows and columns permutation
     * (no indices permutations while mosaic building).
     *
     * @note Currently works only for a single process version.
     *
     * @param[in] y Left-hand-side vector to convolve.
     * @param[in] x Right-hand-side vector to convolve.
     * @param[out] z Resulting convolution vector.
     * @param[in] size Required size of result.
     */
    void convolve(const DataType *y, const DataType *x, DataType *z, uint64_t size = 0) const {
        if (size == 0) {
            size = n_rows + n_cols - 1;
        }
        size = std::min(size, n_rows + n_cols - 1);
        constexpr DataType zero{0.0};
        BLAS::copy(size, &zero, 0, z, 1);

#ifdef MSK_USE_MPI
        // Gather x and y vectors from all processes
        SmartArray<DataType> x_gather(n_cols), y_gather(n_rows);
        _gather_block_vector(1, n_cols, x, x_gather.data(), local_col_sizes.data(), local_col_offs.data());
        _gather_block_vector(1, n_rows, y, y_gather.data(), local_row_sizes.data(), local_row_offs.data());
        const DataType *x_full = x_gather.data();
        const DataType *y_full = y_gather.data();
#else
        const DataType *x_full = x;
        const DataType *y_full = y;
#endif

        // Convolution procedure for a block

        std::function<BlockShape(const MatrixBlock &block, DataType *Y_local)> local_conv;
        std::function<bool(const MatrixBlock &)> block_condition = [](const MatrixBlock &) { return true; };
        switch (msk_trans) {
        case Trans::gen:
            local_conv = [&y_full, &x_full, &size](const MatrixBlock &block, DataType *z_local) {
                uint64_t m{}, n{}, i0{}, j0{};
                std::tie(m, n) = block.size();
                std::tie(i0, j0) = block.corner();
                uint64_t b_off{i0 + j0};
                uint64_t b_size{b_off < size ? size - b_off : 0};

                if (b_size > 0 && z_local != nullptr) {
                    block.convolve(y_full + i0, x_full + j0, z_local, b_size);
                }
                return BlockShape{b_size, b_off};
            };
            break;
        case Trans::symm:
            block_condition = [](const MatrixBlock &block) {
                uint64_t i0{}, j0{};
                std::tie(i0, j0) = block.corner();
                return i0 <= j0;
            };
            local_conv = [&y_full, &x_full, &size](const MatrixBlock &block, DataType *z_local) {
                uint64_t m{}, n{}, i0{}, j0{};
                std::tie(m, n) = block.size();
                std::tie(i0, j0) = block.corner();
                uint64_t b_off{i0 + j0};
                uint64_t b_size{b_off < size ? size - b_off : 0};
                if (b_size > 0 && z_local != nullptr) {
                    block.convolve(y_full + i0, x_full + j0, z_local, b_size);
                    if (i0 != j0) {
                        block.convolve(x_full + i0, y_full + j0, z_local, b_size);
                    }
                }
                return BlockShape{b_size, b_off};
            };
            break;
        case Trans::herm:
            block_condition = [](const MatrixBlock &block) {
                uint64_t i0{}, j0{};
                std::tie(i0, j0) = block.corner();
                return i0 <= j0;
            };
            local_conv = [&y_full, &x_full, &size](const MatrixBlock &block, DataType *z_local) {
                uint64_t m{}, n{}, i0{}, j0{};
                std::tie(m, n) = block.size();
                std::tie(i0, j0) = block.corner();
                uint64_t b_off{i0 + j0};
                uint64_t b_size{b_off < size ? size - b_off : 0};
                if (b_size > 0 && z_local != nullptr) {
                    block.convolve(y_full + i0, x_full + j0, z_local, b_size);
                    if (i0 != j0) {
                        block.conj_convolve(y_full + j0, x_full + i0, z_local, b_size);
                    }
                }
                return BlockShape{b_size, b_off};
            };
            break;
        };

        // Perform convolution
        apply_batched_block_linear(local_conv, z, block_condition);

#ifdef MSK_USE_MPI
        // Summarize convolution parts
        MPI::Allreduce_inplace_sum(z, size, msk_comm);
#endif
    }

    /**
     * @brief Save mosaic-skeleton matrix to a dense array.
     *
     * @details In multiprocess mode the result is saved only on a root process.
     * Non-root output arrays are not used.
     *
     * @param[out] ptr Array to save the matrix.
     * @param[in] ldim Leading dimension of resulting array.
     * @if MSK_USE_MPI
     *   @param[in] root Root process to save the result.
     * @endif
     */
#ifdef MSK_USE_MPI
    void to_dense(DataType *ptr, uint64_t ldim, int root = 0) const
#else
    void to_dense(DataType *ptr, uint64_t ldim) const
#endif
    {
#ifndef MSK_USE_MPI
        int root = 0;
#endif
        if (proc_id == root) {
            // Get all blocks in dense form
            SmartArray<DataType> plain_blocks(n_rows * n_cols);
            uint64_t off{0}, local_off = -1;
#ifdef MSK_USE_MPI
            // Gather blocks from other processes.
            SmartArray<MPI_Request> b_reqs(n_blocks);
            SmartArray<MPI_Status> b_stats(n_blocks);
            for (uint64_t b_id = 0; b_id < n_blocks; ++b_id) {
                auto &block = _manager.get_block(b_id);
                std::pair<uint64_t, uint64_t> b_size{block.size()};
                uint64_t b_elems{b_size.first * b_size.second};
                if (_manager.block_is_local(b_id)) {
                    if (local_off > off) {
                        local_off = off;
                    }
                    b_reqs[b_id] = MPI_REQUEST_NULL;
                } else {
                    uint8_t *recv_ptr = reinterpret_cast<uint8_t *>(plain_blocks.data() + off);
                    uint64_t count{b_elems * sizeof(DataType)};
                    MPI_Irecv(recv_ptr, count, MPI_UINT8_T, MPI_ANY_SOURCE, b_id, msk_comm, &b_reqs[b_id]);
                }
                off += b_elems;
            }
#else
            local_off = 0;
#endif
            off = local_off;
            // Densify blocks on the current process
            for (auto block_it = _manager.begin(true); block_it != _manager.end(true); ++block_it) {
                std::pair<uint64_t, uint64_t> b_size{block_it->size()};
                uint64_t b_elems{b_size.first * b_size.second};
                block_it->to_dense(plain_blocks.data() + off, b_size.first);
                off += b_elems;
            }
#ifdef MSK_USE_MPI
            MPI_Waitall(n_blocks, b_reqs.data(), b_stats.data());
#endif
            DataType *block_ptr{plain_blocks.data()};
            for (auto block_it = _manager.begin(false); block_it != _manager.end(false); ++block_it) {
                auto b_idx = block_it->get_idx();
                std::pair<uint64_t, uint64_t> b_size{block_it->size()};
                uint64_t b_elems{b_size.first * b_size.second};
                uint64_t j{0};
                // Put each block's elemets on their permuted places
                for (auto jdx = b_idx.second.begin(); jdx != b_idx.second.end(); ++jdx, ++j) {
                    uint64_t i{0};
                    for (auto idx = b_idx.first.begin(); idx != b_idx.first.end(); ++idx, ++i) {
                        ptr[*idx + *jdx * ldim] = block_ptr[i + b_size.first * j];
                    }
                }
                block_ptr += b_elems;
            }
        } else {
#ifdef MSK_USE_MPI
            // At non-root process
            // Get total size of dense blocks
            uint64_t total_elems{0};
            for (auto block_it = _manager.begin(true); block_it != _manager.end(true); ++block_it) {
                std::pair<uint64_t, uint64_t> b_size{block_it->size()};
                total_elems += b_size.first * b_size.second;
            }
            SmartArray<DataType> local_blocks(total_elems);
            DataType *put_ptr{local_blocks.data()};
            uint64_t local_n_blocks{_manager.size(true)};
            SmartArray<MPI_Request> send_reqs(local_n_blocks);
            SmartArray<MPI_Status> send_stats(local_n_blocks);
            MPI_Request *send_req_ptr{send_reqs.data()};
            // Densify each block and send it ti the root process
            for (auto block_it = _manager.begin(true); block_it != _manager.end(true); ++block_it) {
                std::pair<uint64_t, uint64_t> b_size{block_it->size()};
                uint64_t b_elems{b_size.first * b_size.second};
                block_it->to_dense(put_ptr, b_size.first);
                uint8_t *send_ptr = reinterpret_cast<uint8_t *>(put_ptr);
                uint64_t count{b_elems * sizeof(DataType)};
                MPI_Isend(send_ptr, count, MPI_UINT8_T, root, block_it->id(), msk_comm, send_req_ptr);
                ++send_req_ptr;
                put_ptr += b_elems;
            }
            MPI_Waitall(local_n_blocks, send_reqs.data(), send_stats.data());
#endif
        }
    }

#ifdef MSK_USE_MPI
    void permute_col(const DataType *x_in, DataType *x_out, bool global = true) const {
        if (!global || n_procs == 1) {
            _p_col_perm->permute(x_in, x_out);
            return;
        }
        SmartArray<DataType> x_in_full(n_cols), x_out_full(n_cols);
        // Gather x from all processes
        _gather_block_vector(1, n_cols, x_in, x_in_full.data(), local_col_sizes.data(), local_col_offs.data());
        _p_col_perm->permute(x_in_full.data(), x_out_full.data());
        BLAS::copy(local_col_sizes[proc_id], x_out_full.data() + local_col_offs[proc_id], 1, x_out, 1);
    }

    inline void permute_col(DataType *x_in_out, bool global = true) const {
        if (!global || n_procs == 1) {
            _p_col_perm->permute(x_in_out);
            return;
        }
        SmartArray<DataType> x_in(local_col_sizes[proc_id]);
        BLAS::copy(local_col_sizes[proc_id], x_in_out, 1, x_in.data(), 1);
        permute_col(x_in.data(), x_in_out, global);
    }

    void unpermute_col(const DataType *x_in, DataType *x_out, bool global = true) const {
        if (!global || n_procs == 1) {
            _p_col_perm->unpermute(x_in, x_out);
            return;
        }
        SmartArray<DataType> x_in_full(n_cols), x_out_full(n_cols);
        // Gather x from all processes
        _gather_block_vector(1, n_cols, x_in, x_in_full.data(), local_col_sizes.data(), local_col_offs.data());
        _p_col_perm->unpermute(x_in_full.data(), x_out_full.data());
        BLAS::copy(local_col_sizes[proc_id], x_out_full.data() + local_col_offs[proc_id], 1, x_out, 1);
    }

    inline void unpermute_col(DataType *x_in_out, bool global = true) const {
        if (!global || n_procs == 1) {
            _p_col_perm->unpermute(x_in_out);
            return;
        }
        SmartArray<DataType> x_in(local_col_sizes[proc_id]);
        BLAS::copy(local_col_sizes[proc_id], x_in_out, 1, x_in.data(), 1);
        unpermute_col(x_in.data(), x_in_out, global);
    }
#else
    /**
     * @brief Permute vector according to matrix columns permutation.
     *
     * @param[in] x_in Input vector to permute.
     * @param[out] x_out Permuted output vector.
     * @if MSK_USE_MPI
     *   @param[in] global If the vector is split by processes.
     * @endif
     */
    inline void permute_col(const DataType *x_in, DataType *x_out) const {
        _p_col_perm->permute(x_in, x_out);
    }

    /**
     * @brief Permute vector according to matrix columns permutation inplace.
     *
     * @param[in, out] x_in_out Vector to permute.
     * @if MSK_USE_MPI
     *   @param[in] global If the vector is split by processes.
     * @endif
     */
    void permute_col(DataType *x_in_out) const {
        _p_col_perm->permute(x_in_out);
    }

    /**
     * @brief Unpermute (inverse permute) vector according to matrix columns permutation.
     *
     * @param[in] x_in Input vector to unpermute.
     * @param[out] x_out Unpermuted output vector.
     * @if MSK_USE_MPI
     *   @param[in] global If the vector is split by processes.
     * @endif
     */
    inline void unpermute_col(const DataType *x_in, DataType *x_out) const {
        _p_col_perm->unpermute(x_in, x_out);
    }

    /**
     * @brief Unpermute (inverse permute) vector according to matrix columns permutation inplace.
     *
     * @param[in, out] x_in_out Vector to unpermute.
     * @if MSK_USE_MPI
     *   @param[in] global If the vector is split by processes.
     * @endif
     */
    void unpermute_col(DataType *x_in_out) const {
        _p_col_perm->unpermute(x_in_out);
    }
#endif

#ifdef MSK_USE_MPI
    void permute_row(const DataType *y_in, DataType *y_out, bool global = true) const {
        if (!global || n_procs == 1) {
            _p_row_perm->permute(y_in, y_out);
            return;
        }
        SmartArray<DataType> y_in_full(n_rows), y_out_full(n_rows);
        // Gather y from all processes
        _gather_block_vector(1, n_rows, y_in, y_in_full.data(), local_row_sizes.data(), local_row_offs.data());
        _p_row_perm->permute(y_in_full.data(), y_out_full.data());
        BLAS::copy(local_row_sizes[proc_id], y_out_full.data() + local_row_offs[proc_id], 1, y_out, 1);
    }

    inline void permute_row(DataType *y_in_out, bool global = true) const {
        if (!global || n_procs == 1) {
            _p_row_perm->permute(y_in_out);
            return;
        }
        SmartArray<DataType> y_in(local_row_sizes[proc_id]);
        BLAS::copy(local_row_sizes[proc_id], y_in_out, 1, y_in.data(), 1);
        permute_row(y_in.data(), y_in_out, global);
    }

    void unpermute_row(const DataType *y_in, DataType *y_out, bool global = true) const {
        if (!global || n_procs == 1) {
            _p_row_perm->unpermute(y_in, y_out);
            return;
        }
        SmartArray<DataType> y_in_full(n_rows), y_out_full(n_rows);
        // Gather y from all processes
        _gather_block_vector(1, n_rows, y_in, y_in_full.data(), local_row_sizes.data(), local_row_offs.data());
        _p_row_perm->unpermute(y_in_full.data(), y_out_full.data());
        BLAS::copy(local_row_sizes[proc_id], y_out_full.data() + local_row_offs[proc_id], 1, y_out, 1);
    }

    inline void unpermute_row(DataType *y_in_out, bool global = true) const {
        if (!global || n_procs == 1) {
            _p_row_perm->unpermute(y_in_out);
            return;
        }
        SmartArray<DataType> y_in(local_row_sizes[proc_id]);
        BLAS::copy(local_row_sizes[proc_id], y_in_out, 1, y_in.data(), 1);
        unpermute_row(y_in.data(), y_in_out, global);
    }
#else
    /**
     * @brief Permute vector according to matrix rows permutation.
     *
     * @param[in] y_in Input vector to permute.
     * @param[out] y_out Permuted output vector.
     * @if MSK_USE_MPI
     *   @param[in] global If the vector is split by processes.
     * @endif
     */
    void permute_row(const DataType *y_in, DataType *y_out) const {
        _p_row_perm->permute(y_in, y_out);
    }

    /**
     * @brief Permute vector according to matrix rows permutation inplace.
     *
     * @param[in, out] y_in_out Vector to permute.
     * @if MSK_USE_MPI
     *   @param[in] global If the vector is split by processes.
     * @endif
     */
    void permute_row(DataType *y_in_out) const {
        _p_row_perm->permute(y_in_out);
    }

    /**
     * @brief Unpermute (inverse permute) vector according to matrix rows permutation.
     *
     * @param[in] y_in Input vector to unpermute.
     * @param[out] y_out Unpermuted output vector.
     * @if MSK_USE_MPI
     *   @param[in] global If the vector is split by processes.
     * @endif
     */
    void unpermute_row(const DataType *y_in, DataType *y_out) const {
        _p_row_perm->unpermute(y_in, y_out);
    }

    /**
     * @brief Unpermute (inverse permute) vector according to matrix rows permutation inplace.
     *
     * @param[in, out] y_in_out Vector to unpermute.
     * @if MSK_USE_MPI
     *   @param[in] global If the vector is split by processes.
     * @endif
     */
    void unpermute_row(DataType *y_in_out) const {
        _p_row_perm->unpermute(y_in_out);
    }
#endif

    /**
     * @brief Get Frobenius norm af the matrix.
     *
     * @if MSK_USE_MPI
     *   @param[in] global If to return the whole matrix norm or for stored at the current process only.
     * @endif
     * @return Norm value.
     */
#ifdef MSK_USE_MPI
    RealType frob_norm(bool global = true) const {
        if (!global) {
            return std::sqrt(f_norm_sq);
        }
        RealType result{f_norm_sq};
        MPI::Allreduce_inplace_sum(&result, 1, msk_comm);
        return std::sqrt(result);
    }
#else
    RealType frob_norm() const {
        return std::sqrt(f_norm_sq);
    }
#endif

    /**
     * @brief Get the memory amount used to store the representation.
     *
     * @if MSK_USE_MPI
     *   @param[in] global If to return the whole matrix memory or for stored at the current process only.
     * @endif
     * @return Number of DataType memory cells.
     */
#ifdef MSK_USE_MPI
    uint64_t memory(bool global = true) const {
        uint64_t result = _memory;
        if (global) {
            MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_UINT64_T, MPI_SUM, msk_comm);
        }
        return result;
    }
#else
    uint64_t memory() const {
        return _memory;
    }
#endif

    /**
     * @brief Perform recompression procedure.
     *
     * @param[in] rel_tol Error relative tolerance.
     */
    void recompress(double rel_tol) {
        uint64_t compressible_mem{0};
        for (auto block_it = _manager.begin(true); block_it != _manager.end(true); ++block_it) {
            if(!block_it->exact()) {
                uint64_t m, n;
                std::tie(m, n) = block_it->size();
                compressible_mem += m * n;
            }
        }
#ifdef MSK_USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &compressible_mem, 1, MPI_UINT64_T, MPI_SUM, msk_comm);
#endif
        double total_norm{frob_norm()};
        auto b_recompress = [rel_tol, compressible_mem, total_norm](MatrixBlock &block, DataType *) {
            uint64_t m{}, n{};
            std::tie(m, n) = block.size();
            double part{(1.0 * m * n) / compressible_mem};
            // Each block's is set according to the part of the matrix it takes
            comp::Parameters params{rel_tol * std::sqrt(part) * total_norm, 0, 0};
            block.recompress(params);
            if (block.memory() >= m * n) {
                block.densify();
            }
            return BlockShape{0, 0}; // dummy
        };

        // Perform recompression
        apply_batched_block_linear(b_recompress);

#ifdef MSK_USE_MPI
        MPI_Barrier(msk_comm);
        // Rebalance blocks with the new memory distribution
        _manager.balance_blocks();
#endif
        refresh_stat();
        /*
        uint64_t max_rank{0};
        for (auto block_it = _manager.begin(true); block_it != _manager.end(true); ++block_it) {
            if(!block_it->exact()) {
                uint64_t m, n;
                std::tie(m, n) = block_it->size();
                max_rank = std::max(max_rank, block_it->memory() / (m + n));
            }
        }
        if (msk_log && proc_id == 0) msk_log->info("Maximal rank = {}", max_rank);
        */
    }

    /**
     * @brief Calculate true approximation error.
     *
     * @details Error is calculated blockwise and then accumulated.
     *
     * @tparam Matrix Comparing matrix type.
     * @param[in] matrix Matrix to compare to.
     * @param[in] relative If to calculate relative error or absolute.
     * @if MSK_USE_MPI
     *   @param[in] global If to return the whole matrix error or for stored at the current process only.
     * @endif
     * @return Error value.
     */
#ifdef MSK_USE_MPI
    template<class Matrix>
    double check_error(Matrix &matrix, bool relative = true, bool global = true) const {
#else
    template<class Matrix>
    double check_error(Matrix &matrix, bool relative = true) const {
#endif
        auto b_check_error = [&matrix, relative](const MatrixBlock &block, DataType *err_nrm) {
            if (err_nrm != nullptr) {
                double err{}, nrm{};
                if (block.exact()) {
                    err = 0;
                    if (relative) {
                        nrm = block.frob_norm();
                    }
                } else {
                    uint64_t m, n;
                    std::tie(m, n) = block.size();
                    SmartArray<DataType> approx_block(m * n), true_block(m * n);
                    block.to_dense(approx_block.data(), m);
                    matrix.block(block.rows(), block.columns(), true_block.data());
                    BLAS::axpy(m * n, DataType{-1.0}, true_block.data(), 1, approx_block.data(), 1);
                    err = BLAS::nrm2(m * n, approx_block.data(), 1);
                    if (relative) {
                        nrm = BLAS::nrm2(m * n, true_block.data(), 1);
                    }
                }
                err_nrm[0] = err * err;
                if (relative) {
                    err_nrm[1] = nrm * nrm;
                }
            }
            return BlockShape{relative ? 2ul : 1ul, 0};
        };

        // Accumulate error, and norm if needed
        DataType err_nrm_sq[2]{0, 0};
        apply_batched_block_linear(b_check_error, &err_nrm_sq[0]);

        double local_nrm{0}, local_err{0};
        local_err = std::real(err_nrm_sq[0]);
        if (relative) {
            local_nrm = std::real(err_nrm_sq[1]);;
        }

#ifdef  MSK_USE_MPI
        if (global) {
            MPI_Allreduce(MPI_IN_PLACE, &local_err, 1, MPI_DOUBLE, MPI_SUM, msk_comm);
            if (relative) {
                MPI_Allreduce(MPI_IN_PLACE, &local_nrm, 1, MPI_DOUBLE, MPI_SUM, msk_comm);
            }
        }
#endif
        local_err = std::sqrt(local_err);
        if (relative) {
            local_err /= std::sqrt(local_nrm);
        }
        return local_err;
    }

    /**
     * @brief Calculate true approximation error in max (Chebyshev) norm.
     *
     * @details Error is calculated blockwise and then accumulated.
     *
     * @tparam Matrix Comparing matrix type.
     * @param[in] matrix Matrix to compare to.
     * @param[in] relative If to calculate relative error or absolute.
     * @if MSK_USE_MPI
     *   @param[in] global If to return the whole matrix error or for stored at the current process only.
     * @endif
     * @return Error value.
     */
#ifdef MSK_USE_MPI
    template<class Matrix>
    double check_max_error(Matrix &matrix, bool relative = true, bool global = true) const {
#else
    template<class Matrix>
    double check_max_error(Matrix &matrix, bool relative = true) const {
#endif
        RealType local_err{0}, local_nrm{0};
        for (auto it = _manager.begin(true); it < _manager.end(true); ++it) {
            const auto &block = *it;
            if (!block.exact()) {
                uint64_t m, n;
                std::tie(m, n) = block.size();
                SmartArray<DataType> approx_block(m * n), true_block(m * n);
                block.to_dense(approx_block.data(), m);
                matrix.block(block.rows(), block.columns(), true_block.data());
                BLAS::axpy(m * n, DataType{-1.0}, true_block.data(), 1, approx_block.data(), 1);
                for (uint64_t i = 0; i < m * n; ++i) {
                    local_err = std::max(local_err, std::abs(approx_block[i]));
                }
                if (relative) {
                    for (uint64_t i = 0; i < m * n; ++i) {
                        local_nrm = std::max(local_nrm, std::abs(true_block[i]));
                    }
                }
            } else if (relative) {
                uint64_t m, n;
                std::tie(m, n) = block.size();
                SmartArray<DataType> true_block(m * n);
                matrix.block(block.rows(), block.columns(), true_block.data());
                for (uint64_t i = 0; i < m * n; ++i) {
                    local_nrm = std::max(local_nrm, std::abs(true_block[i]));
                }
            }
        }

#ifdef  MSK_USE_MPI
        if (global) {
            MPI_Allreduce(MPI_IN_PLACE, &local_err, 1, MPI_DOUBLE, MPI_MAX, msk_comm);
            if (relative) {
                MPI_Allreduce(MPI_IN_PLACE, &local_nrm, 1, MPI_DOUBLE, MPI_MAX, msk_comm);
            }
        }
#endif
        if (relative) {
            local_err /= local_nrm;
        }
        return local_err;
    }


    /**
     * @brief Calculate true approximation error in max (Chebyshev) norm for separate blocks.
     *
     * @details Error is calculated blockwise and gathered among processes.
     *
     * @tparam Matrix Comparing matrix type.
     * @param[in] matrix Matrix to compare to.
     * @param[out] errors Vector of errors with respect to blocks ids.
     * @param[out] norms Vector of norms with respect to blocks ids.
     */
    template<class Matrix>
    void check_max_error_blockwise(Matrix &matrix, std::vector<double> &errors, std::vector<double> &norms) const {
        errors.resize(n_blocks);
        norms.resize(n_blocks);
        for (auto it = _manager.begin(true); it < _manager.end(true); ++it) {
            RealType local_err{0}, local_nrm{0};
            const auto &block = *it;uint64_t m, n;
            std::tie(m, n) = block.size();
            SmartArray<DataType> true_block(m * n);
            matrix.block(block.rows(), block.columns(), true_block.data());
            for (uint64_t i = 0; i < m * n; ++i) {
                local_nrm = std::max(local_nrm, std::abs(true_block[i]));
            }
            if (!block.exact()) {
                SmartArray<DataType> approx_block(m * n);
                block.to_dense(approx_block.data(), m);
                BLAS::axpy(m * n, DataType{-1.0}, true_block.data(), 1, approx_block.data(), 1);
                for (uint64_t i = 0; i < m * n; ++i) {
                    local_err = std::max(local_err, std::abs(approx_block[i]));
                }
            }
            errors[block.id()] = local_err;
            norms[block.id()] = local_nrm;
        }
#ifdef  MSK_USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, errors.data(), n_blocks, MPI_DOUBLE, MPI_MAX, msk_comm);
        MPI_Allreduce(MPI_IN_PLACE, norms.data(), n_blocks, MPI_DOUBLE, MPI_MAX, msk_comm);
#endif
    }

    /**
     * @brief Write blocks information to a file.
     *
     * @details Format: "ID M N I0 J0 [K]\n" for each block, where Id is the unique id of a block,
     * M and N are its row and column sizes, I0 and J0 are the indices of its first row and columns
     * and (optionally) K is the amount of memory it uses.
     *
     * @param[in] filename File name to write to.
     * @param[in] memory If to write memory of blocks or not.
     */
    void dump_blocks(const char *filename, bool memory = false) {
        SmartArray<uint64_t> mems;
        if (memory) {
            mems.resize(n_blocks);
            uint64_t zero{0};
            BLAS::copy(n_blocks, &zero, 0, mems.data(), 1);
            for (auto block_it = _manager.begin(true); block_it != _manager.end(true); ++block_it) {
                mems[block_it->id()] = block_it->memory();
            }
            // Gather memory information from other processes
#ifdef MSK_USE_MPI
            if (proc_id == 0) {
                MPI_Reduce(MPI_IN_PLACE, mems.data(), n_blocks, MPI_UINT64_T, MPI_SUM, 0, msk_comm);
            } else {
                MPI_Reduce(mems.data(), mems.data(), n_blocks, MPI_UINT64_T, MPI_SUM, 0, msk_comm);
            }
#endif
        }
        std::vector<int> flags(n_blocks, 0);
        for (auto block_it = _manager.begin(true); block_it != _manager.end(true); ++block_it) {
            flags[block_it->id()] = block_it->exact();
        }
#ifdef MSK_USE_MPI
        if (proc_id == 0) {
            MPI_Reduce(MPI_IN_PLACE, flags.data(), n_blocks, MPI_INT, MPI_SUM, 0, msk_comm);
        } else {
            MPI_Reduce(flags.data(), flags.data(), n_blocks, MPI_INT, MPI_SUM, 0, msk_comm);
        }
#endif
        if (proc_id == 0) {
            std::ofstream f(filename);
            uint64_t block_id{0};
            for (auto block_it = _manager.begin(false); block_it != _manager.end(false); ++block_it, ++block_id) {
                uint64_t m, n, i0, j0;
                std::tie(m, n) = block_it->size();
                std::tie(i0, j0) = block_it->corner();
                f << block_id << " " << m << " " << n << " " << i0 << " " << j0;
                if (memory) {
                    f << " " << mems[block_id];
                }
                f << " " << flags[block_id];
                f << "\n";
            }
            f.close();
        }
    }

    /**
     * @brief Extract all blocks' data in low-rank format.
     *
     * @details If not all blocks are low-rank then corresponding pointers are null.
     *
     * @param[out] blocks Resulting vector of pointers to low-rank blocks.
     */
    void extract_blocks(std::vector<const block::LowRank<DataType> *> &blocks) {
        blocks.resize(n_blocks);
        using BType = std::tuple<uint64_t, uint64_t, const block::LowRank<DataType> *>;
        std::vector<BType> unsorted_blocks;
        unsorted_blocks.reserve(n_blocks);
        for (auto block_it = _manager.begin(false); block_it != _manager.end(false); ++block_it) {
            uint64_t /*m, n,*/ i0, j0;
            // std::tie(m, n) = block_it->size();
            std::tie(i0, j0) = block_it->corner();
            auto ptr = dynamic_cast<const block::LowRank<DataType> *>(block_it->data());
            unsorted_blocks.emplace_back(i0, j0, ptr);
        }
        std::sort(unsorted_blocks.begin(), unsorted_blocks.end());
        for (uint64_t i = 0; i < n_blocks; ++i) {
            blocks[i] = std::get<2>(unsorted_blocks[i]);
        }
    }

    /**
     * @brief Extract all blocks.
     *
     * @param[out] blocks Resulting vector of pointers to blocks.
     */
    void blocks(std::vector<block::Block<DataType> *> &blocks) {
        blocks.resize(n_blocks);
        using BType = std::tuple<uint64_t, uint64_t, block::Block<DataType> *>;
        std::vector<BType> unsorted_blocks;
        unsorted_blocks.reserve(n_blocks);
        for (auto block_it = _manager.begin(false); block_it != _manager.end(false); ++block_it) {
            uint64_t /*m, n,*/ i0, j0;
            // std::tie(m, n) = block_it->size();
            std::tie(i0, j0) = block_it->corner();
            unsorted_blocks.emplace_back(i0, j0, &(*block_it));
        }
        std::sort(unsorted_blocks.begin(), unsorted_blocks.end());
        for (uint64_t i = 0; i < n_blocks; ++i) {
            blocks[i] = std::get<2>(unsorted_blocks[i]);
        }
    }

    decltype(_manager.begin(true)) begin(bool local = true) const {
        return _manager.begin(local);
    }

    decltype(_manager.end(true)) end(bool local = true) const {
        return _manager.end(local);
    }

    /**
     * @brief Write rows permutation array to a file.
     *
     * @param[in] f File object to write to.
     */
    void dump_row_perm(std::ostream &f) const {
        for (uint64_t i = 0; i < n_rows; ++i) {
            f << (*_p_row_perm)[i] << " ";
        }
        f << "\n";
    }

    /**
     * @brief Write columns permutation array to a file.
     *
     * @param[in] f File object to write to.
     */
    void dump_col_perm(std::ostream &f) const {
        for (uint64_t i = 0; i < n_cols; ++i) {
            f << (*_p_col_perm)[i] << " ";
        }
        f << "\n";
    }

    /**
     * @brief Get the matrix size (shape).
     *
     * @return Rows number, columns number.
     */
    std::pair<uint64_t, uint64_t> size() const {
        return {n_rows, n_cols};
    }

    /**
     * @brief Get rows related offset index for the current process.
     *
     * @details Routine for multiprocessing purposes, otherwise returns 0.
     * It is the starting index of an output vector at matvec procedure and
     * the starting index of an input vector at transp_matvec procedure.
     *
     * @return Rows offset.
     */
    uint64_t local_row_off() const {
        return local_row_offs[proc_id];
    }

    /**
     * @brief Get columns related offset index for the current process.
     *
     * @details Routine for multiprocessing purposes, otherwise returns 0.
     * It is the starting index of an input vector at matvec procedure and
     * the starting index of an output vector at transp_matvec procedure.
     *
     * @return Cols offset.
     */
    uint64_t local_col_off() const {
        return local_col_offs[proc_id];
    }

    /**
     * @brief Get the number of rows belonging to the current process.
     *
     * @details Routine for multiprocessing purposes, otherwise returns matrix rows number.
     * It is the size of an output vector at matvec procedure and
     * the size of an input vector at transp_matvec procedure.
     *
     * @return Number of local rows.
     */
    uint64_t local_row_size() const {
        return local_row_sizes[proc_id];
    }

    /**
     * @brief Get the number of columns belonging to the current process.
     *
     * @details Routine for multiprocessing purposes, otherwise returns matrix columns number.
     * It is the size of an input vector at matvec procedure and
     * the size of an output vector at transp_matvec procedure.
     *
     * @return Number of local rows.
     */
    uint64_t local_col_size() const {
        return local_col_sizes[proc_id];
    }

    /**
     * @brief Get matrix local rows number, see local_row_size() for details.
     *
     * @return Number of local rows.
     */
    uint64_t get_row_number() const {
        return local_row_size();
    }

    /**
     * @brief Get matrix local columns number, see local_col_size() for details.
     *
     * @return Number of local columns.
     */
    uint64_t get_column_number() const {
        return local_col_size();
    }

    /**
     * @brief Get matrix number of rows.
     *
     * @return Number of rows.
     */
    uint64_t get_global_row_number() const {
        return n_rows;
    }

    /**
     * @brief Get matrix number of columns.
     *
     * @return Number of columns.
     */
    uint64_t get_global_column_number() const {
        return n_cols;
    }

#ifdef MSK_USE_MPI
    /**
     * @if MSK_USE_MPI
     *   @brief Get MPI communicator of the matrix.
     * @endif
     */
    MPI_Comm get_comm() const {
        return msk_comm;
    }
#endif

#ifdef MSK_USE_MPI
    /**
     * @if MSK_USE_MPI
     *   @brief Get MPI communicator of the matrix.
     * @endif
     */
    MPI_Comm comm() const {
        return msk_comm;
    }
#endif

#ifdef MSK_USE_MPI
    void set_col_split_rule(const int *offs, const int *sizes = nullptr) {
        set_split_rule(n_procs, n_cols, offs, sizes, local_col_sizes, local_col_offs);
    }

    void set_row_split_rule(const int *offs, const int *sizes = nullptr) {
        set_split_rule(n_procs, n_rows, offs, sizes, local_row_sizes, local_row_offs);
    }
#endif

    /**
     * @brief Get the logger.
     *
     * @return Logger instance.
     */
    std::shared_ptr<log::Logger> logger() {
        return msk_log;
    }
};


extern template class MosaicSkeleton<float>;

extern template class MosaicSkeleton<double>;

extern template class MosaicSkeleton<std::complex<float>>;

extern template class MosaicSkeleton<std::complex<double>>;


extern template MosaicSkeleton<float> MosaicSkeleton<float>::approximate(matrix::Elementary<float> &matrix,
        oracle::Elementary<float> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<float> MosaicSkeleton<float>::approximate(matrix::Dense<float, 2> &matrix,
        oracle::Distant<float, matrix::Dense<float, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<float> MosaicSkeleton<float>::approximate(matrix::Dense<float, 3> &matrix,
        oracle::Distant<float, matrix::Dense<float, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<float> MosaicSkeleton<float>::approximate(matrix::BiPoint<float, float, 2> &matrix,
        oracle::Distant<float, matrix::BiPoint<float, float, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<float> MosaicSkeleton<float>::approximate(matrix::BiPoint<float, float, 3> &matrix,
        oracle::Distant<float, matrix::BiPoint<float, float, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);


extern template MosaicSkeleton<double> MosaicSkeleton<double>::approximate(matrix::Elementary<double> &matrix,
        oracle::Elementary<double> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<double> MosaicSkeleton<double>::approximate(matrix::Dense<double, 2> &matrix,
        oracle::Distant<double, matrix::Dense<double, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<double> MosaicSkeleton<double>::approximate(matrix::Dense<double, 3> &matrix,
        oracle::Distant<double, matrix::Dense<double, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<double> MosaicSkeleton<double>::approximate(matrix::BiPoint<double, double, 2> &matrix,
        oracle::Distant<double, matrix::BiPoint<double, double, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<double> MosaicSkeleton<double>::approximate(matrix::BiPoint<double, double, 3> &matrix,
        oracle::Distant<double, matrix::BiPoint<double, double, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);


extern template MosaicSkeleton<std::complex<float>> MosaicSkeleton<std::complex<float>>::approximate(matrix::Elementary<std::complex<float>> &matrix,
        oracle::Elementary<std::complex<float>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<std::complex<float>> MosaicSkeleton<std::complex<float>>::approximate(matrix::Dense<std::complex<float>, 2> &matrix,
        oracle::Distant<std::complex<float>, matrix::Dense<std::complex<float>, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<std::complex<float>> MosaicSkeleton<std::complex<float>>::approximate(matrix::Dense<std::complex<float>, 3> &matrix,
        oracle::Distant<std::complex<float>, matrix::Dense<std::complex<float>, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<std::complex<float>> MosaicSkeleton<std::complex<float>>::approximate(matrix::BiPoint<std::complex<float>, float, 2> &matrix,
        oracle::Distant<std::complex<float>, matrix::BiPoint<std::complex<float>, float, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<std::complex<float>> MosaicSkeleton<std::complex<float>>::approximate(matrix::BiPoint<std::complex<float>, float, 3> &matrix,
        oracle::Distant<std::complex<float>, matrix::BiPoint<std::complex<float>, float, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);


extern template MosaicSkeleton<std::complex<double>> MosaicSkeleton<std::complex<double>>::approximate(matrix::Elementary<std::complex<double>> &matrix,
        oracle::Elementary<std::complex<double>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<std::complex<double>> MosaicSkeleton<std::complex<double>>::approximate(matrix::Dense<std::complex<double>, 2> &matrix,
        oracle::Distant<std::complex<double>, matrix::Dense<std::complex<double>, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<std::complex<double>> MosaicSkeleton<std::complex<double>>::approximate(matrix::Dense<std::complex<double>, 3> &matrix,
        oracle::Distant<std::complex<double>, matrix::Dense<std::complex<double>, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<std::complex<double>> MosaicSkeleton<std::complex<double>>::approximate(matrix::BiPoint<std::complex<double>, double, 2> &matrix,
        oracle::Distant<std::complex<double>, matrix::BiPoint<std::complex<double>, double, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);

extern template MosaicSkeleton<std::complex<double>> MosaicSkeleton<std::complex<double>>::approximate(matrix::BiPoint<std::complex<double>, double, 3> &matrix,
        oracle::Distant<std::complex<double>, matrix::BiPoint<std::complex<double>, double, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<log::Logger> logger_ptr);


extern template void MosaicSkeleton<float>::build_mosaic(oracle::Elementary<float> &oracle);
extern template void MosaicSkeleton<float>::build_mosaic(oracle::Distant<float, matrix::Dense<float, 2>> &oracle);
extern template void MosaicSkeleton<float>::build_mosaic(oracle::Distant<float, matrix::Dense<float, 3>> &oracle);
extern template void MosaicSkeleton<float>::build_mosaic(oracle::Distant<float, matrix::BiPoint<float, float, 2>> &oracle);
extern template void MosaicSkeleton<float>::build_mosaic(oracle::Distant<float, matrix::BiPoint<float, float, 3>> &oracle);

extern template void MosaicSkeleton<float>::approximate(matrix::Elementary<float> &matrix);
extern template void MosaicSkeleton<float>::approximate(matrix::Dense<float, 2> &matrix);
extern template void MosaicSkeleton<float>::approximate(matrix::Dense<float, 3> &matrix);
extern template void MosaicSkeleton<float>::approximate(matrix::BiPoint<float, float, 2> &matrix);
extern template void MosaicSkeleton<float>::approximate(matrix::BiPoint<float, float, 3> &matrix);


extern template void MosaicSkeleton<double>::build_mosaic(oracle::Elementary<double> &oracle);
extern template void MosaicSkeleton<double>::build_mosaic(oracle::Distant<double, matrix::Dense<double, 2>> &oracle);
extern template void MosaicSkeleton<double>::build_mosaic(oracle::Distant<double, matrix::Dense<double, 3>> &oracle);
extern template void MosaicSkeleton<double>::build_mosaic(oracle::Distant<double, matrix::BiPoint<double, double, 2>> &oracle);
extern template void MosaicSkeleton<double>::build_mosaic(oracle::Distant<double, matrix::BiPoint<double, double, 3>> &oracle);

extern template void MosaicSkeleton<double>::approximate(matrix::Elementary<double> &matrix);
extern template void MosaicSkeleton<double>::approximate(matrix::Dense<double, 2> &matrix);
extern template void MosaicSkeleton<double>::approximate(matrix::Dense<double, 3> &matrix);
extern template void MosaicSkeleton<double>::approximate(matrix::BiPoint<double, double, 2> &matrix);
extern template void MosaicSkeleton<double>::approximate(matrix::BiPoint<double, double, 3> &matrix);


extern template void MosaicSkeleton<std::complex<float>>::build_mosaic(oracle::Elementary<std::complex<float>> &oracle);
extern template void MosaicSkeleton<std::complex<float>>::build_mosaic(oracle::Distant<std::complex<float>, matrix::Dense<std::complex<float>, 2>> &oracle);
extern template void MosaicSkeleton<std::complex<float>>::build_mosaic(oracle::Distant<std::complex<float>, matrix::Dense<std::complex<float>, 3>> &oracle);
extern template void MosaicSkeleton<std::complex<float>>::build_mosaic(oracle::Distant<std::complex<float>, matrix::BiPoint<std::complex<float>, float, 2>> &oracle);
extern template void MosaicSkeleton<std::complex<float>>::build_mosaic(oracle::Distant<std::complex<float>, matrix::BiPoint<std::complex<float>, float, 3>> &oracle);

extern template void MosaicSkeleton<std::complex<float>>::approximate(matrix::Elementary<std::complex<float>> &matrix);
extern template void MosaicSkeleton<std::complex<float>>::approximate(matrix::Dense<std::complex<float>, 2> &matrix);
extern template void MosaicSkeleton<std::complex<float>>::approximate(matrix::Dense<std::complex<float>, 3> &matrix);
extern template void MosaicSkeleton<std::complex<float>>::approximate(matrix::BiPoint<std::complex<float>, float, 2> &matrix);
extern template void MosaicSkeleton<std::complex<float>>::approximate(matrix::BiPoint<std::complex<float>, float, 3> &matrix);


extern template void MosaicSkeleton<std::complex<double>>::build_mosaic(oracle::Elementary<std::complex<double>> &oracle);
extern template void MosaicSkeleton<std::complex<double>>::build_mosaic(oracle::Distant<std::complex<double>, matrix::Dense<std::complex<double>, 2>> &oracle);
extern template void MosaicSkeleton<std::complex<double>>::build_mosaic(oracle::Distant<std::complex<double>, matrix::Dense<std::complex<double>, 3>> &oracle);
extern template void MosaicSkeleton<std::complex<double>>::build_mosaic(oracle::Distant<std::complex<double>, matrix::BiPoint<std::complex<double>, double, 2>> &oracle);
extern template void MosaicSkeleton<std::complex<double>>::build_mosaic(oracle::Distant<std::complex<double>, matrix::BiPoint<std::complex<double>, double, 3>> &oracle);

extern template void MosaicSkeleton<std::complex<double>>::approximate(matrix::Elementary<std::complex<double>> &matrix);
extern template void MosaicSkeleton<std::complex<double>>::approximate(matrix::Dense<std::complex<double>, 2> &matrix);
extern template void MosaicSkeleton<std::complex<double>>::approximate(matrix::Dense<std::complex<double>, 3> &matrix);
extern template void MosaicSkeleton<std::complex<double>>::approximate(matrix::BiPoint<std::complex<double>, double, 2> &matrix);
extern template void MosaicSkeleton<std::complex<double>>::approximate(matrix::BiPoint<std::complex<double>, double, 3> &matrix);


}

#endif //MOSAIC_SKELETON_MOSAIC_SKELETON_HPP


