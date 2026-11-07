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
// Created by bulat_v on 06.09.22.
//


/**
 * @file
 *
 * @brief Contains class block::Manager.
 */


#ifndef MOSAIC_SKELETON_BLOCK_MANAGER_HPP
#define MOSAIC_SKELETON_BLOCK_MANAGER_HPP


#include "msk_config.h"

#include <cstdint>
#include <exception>
#include <mutex>
#include <condition_variable>
#include <map>
#include <set>

#include "smart_array/smart_array.hpp"
#ifdef MSK_USE_MPI
#include "cpp_mpi/MPI.hpp"
#endif

#include "msk/computation/computer.hpp"
#include "msk/block/block.hpp"
#include "msk/thread_pool/thread_pool.hpp"
#include "msk/utility/logging.hpp"


namespace MSk {

namespace eval {
template<typename DataType>
class Queue;
}

namespace block {


/**
 * @brief Manager of a set of block::Block entities.
 *
 * @details Manages evaluation procedure: queueing tasks and requesting elements,
 * provides iteration over the set and some synchronization.
 *
 * @tparam DataType Blocks elements' type
 */
template<class DataType>
class Manager {
private:
    /**
     * @brief Structure with block's statistics for sorting reasons.
     */
    struct BSortInfo {
        uint64_t id{}; ///< Unique identificator.
        int old_proc_id{0}; ///< Current (MPI) process id.
        int new_proc_id{0}; ///< New (MPI) process id to migrate to.
        uint64_t mem{0}; ///< Memory required to store the block.
        uint64_t prior{0}; ///< Priority flag of the block.
    public:
        BSortInfo(uint64_t b_id = 0, int old_h = 0, int new_h = 0, uint64_t memory = 0, uint64_t priority = 0) :
                id{b_id}, old_proc_id{old_h}, new_proc_id{new_h}, mem{memory}, prior{priority} {}
    };

    /**
     * @brief Status of Request of a Block.
     */
    struct RequestInfo {
        bool to_step_not_eval{true}; ///< Flag of performing "step" requirement.
        uint64_t enqueued{0}; ///< Number of elements left enqueued for evaluation.
        eval::Request<DataType> request{}; ///< Request for matrix elements.
    };

    /// Synchronization tools for approximation process.
    std::map<uint64_t, RequestInfo> block_requested; ///< Request of each block.
    SmartArray<uint64_t> _ready_blocks; ///< Blocks' ids whose requests need to be checked.
    uint64_t _put_done_idx{0}; ///< Queue of ids' to be checked.
    std::mutex _m; ///< The queue lock.
    std::condition_variable cv; ///< Condition variable to control the queue's non-emptiness.
    std::exception_ptr approximation_error{}; ///< First failed block step, reported after workers finish.

    /// Set of blocks stored on the current (MPI) process.
    using TBlock = Block<DataType>;
    uint64_t n_blocks{};
    std::vector<TBlock> _blocks; ///< The blocks themselves.
    uint64_t local_n_blocks{}, local_start_b_id{}, local_end_b_id{}; ///< Global ids of locally stored blocks.

    /// Multiprocessing (MPI) data
    int n_procs{1}, proc_id{0};
#ifdef MSK_USE_MPI
    MPI_Comm _comm;
    SmartArray<MPI_Comm> block_comms;
    SmartArray<int> block_procs; ///< MPI ranks of blocks.
#endif

    std::shared_ptr<log::Logger> log{nullptr};


public:
    Manager(
#ifdef MSK_USE_MPI
            MPI_Comm comm,
#endif
            std::shared_ptr<log::Logger> log_ptr = nullptr) :
#ifdef MSK_USE_MPI
            _comm{comm},
#endif
            log{std::move(log_ptr)} {
#ifdef MSK_USE_MPI
        if (_comm != MPI_COMM_NULL) {
            MPI_Comm_size(_comm, &n_procs);
            MPI_Comm_rank(_comm, &proc_id);
        }
#endif
    }

    ~Manager() {
#ifdef MSK_USE_MPI
        int flag{0};
        MPI_Finalized(&flag);
        if (!flag) {
            for (uint64_t i = 0; i < n_blocks; ++i) {
                MPI_Comm_free(&block_comms[i]);
            }
        }
#endif
    }

    Manager(Manager &&other) noexcept :
            n_blocks{other.n_blocks}, _blocks{std::move(other._blocks)}, local_n_blocks{other.local_n_blocks},
            local_start_b_id{other.local_start_b_id}, local_end_b_id{other.local_end_b_id},
            n_procs{other.n_procs}, proc_id{other.proc_id},
#ifdef MSK_USE_MPI
            _comm{std::move(other._comm)}, block_comms{std::move(other.block_comms)},
            block_procs{std::move(other.block_procs)},
#endif
            log{std::move(other.log)} {
        other.n_blocks = other.local_n_blocks = 0;
        other.local_start_b_id = other.local_end_b_id = 0;
        other._blocks.clear();
    }

    Manager &operator=(Manager &&other) noexcept {
        if (this == &other) {
            return *this;
        }
        n_blocks = other.n_blocks, other.n_blocks = 0;
        local_n_blocks = other.local_n_blocks, other.local_n_blocks = 0;
        _blocks = std::move(other._blocks), other._blocks.clear();

        local_start_b_id = other.local_start_b_id, other.local_start_b_id = 0;
        local_end_b_id = other.local_end_b_id, other.local_end_b_id = 0;

        n_procs = other.n_procs, other.n_procs = 0;
        proc_id = other.proc_id, other.proc_id = 0;
#ifdef MSK_USE_MPI
        _comm = std::move(other._comm), other._comm = MPI_COMM_NULL;
        block_comms = std::move(other.block_comms);
        block_procs = std::move(other.block_procs);
#endif

        log = std::move(other.log);

        return *this;
    }

    /**
     * @brief Transfer blocks to Manager.
     *
     * Should be called once before the computation.
     * Procedure requires all blocks at every process.
     * They may initially be ordered differently but sorting is performed uniformly.
     *
     * @param[in, out] blocks Blocks to manage.
     */
    void init_blocks(std::vector<TBlock> &&blocks) {
        _blocks = std::move(blocks);
        _blocks.shrink_to_fit();
        n_blocks = _blocks.size();
        std::function<bool(const TBlock &, const TBlock &)> comparator =
                [](const TBlock &lhs, const TBlock &rhs) {
                    return lhs.is_col_less(rhs);
                };
        // Sort blocks the same way at every process
        sort(comparator);
        local_start_b_id = 0;
        local_end_b_id = local_n_blocks = n_blocks;
#ifdef MSK_USE_MPI
        block_comms.resize(n_blocks);
        for (uint64_t i = 0; i < n_blocks; ++i) {
            MPI_Comm_dup(_comm, &block_comms[i]);
            _blocks[i].set_comm(block_comms[i]);
        }
        if (proc_id != 0) {
            local_start_b_id = local_end_b_id = local_n_blocks = 0;
        }
        block_procs.resize(n_blocks);
        // Balance blocks between processes
        balance_cold_start(false, true);
#endif
    }

    /**
     * @brief Number of blocks.
     *
     * @param[in] local Flag of getting the number of locally stored blocks or in total.
     * @return Number of blocks
     */
    uint64_t size(bool local) const {
        if (local) {
            return local_n_blocks;
        } else {
            return n_blocks;
        }
    }

    /**
     * @brief Beginning iterator over the set of blocks.
     *
     * @param[in] local Flag: iterate over locally stored blocks only or the whole set.
     * @return Iterator.
     */
    typename decltype(_blocks)::iterator begin(bool local) {
        if (local) {
            return _blocks.begin() + local_start_b_id;
        } else {
            return _blocks.begin();
        }
    }

    /**
     * @brief Ending iterator over the set of blocks.
     *
     * @param[in] local Flag: iterate over locally stored blocks only or the whole set.
     * @return Iterator.
     */
    typename decltype(_blocks)::iterator end(bool local) {
        if (local) {
            return _blocks.begin() + local_end_b_id;
        } else {
            return _blocks.end();
        }
    }

    /**
     * @brief Constant beginning iterator over the set of blocks.
     *
     * @param[in] local Flag: iterate over locally stored blocks only or the whole set.
     * @return Const Iterator.
     */
    typename decltype(_blocks)::const_iterator begin(bool local) const {
        if (local) {
            return _blocks.cbegin() + local_start_b_id;
        } else {
            return _blocks.cbegin();
        }
    }

    /**
     * @brief Constant ending iterator over the set of blocks.
     *
     * @param[in] local Flag: iterate over locally stored blocks only or the whole set.
     * @return Const Iterator.
     */
    typename decltype(_blocks)::const_iterator end(bool local) const {
        if (local) {
            return _blocks.cbegin() + local_end_b_id;
        } else {
            return _blocks.cend();
        }
    }

    /**
     * @brief Get the block in global set.
     *
     * @param[in] idx Global id of the block.
     * @return Block.
     */
    const TBlock &get_block(uint64_t idx) const {
        if (idx >= n_blocks) {
            throw std::logic_error("Block index out of range");
        }
        return _blocks[idx];
    }

    /**
     * @brief Check if block is stored locally.
     *
     * @param[in] block_id Global block id.
     * @return Boolean locality flag.
     */
    bool block_is_local(uint64_t block_id) const {
        return local_start_b_id <= block_id && block_id < local_end_b_id;
    }

#ifdef MSK_USE_MPI
    /**
     * @brief Balance blocks between processes.
     *
     * @param[in] warm Use blocks' memory information or sizes and position only.
     */
    void balance_blocks(bool warm = true) {
        if (warm) {
            balance_warm_start();
        } else {
            balance_cold_start();
        }
    }
#else
    void balance_blocks() {
        balance_cold_start(true, false);
    }
#endif

    /**
     * @brief Sort the whole set of blocks.
     *
     * @param[in] comparator Comparator used for sorting: "less" blocks will be palced ahead.
     *
     * @note All local blocks should go sequentianlly.
     */
    void sort(const std::function<bool(const TBlock &, const TBlock &)> &comparator) {
        std::sort(_blocks.data(), _blocks.data() + n_blocks, comparator);
        for (uint64_t b_id = 0; b_id < n_blocks; ++b_id) {
            _blocks[b_id].set_id(b_id);
        }
#ifdef MSK_USE_MPI
        MPI_Barrier(_comm);
#endif
    }

    /**
     * @brief Start blocks' approximation (evaluation) procedure.
     *
     * @param[in] pool ThreadPool for tasks (approximation steps).
     * @param[in] elements_queue Queue for requesting elements of blocks (matrix).
     */
    void start_approximation(exec::ThreadPool *pool, eval::Queue<DataType> *elements_queue) {
        approximation_error = nullptr;
        // Every block is considered to have only one request to be looked at the moment
        _ready_blocks.resize(2 * local_n_blocks + 1);
        // Firstly, every block is ging to make a step
        for (uint64_t i = local_start_b_id; i < local_end_b_id; ++i) {
            block_requested[i].to_step_not_eval = true;
            block_requested[i].enqueued = 0;
            block_requested[i].request = {};
            _ready_blocks[i - local_start_b_id] = i;
        }

        uint64_t _look_done_idx{0};
        _put_done_idx = local_n_blocks;
        uint64_t next_blocks_verbose{local_n_blocks};
#ifdef MSK_USE_MPI
        if (log) log->info("{} :: start approximation of {} blocks, id = [{}, {})", proc_id,
                           local_n_blocks, local_start_b_id, local_end_b_id);
#else
        if (log) log->info("Start approximation of {} blocks, id = [{}, {})",
                           local_n_blocks, local_start_b_id, local_end_b_id);
#endif
        uint64_t blocks_left{local_n_blocks}, iteration{0};
        std::vector<eval::Request<DataType>> req_to_enq(local_n_blocks);
        while (blocks_left > 0) {
            // Log at every n_local_blocks/5 blocks completion
            if (blocks_left <= next_blocks_verbose) {
                if (log) log->debug("Iter {}, blocks left: {}", iteration, blocks_left);
                next_blocks_verbose -= std::min(local_n_blocks / 5, next_blocks_verbose);
            }
            // Retrieve ids of blocks whose requests need to be checked
            std::unique_lock<std::mutex> lock(_m);
            uint64_t n_look;
            if (_look_done_idx <= _put_done_idx) {
                n_look = _put_done_idx - _look_done_idx;
            } else {
                n_look = _ready_blocks.size() + _put_done_idx - _look_done_idx;
            }
            if (n_look == 0) { // Sleep if no blocks to check
                cv.wait(lock);
                continue;
            }
            // n_look > 0: process blocks withoul locking the queue of ready ones
            lock.unlock();
            uint64_t to_enq{0};
            for (uint64_t i = 0; i < n_look; ++i) {
                uint64_t block_id{_ready_blocks[_look_done_idx]};
                auto &req_info = block_requested.at(block_id);
                if (req_info.to_step_not_eval) { // Set a block to make a step
                    pool->execute(&Manager::block_step, this, block_id);
                } else if (!req_info.request.empty()) { // Elements to evaluate and no step required -> request
                    auto &request = req_info.request;
                    req_info.enqueued += request.size();
                    req_to_enq[to_enq] = std::move(request); // Gather all current requests of elements
                    ++to_enq;
                } else { // Approximation of a Block is finished
                    --blocks_left;
                }
                ++_look_done_idx;
                if (_look_done_idx == _ready_blocks.size()) {
                    _look_done_idx = 0;
                }
            }
            if (to_enq > 0 && !elements_queue->add(to_enq, req_to_enq.data())) {
                if (log) log->error("Could not add elements to evaluation queue");
            }
            to_enq = 0;
            ++iteration;
        }
        _ready_blocks.free();
        if (approximation_error) std::rethrow_exception(approximation_error);
#ifdef MSK_USE_MPI
        if (log) log->info("{} :: manager finished approximation", proc_id);
#else
        if (log) log->debug("Manager finished approximation");
#endif
    }

    /**
     * @brief Mark requests as complete and signalize corresponding blocks.
     *
     * @param[in] n Number of requests.
     * @param[in] eval_reqs Requests array.
     */
    void transfer_elements(const uint64_t n, const eval::Request<DataType> *eval_reqs) {
        for (uint64_t i = 0; i < n; ++i) {
            auto &req = eval_reqs[i];
            uint64_t block_id = req.block_id();
            auto &req_info = block_requested.at(block_id);
            req_info.enqueued -= req.size();
            if (log) log->debug("Transfer to block {}: {} elems, {} left", block_id, req.size(), req_info.enqueued);
            if (req_info.enqueued == 0) {
                mark_task_done(block_id);
            }
        }
    }


private:
    /**
     * @brief Perform one step of a block approximation procedure.
     * 
     * @param[in] block_id Global block id.
     */
    void block_step(uint64_t block_id) {
        if (log) log->debug("Step on block {}", block_id);
        try {
            comp::Result res = _blocks[block_id].step_enqueue(std::ref(block_requested.at(block_id).request));
            if (log && res.finished) log->debug("Finished block {}", block_id);
        } catch (...) {
            // A failed worker must still release its slot in the manager queue;
            // otherwise the manager waits forever for a completion it will not see.
            block_requested.at(block_id).request.free();
            std::lock_guard<std::mutex> lock(_m);
            if (!approximation_error) approximation_error = std::current_exception();
        }
        mark_task_done(block_id);
    }

    /**
     * @brief Mark block's request as complete, change its state and put to the queue to be checked.
     * 
     * @param[in] block_id Global block id.
     */
    void mark_task_done(uint64_t block_id) {
        block_requested.at(block_id).to_step_not_eval = !block_requested.at(block_id).to_step_not_eval;
        {
            std::unique_lock<std::mutex> lock(_m);
            _ready_blocks[_put_done_idx] = block_id;
            ++_put_done_idx;
            if (_put_done_idx == _ready_blocks.size()) {
                _put_done_idx = 0;
            }
            lock.unlock();
            cv.notify_one();
        }
    }

    /**
     * @brief Count statistics all blocks.
     *
     * @param[out] info Array of stat structures.
     * @param[in] use_mem Include info about memory utilization.
     * @param[in] exact_prior Prioritize exact blocks.
     */
    void count_block_costs(BSortInfo *info, bool use_mem = true, bool exact_prior = true) {
        for (uint64_t i = 0; i < n_blocks; ++i) {
            info[i] = {0, 0, 0, 0};
        }
        uint64_t off = local_start_b_id;
        for (auto block_it = begin(true); block_it != end(true); ++block_it, ++off) {
            uint64_t priority{0};
            if (exact_prior && block_it->exact()) {
                priority = 1;
            }
            uint64_t mem{0};
            if (use_mem) {
                mem = block_it->memory();
            } else {
                uint64_t m, n;
                std::tie(m, n) = block_it->size();
                mem = m * n;
            }
            info[off] = {block_it->id(), proc_id, proc_id, mem, priority};
        }
#ifdef MSK_USE_MPI
        // Values for non-local are expected to zeroes for plain addition
        MPI_Allreduce(MPI_IN_PLACE, info, sizeof(info[0]) * n_blocks, MPI_UINT8_T, MPI_SUM, _comm);
#endif
    }


// #ifdef MSK_USE_MPI
    /**
     * @brief Balance blocks between processes ignoring current distribution.
     *
     * @details First, all processes are set as empty boxes. Then all blocks
     * starting from the largest (and prioritized) are assigned to the least filled box.
     *
     * @param[in] use_mem Use information about memory utilization.
     * @param[in] exact_prior Prioritize exact blocks.
     */
    void balance_cold_start(bool use_mem = true, bool exact_prior = true) {
        SmartArray<BSortInfo> block_info(n_blocks);
        count_block_costs(block_info.data(), use_mem, exact_prior);
        auto mem_comp = [](const BSortInfo &a, const BSortInfo &b) {
            return std::tie(a.prior, a.mem) > std::tie(b.prior, b.mem);
        };
        // Sort all blocks according to priority and memory utilization (if corresponding flag is true)
        std::stable_sort(block_info.data(), block_info.data() + n_blocks, mem_comp);
        std::vector<uint64_t> proc_fills(n_procs, 0);
        for (uint64_t old_b_id = 0; old_b_id < n_blocks; ++old_b_id) {
            uint64_t min_proc{0};
            for (int p = 0; p < n_procs; ++p) {
                if (proc_fills[p] < proc_fills[min_proc]) {
                    min_proc = p;
                }
            }
            // Assign largest block to the least filled process
            proc_fills[min_proc] += block_info[old_b_id].mem;
            block_info[old_b_id].new_proc_id = min_proc;
        }
        // Now redistribute blocks
        redistribute_blocks(block_info);
    }
// #endif

// Iteratively rebalance blocks, trying to get procs to mean capacity, two procs exchange at step
#ifdef MSK_USE_MPI
    /**
     * @brief Balance blocks between processes with respect to the current distribution.
     *
     * @details Each process counts its imbalance = |1 - M P / G|,
     * where M is local memory, G is total memory, P is the number of processes.
     * While some process's imbalance is above the tolerance, we choose the least and
     * the most filled processes. Then we find the block, whose transfer would reduce
     * both imbalances closer to the mean. The number of such iterations is limited to 1000.
     *
     * @param[in] eps Parameter of relative acceptable imbalance.
     */
    void balance_warm_start(double eps = 0.005) {
        if (eps < 0) {
            return;
        }
        SmartArray<BSortInfo> block_info(n_blocks);
        count_block_costs(block_info.data());

        /// TODO: get global memory from BSortInfo array
        uint64_t global_mem{};
        for (auto block_it = begin(true); block_it != end(true); ++block_it) {
            global_mem += block_it->memory();
        }
        MPI_Allreduce(MPI_IN_PLACE, &global_mem, 1, MPI_UINT64_T, MPI_SUM, _comm);

        // Count starting processes filling
        // Sets of blocks at aech process at the moment
        std::vector<std::set<std::pair<uint64_t, uint64_t>>> local_blocks(n_procs); // {mem, _id}
        std::vector<uint64_t> proc_fills(n_procs, 0);
        for (uint64_t b_id{0}; b_id < n_blocks; ++b_id) {
            local_blocks[block_info[b_id].old_proc_id].emplace(block_info[b_id].mem, b_id);
            proc_fills[block_info[b_id].old_proc_id] += block_info[b_id].mem;
        }

        // Count starting maximal imbalance
        double max_eps{0};
        for (int p{0}; p < n_procs; ++p) {
            double p_eps{(1.0 * proc_fills[p] * n_procs) / global_mem - 1.0};
            max_eps = std::max(max_eps, p_eps);
        }

        // Main iteration loop
        uint64_t iter{0}, max_iter{1000};
        while (max_eps > eps && iter < max_iter) {
            ++iter;
            // Find least and most filled processes
            uint64_t mem_min{proc_fills[0]}, mem_max{proc_fills[0]};
            int p_min{0}, p_max{0};
            for (int p{0}; p < n_procs; ++p) {
                if (proc_fills[p] > mem_max) {
                    p_max = p;
                    mem_max = proc_fills[p];
                }
                if (proc_fills[p] < mem_min) {
                    p_min = p;
                    mem_min = proc_fills[p];
                }
            }
            if (p_min == p_max) {
                break;
            }

            // Mean imbalance
            auto send_pos{local_blocks[p_max].end()};
            double best_eps{0};
            for (uint64_t p : {p_max, p_min}) {
                double p_eps{std::abs((1.0 * proc_fills[p] * n_procs) / global_mem - 1.0)};
                best_eps = std::max(best_eps, p_eps);
            }
            // Find the best block to transfer
            for (auto to_send = local_blocks[p_max].begin(); to_send != local_blocks[p_max].end(); ++to_send) {
                uint64_t mem_at_max{mem_max - to_send->first}, mem_at_min{mem_min + to_send->first};
                double eps_at_max{(1.0 * mem_at_max * n_procs) / global_mem - 1.0};
                double eps_at_min{(1.0 * mem_at_min * n_procs) / global_mem - 1.0};
                if (eps_at_max < -best_eps || eps_at_min > best_eps) { // too big block, makes balance worse
                    break;
                }
                double new_eps{std::max(std::abs(eps_at_min), std::abs(eps_at_max))};
                if (new_eps > best_eps) { // Further block would be worse
                    break;
                }
                send_pos = to_send;
                best_eps = new_eps;
                // "Little" process' filling is above the mean and "big"'s is below the mean
                if (eps_at_min > 0 || eps_at_max < 0) {
                    break;
                }
            }
            if (send_pos == local_blocks[p_max].end() || best_eps >= max_eps) { // no good block to send
                break;
            }
            // Apply chosen transfer (without actual data sending)
            proc_fills[p_max] -= send_pos->first;
            proc_fills[p_min] += send_pos->first;
            local_blocks[p_min].insert(*send_pos);
            assert(block_info[send_pos->second].old_proc_id == p_max && "wrong block sending");
            block_info[send_pos->second].old_proc_id = block_info[send_pos->second].new_proc_id;
            block_info[send_pos->second].new_proc_id = p_min;
            local_blocks[p_max].erase(send_pos);
            // Update current maximal imbalance
            double max_new_eps{0};
            for (uint64_t p : {p_max, p_min}) {
                double p_eps{std::abs((1.0 * proc_fills[p] * n_procs) / global_mem - 1.0)};
                max_new_eps = std::max(max_new_eps, p_eps);
            }
            max_eps = std::min(max_eps, max_new_eps);
        }
        // Perform the actual blocks' data transfer
        redistribute_blocks(block_info);
    }
#endif

// #ifdef MSK_USE_MPI
    /**
     * @brief Redistribute blocks according to a new partition.
     *
     * @param[in, out] block_info Array of transfer info.
     */
    void redistribute_blocks(SmartArray<BSortInfo> &block_info) {
        // Sort according new partition so blocks at every process will lie sequentially
        auto comp = [](const BSortInfo &a, const BSortInfo &b) {
            return std::tie(a.new_proc_id, b.mem, a.old_proc_id, a.id) < 
                   std::tie(b.new_proc_id, a.mem, b.old_proc_id, b.id);
        };
        std::stable_sort(block_info.data(), block_info.data() + n_blocks, comp);

        std::vector<TBlock> new_blocks(n_blocks);
        for (uint64_t new_id = 0; new_id < n_blocks; ++new_id) {
            uint64_t old_id{block_info[new_id].id};
            new_blocks[new_id] = std::move(_blocks[old_id]);
        }
        _blocks = std::move(new_blocks);

#ifdef MSK_USE_MPI
        // Actual blocks transfer
        uint64_t n_send{0};
#pragma omp parallel for default(none) shared(n_send, block_info)
        for (uint64_t new_b_id = 0; new_b_id < n_blocks; ++new_b_id) {
            BSortInfo &info{block_info[new_b_id]};
            if (info.new_proc_id != info.old_proc_id) {
                if (proc_id == info.new_proc_id) { // receive
                    auto &b = _blocks[new_b_id];
                    if (log) log->debug("MIGRATE {} - recv: {} -> {}: corner = ({}, {}), size = ({}, {})",
                                        new_b_id, info.old_proc_id, proc_id,
                                        b.corner().first, b.corner().second,
                                        b.size().first, b.size().second);
                    _blocks[new_b_id].migrate(false, info.old_proc_id, _comm);
                }
                if (proc_id == info.old_proc_id) { // send
                    ++n_send;
                    auto &b = _blocks[new_b_id];
                    if (log) log->debug("MIGRATE {} - send: {} -> {}: corner = ({}, {}), size = ({}, {})",
                                        new_b_id, proc_id, info.new_proc_id,
                                        b.corner().first, b.corner().second,
                                        b.size().first, b.size().second);
                    _blocks[new_b_id].migrate(true, info.new_proc_id, _comm);
                }
            }
        }
        if (log) log->debug("migrated, sent {} blocks", n_send);
        if (proc_id == 0) {
            MPI_Reduce(MPI_IN_PLACE, &n_send, 1, MPI_UINT64_T, MPI_SUM, 0, _comm);
        } else {
            MPI_Reduce(&n_send, &n_send, 1, MPI_UINT64_T, MPI_SUM, 0, _comm);
        }
        if (log && proc_id == 0) log->debug("total sent {} blocks", n_send);

        // Update borders of locally stored blocks
        for (local_start_b_id = 0; 
             local_start_b_id < n_blocks && block_info[local_start_b_id].new_proc_id != proc_id; 
             ++local_start_b_id) {}
        for (local_end_b_id = local_start_b_id; 
             local_end_b_id < n_blocks && block_info[local_end_b_id].new_proc_id == proc_id; 
             ++local_end_b_id) {}
        local_n_blocks = local_end_b_id - local_start_b_id;
        if (log) log->debug("local borders: [{}, {})", local_start_b_id, local_end_b_id);
#endif
        // Set new ids (indices) to blocks
        for (uint64_t b_id = 0; b_id < n_blocks; ++b_id) {
            _blocks[b_id].set_id(b_id);
        }
        // Synchronize blocks distribution info
#ifdef MSK_USE_MPI
        sync_blocks_distribution();
#endif
    }
// #endif

#ifdef MSK_USE_MPI
    /**
     * @brief Synchronize blocks distribution information.
     */
    void sync_blocks_distribution() {
        int zero{0};
        BLAS::copy(n_blocks, &zero, 0, block_procs.data(), 1);
        for (uint64_t b_id = local_start_b_id; b_id < local_end_b_id; ++b_id) {
            block_procs[b_id] = proc_id;
        }
        MPI_Allreduce(MPI_IN_PLACE, block_procs.data(), n_blocks, MPI_INT, MPI_SUM, _comm);
    }
#endif

};


extern template class Manager<float>;

extern template class Manager<double>;

extern template class Manager<std::complex<float>>;

extern template class Manager<std::complex<double>>;


}
} // namespace MSk



#endif //MOSAIC_SKELETON_BLOCK_MANAGER_HPP
