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
// Created by b.valiakhmetov on 4/27/22.
//


/**
 * @file
 *
 * @brief Contains class eval::Queue for asinchronous matrix elements evaluation.
 */


#ifndef MOSAIC_SKELETON_EVALUATION_QUEUE_HPP
#define MOSAIC_SKELETON_EVALUATION_QUEUE_HPP


#include "msk_config.h"

#include <cstdint>
#include <deque>
#include <utility>
#include <vector>
#include <unordered_map>
#include <condition_variable>
#include <mutex>
#include <functional>

#include "request.hpp"
#include "smart_array/smart_array.hpp"

#include "msk/thread_pool/thread_pool.hpp"
#include "msk/evaluation/request.hpp"
#include "msk/utility/logging.hpp"
#include "msk/block/manager.hpp"

#include "msk/matrix/basic.hpp"


namespace MSk {

namespace block {
template<class DataType>
class Manager;
}


namespace eval {


/**
 * @brief Queue of matrix element Request.
 *
 * @details Manages matrix elements asynchronous evaluation.
 *
 * @tparam DataType Matrix element type.
 */
template<typename DataType>
class Queue {
private:
    /**
     * @brief Queue state.
     */
    enum class State {
        running, ///< Accept and evaluate elements; can't finish
        stopping, ///< Don't accept elements, evaluate already enqueued; on the way to finish
        stopped ///< Don't accept or evaluate elements; finished
    };
    State state; ///< Queue's current state

    /**
     * \brief Request that is currently being evaluated.
     */
    struct OnlineRequest {
        Request<DataType> request{}; ///< The request
        bool complete{false}; ///< Completion status
    #ifdef MSK_USE_MPI
        int target_proc_id{0}; ///< Process that hosted the request
        int eval_proc_id{0}; ///< Process that performs evaluation
    #endif

    public:
        OnlineRequest() {}
    #ifdef MSK_USE_MPI
        OnlineRequest(Request<DataType> &&r, int tp = 0, int ep = 0, bool comp = false) :
                request{std::move(r)}, complete{comp}, target_proc_id{tp}, eval_proc_id{ep} {}
    #else
        OnlineRequest(Request<DataType> &&r, bool comp = false) :
                request{std::move(r)}, complete{comp} {}
    #endif
        ~OnlineRequest() = default;
        OnlineRequest(const OnlineRequest &) = delete;
        OnlineRequest(OnlineRequest &&other) : request{std::move(other.request)},
                complete{other.complete}
#ifdef MSK_USE_MPI
                , target_proc_id{other.target_proc_id}, eval_proc_id{other.eval_proc_id}
#endif
                {
                    other.complete = false;
#ifdef MSK_USE_MPI
                    other.target_proc_id = -1;
                    other.eval_proc_id = -1;
#endif
                }
        OnlineRequest &operator=(const OnlineRequest &) = delete;
        OnlineRequest &operator=(OnlineRequest &&other) {
            request = std::move(other.request);
            complete = other.complete;
            other.complete = false;
#ifdef MSK_USEMPI
            target_proc_id = other.target_proc_id;
            other.target_proc_id = -1;
            eval_proc_id = other.eval_proc_id;
            other.eval_proc_id = -1;
#endif
            return *this;
        }
    };

    std::deque<Request<DataType>> _stock_reqs; ///< Enqueued requests to be evaluated


    std::unordered_map<int, OnlineRequest> _online_reqs; ///< Requests currently being evaluated
    int _next_req_id{0}; ///< Id to assign to request that next will be put to evaluate
    std::deque<int> _complete_online_ids; ///< // Queue of evaluated but not completed OnlineRequest ids

    int proc_id{0}, n_proc{1}; ///< Process id and total number of processes

    // Containers synchronization
    std::mutex q_m, compl_m, state_m; ///< Mutexes to lock containers
    std::condition_variable q_cv, compl_cv, state_cv; ///< Condvars to wait on while containers are empty
    std::thread _extract_master, _complete_master; ///< Working threads

    uint64_t n_ready{0}; ///< Current number of enqueued elements, ready to evaluate
    uint64_t batch_s{0}; ///< Batch size to take from queue and set for evaluation
    block::Manager<DataType> *_manager; ///< Blocks manager to transfer evaluated elements

    std::shared_ptr<log::Logger> log{nullptr}; ///< Logger

public:
    /**
     * @brief Set Queue's parameters.
     *
     * @param[in] manager Manager halding blocks.
     * @param[in] batch_size Portion of matrix elements to form as task for evaluation.
     * @param[in] log_ptr Logger
     */
    Queue(block::Manager<DataType> *manager, uint64_t batch_size = 0,
          std::shared_ptr<log::Logger> log_ptr = nullptr) :
          state{State::stopped}, batch_s{batch_size},
          _manager{manager},
          log{std::move(log_ptr)} {  }

    ~Queue() {
        stop();
    }

    /**
     * @brief Start accepting requests and the evaluation process on a matrix.
     *
     * @tparam Matrix Abstract matrix type, required to have at least element(i, j) method.
     * @param[in] matrix Matrix.
     * @param[in] pool Thread pool to put evaluation tasks.
     */
    template<class Matrix>
    void start(Matrix &matrix, exec::ThreadPool *pool) {
        if (log) log->info("Queue start");
        State st;
        {
            std::lock_guard<std::mutex> lock(state_m);
            st = state;
        }
        if (st != State::stopped) {
            stop();
        }
        {
            std::lock_guard<std::mutex> lock(state_m);
            state = State::running;
            _extract_master = std::thread(&Queue::template run_extract<Matrix>, this, std::ref(matrix), pool);
            _complete_master = std::thread(&Queue::run_complete, this);
        }
    }

    /**
     * @brief Stop accepting requests, enqueued elements will be evaluated.
     * 
     * @param[in] await If true then blocks until evaluation process finishes.
     */
    void stop(bool await = true) {
        if (log) log->info("Queue stop");
        std::unique_lock<std::mutex> lock(state_m);
        if (state == State::stopped) {
            return;
        }
        State old = state;
        state = State::stopping;
        lock.unlock();
        if (old == State::running) {
            q_cv.notify_all();
            compl_cv.notify_all();
        }
        if (await) {
            _extract_master.join();
            _complete_master.join();
            lock.lock();
            state = State::stopped;
        }
    }

    /**
     * @brief Enqueue elements to be evaluated.
     *
     * @param[in, out] request Request of elements, freed on output.
     * @return Enqueuement success.
     */
    bool add(Request<DataType> &&request) {
        std::unique_lock<std::mutex> lock(q_m);
        n_ready += request.size();
        if (batch_s == 0 || request.size() <= batch_s) {
            _stock_reqs.push_back(std::move(request));
        } else { // Split too large request
            while (!request.empty()) {
                Request<DataType> short_req{request.split(batch_s)};
                _stock_reqs.push_back(std::move(short_req));
            }
        }
        lock.unlock();
        // Notify evaluation managing thread
        q_cv.notify_one();
        return true;
    }

    /**
     * @brief Enqueue array of requests to be evaluated.
     *
     * @param[in] n Number of requests.
     * @param[in, out] requests Array of requests, all are freed on output.
     * @return Enqueuement success.
     */
    bool add(const uint64_t n, Request<DataType> *requests) {
        std::unique_lock<std::mutex> lock(q_m);
        for (uint64_t i = 0; i < n; ++i) {
            auto &request = requests[i];
            n_ready += request.size();
            if (batch_s == 0 || request.size() <= batch_s) {
                _stock_reqs.push_back(std::move(request));
            } else { // Split too large request
                while (!request.empty()) {
                    Request<DataType> short_req{request.split(batch_s)};
                    _stock_reqs.push_back(std::move(short_req));
                }
            }
        }
        lock.unlock();
        // Notify evaluation managing thread
        q_cv.notify_one();
        return true;
    }

private:
    /**
     * @brief Finalize evaluation for array of requests.
     *
     * @param[in] n Nmber of requests.
     * @param[in] requests Array of requests.
     */
    void complete_requests(uint64_t n, const Request<DataType> *requests) {
        _manager->transfer_elements(n, requests);
    }

    /**
     * @brief Function for the queue managing thread to run.
     *
     * @details Extrancts batches of elements from the queue and sets evaluation tasks.
     *
     * @tparam Matrix Matrix type.
     * @param[in] matrix Matrix.
     * @param[in] pool Thread pool to put evaluation tasks.
     */
    template<class Matrix>
    void run_extract(Matrix &matrix, exec::ThreadPool *pool) {
        if (log) log->info("Queue start extraction");
        uint64_t total_iter{0}, total_eval{0}, total_req{0};
        while (true) {
            constexpr uint64_t N_BATCH = 1; // Currently extract requests by one batch
            std::vector<OnlineRequest> eval_online_reqs;
            std::unique_lock<std::mutex> lock(q_m);
            uint64_t n_reqs{0}, n_eval{0};
            if (n_ready > 0) { // Extract requests from the queue's head
                // Count number of requests to meet batch size requrement
                if (batch_s == 0 || n_ready <= batch_s * N_BATCH) {
                    n_eval = n_ready;
                    n_reqs = _stock_reqs.size();
                } else {
                    for (auto req_it = _stock_reqs.begin(); req_it < _stock_reqs.end(); ++req_it) {
                        if (n_eval + req_it->size() > batch_s * N_BATCH) {
                            break;
                        }
                        n_eval += req_it->size();
                        ++n_reqs;
                    }
                }
                n_ready -= n_eval;
                eval_online_reqs.resize(n_reqs);
                auto req_it = _stock_reqs.begin();
                for (uint64_t i = 0; i < n_reqs; ++i) {
                    eval_online_reqs[i] = {std::move(*req_it)};
                    ++req_it;
                }
                _stock_reqs.erase(_stock_reqs.begin(), _stock_reqs.begin() + n_reqs);
                std::vector<int> eval_online_reqs_ids(n_reqs);
                // Set extracted requests as "online"
                lock.unlock();
                {
                    std::lock_guard<std::mutex> compl_lock{compl_m};
                    for (uint64_t i = 0; i < n_reqs; ++i) {
                        eval_online_reqs_ids[i] = _next_req_id++;
                        _online_reqs.emplace(eval_online_reqs_ids[i], std::move(eval_online_reqs[i]));
                    }
                }
                // Set evaluation task
                pool->execute(&Queue::evaluate_batch<Matrix>, this, std::ref(matrix), std::move(eval_online_reqs_ids));
                total_eval += n_eval;
                total_req += n_reqs;
                ++total_iter;
                if (log) log->debug("Eval: {} elems in {} reqs", n_eval, n_reqs);
            } else { // No elements to evaluate
                {
                    // Check if need to stop evaluation
                    std::lock_guard<std::mutex> s_lock(state_m);
                    if (state == State::stopping || state == State::stopped) {
                        break;
                    }
                }
                // Wait for elements
                q_cv.wait(lock);
            }
        }
        // Signalize about evaluation finish
        state_cv.notify_all();
        if (log) log->info("Queue end run");
        if (log) log->info("Total: {} elems, {} requests, {} batches", total_eval, total_req, total_iter);
    }

    /**
     * @brief Function for the requests completing thread to run.
     */
    void run_complete() {
        if (log) log->debug("Queue start completion");
        std::vector<Request<DataType>> complete_local;
        while (true) {
            std::unique_lock<std::mutex> lock{compl_m};
            // Check requests that were evaluated but not marked as complete
            uint64_t n_complete_local{_complete_online_ids.size()};
            if (n_complete_local > 0) {
                complete_local.resize(n_complete_local);
                // Extract the requests from the list of "online"
                for (uint64_t i = 0; i < n_complete_local; ++i) {
                    auto it = _online_reqs.find(_complete_online_ids[i]);
                    complete_local[i] = std::move(it->second.request);
                    _online_reqs.erase(it);
                }
                _complete_online_ids.clear();
                lock.unlock();
                // Finalize requests evaluation
                complete_requests(n_complete_local, complete_local.data());
            } else { // No requests to complete
                {
                    // Check if need to stop completion
                    std::lock_guard<std::mutex> s_lock(state_m);
                    if (state == State::stopping || state == State::stopped) {
                        break;
                    }
                }
                // Wait for requests to complete
                compl_cv.wait(lock);
            }
        }
    }

    /**
     * @brief Evaluate a batch of requests.
     *
     * @details This routine is passed as task to exec::ThreadPool.
     *
     * @tparam Matrix Matrix type.
     * @param[in] matrix Matrix that provides elements evaluation.
     * @param[in] eval_req_ids Ids of "online" requests to be evaluated at current call.
     */
    template<class Matrix>
    void evaluate_batch(Matrix &matrix, std::vector<int> eval_req_ids) {
        // Extract requests to evaluate
        uint64_t n_reqs{eval_req_ids.size()};
        SmartArray<const Request<DataType> *> eval_reqs{n_reqs};
        std::unique_lock<std::mutex> lock{compl_m};
        for (uint64_t i = 0; i < n_reqs; ++i) {
            eval_reqs[i] = &_online_reqs.at(eval_req_ids[i]).request;
        }
        lock.unlock();
        // Perform the actual elements evaluation
        for (uint64_t i = 0; i < n_reqs; ++i) {
            const auto &req = *eval_reqs[i];
            DataType *result = req.result_ptr();
            std::vector<uint64_t> rows, cols;
            req.extract_compact(rows, cols);
            switch (req.type()) {
            case Request<DataType>::Type::element:
                result[0] = matrix.element(rows[0], cols[0]);
                break;
            case Request<DataType>::Type::row:
                matrix.row(rows[0], cols, result);
                break;
            case Request<DataType>::Type::column:
                matrix.column(rows, cols[0], result);
                break;
            case Request<DataType>::Type::block:
                matrix.block(rows, cols, result);
                break;
            case Request<DataType>::Type::list:
                matrix.list(rows, cols, result);
            default:
                rows.clear();
                cols.clear();
                req.extract_list(rows, cols);
                for (uint64_t j = 0; j < req.size(); ++j) {
                    result[j] = matrix.element(rows[j], cols[j]);
                }
            };
        }
        // Notify completion thread that current requests are ready to finalize
        lock.lock();
        _complete_online_ids.insert(_complete_online_ids.cend(), eval_req_ids.cbegin(), eval_req_ids.cend());
        lock.unlock();
        compl_cv.notify_one();
    }
};


extern template class Queue<float>;

extern template class Queue<double>;

extern template class Queue<std::complex<float>>;

extern template class Queue<std::complex<double>>;


extern template void Queue<float>::start(matrix::Elementary<float> &matrix, exec::ThreadPool *pool);
extern template void Queue<float>::start(matrix::Dense<float, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<float>::start(matrix::Dense<float, 3> &matrix, exec::ThreadPool *pool);
extern template void Queue<float>::start(matrix::BiPoint<float, float, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<float>::start(matrix::BiPoint<float, float, 3> &matrix, exec::ThreadPool *pool);

extern template void Queue<double>::start(matrix::Elementary<double> &matrix, exec::ThreadPool *pool);
extern template void Queue<double>::start(matrix::Dense<double, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<double>::start(matrix::Dense<double, 3> &matrix, exec::ThreadPool *pool);
extern template void Queue<double>::start(matrix::BiPoint<double, double, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<double>::start(matrix::BiPoint<double, double, 3> &matrix, exec::ThreadPool *pool);

extern template void Queue<std::complex<float>>::start(matrix::Elementary<std::complex<float>> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<float>>::start(matrix::Dense<std::complex<float>, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<float>>::start(matrix::Dense<std::complex<float>, 3> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<float>>::start(matrix::BiPoint<std::complex<float>, float, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<float>>::start(matrix::BiPoint<std::complex<float>, float, 3> &matrix, exec::ThreadPool *pool);

extern template void Queue<std::complex<double>>::start(matrix::Elementary<std::complex<double>> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<double>>::start(matrix::Dense<std::complex<double>, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<double>>::start(matrix::Dense<std::complex<double>, 3> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<double>>::start(matrix::BiPoint<std::complex<double>, double, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<double>>::start(matrix::BiPoint<std::complex<double>, double, 3> &matrix, exec::ThreadPool *pool);


extern template void Queue<float>::run_extract(matrix::Elementary<float> &matrix, exec::ThreadPool *pool);
extern template void Queue<float>::run_extract(matrix::Dense<float, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<float>::run_extract(matrix::Dense<float, 3> &matrix, exec::ThreadPool *pool);
extern template void Queue<float>::run_extract(matrix::BiPoint<float, float, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<float>::run_extract(matrix::BiPoint<float, float, 3> &matrix, exec::ThreadPool *pool);

extern template void Queue<double>::run_extract(matrix::Elementary<double> &matrix, exec::ThreadPool *pool);
extern template void Queue<double>::run_extract(matrix::Dense<double, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<double>::run_extract(matrix::Dense<double, 3> &matrix, exec::ThreadPool *pool);
extern template void Queue<double>::run_extract(matrix::BiPoint<double, double, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<double>::run_extract(matrix::BiPoint<double, double, 3> &matrix, exec::ThreadPool *pool);

extern template void Queue<std::complex<float>>::run_extract(matrix::Elementary<std::complex<float>> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<float>>::run_extract(matrix::Dense<std::complex<float>, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<float>>::run_extract(matrix::Dense<std::complex<float>, 3> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<float>>::run_extract(matrix::BiPoint<std::complex<float>, float, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<float>>::run_extract(matrix::BiPoint<std::complex<float>, float, 3> &matrix, exec::ThreadPool *pool);

extern template void Queue<std::complex<double>>::run_extract(matrix::Elementary<std::complex<double>> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<double>>::run_extract(matrix::Dense<std::complex<double>, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<double>>::run_extract(matrix::Dense<std::complex<double>, 3> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<double>>::run_extract(matrix::BiPoint<std::complex<double>, double, 2> &matrix, exec::ThreadPool *pool);
extern template void Queue<std::complex<double>>::run_extract(matrix::BiPoint<std::complex<double>, double, 3> &matrix, exec::ThreadPool *pool);


extern template void Queue<float>::evaluate_batch(matrix::Elementary<float> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<float>::evaluate_batch(matrix::Dense<float, 2> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<float>::evaluate_batch(matrix::Dense<float, 3> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<float>::evaluate_batch(matrix::BiPoint<float, float, 2> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<float>::evaluate_batch(matrix::BiPoint<float, float, 3> &matrix, std::vector<int> eval_req_ids);

extern template void Queue<double>::evaluate_batch(matrix::Elementary<double> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<double>::evaluate_batch(matrix::Dense<double, 2> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<double>::evaluate_batch(matrix::Dense<double, 3> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<double>::evaluate_batch(matrix::BiPoint<double, double, 2> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<double>::evaluate_batch(matrix::BiPoint<double, double, 3> &matrix, std::vector<int> eval_req_ids);

extern template void Queue<std::complex<float>>::evaluate_batch(matrix::Elementary<std::complex<float>> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<std::complex<float>>::evaluate_batch(matrix::Dense<std::complex<float>, 2> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<std::complex<float>>::evaluate_batch(matrix::Dense<std::complex<float>, 3> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<std::complex<float>>::evaluate_batch(matrix::BiPoint<std::complex<float>, float, 2> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<std::complex<float>>::evaluate_batch(matrix::BiPoint<std::complex<float>, float, 3> &matrix, std::vector<int> eval_req_ids);

extern template void Queue<std::complex<double>>::evaluate_batch(matrix::Elementary<std::complex<double>> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<std::complex<double>>::evaluate_batch(matrix::Dense<std::complex<double>, 2> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<std::complex<double>>::evaluate_batch(matrix::Dense<std::complex<double>, 3> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<std::complex<double>>::evaluate_batch(matrix::BiPoint<std::complex<double>, double, 2> &matrix, std::vector<int> eval_req_ids);
extern template void Queue<std::complex<double>>::evaluate_batch(matrix::BiPoint<std::complex<double>, double, 3> &matrix, std::vector<int> eval_req_ids);


}
}



#endif //MOSAIC_SKELETON_EVALUATION_QUEUE_HPP
