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
 * @brief Defines class ThreadPool.
 */


#ifndef MOSAIC_SKELETON_THREAD_POOL_THREAD_POOL_HPP
#define MOSAIC_SKELETON_THREAD_POOL_THREAD_POOL_HPP


#include <omp.h>

#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <chrono>
#include <exception>
#include <iostream>


namespace MSk {

namespace exec {


class ThreadPool;

/**
 * @brief Function that all ThreadPool workers are running.
 *
 * @details Extracts tasks from the queue and executes them.
 * If there are no tasks and too many threads, then finishes.
 *
 * @param[in] pool ThreadPool for thread to work on.
 */
void perform(ThreadPool *pool);

/**
 * @brief Pool of C++ threads for arbitrary functions parallel execution.
 */
class ThreadPool {
    /**
     * @brief Pool state
     */
    enum class State {
        running, ///< Receiving and executing tasks
        stopping, ///< Not accepting task and finalizing existing ones
        stopped ///< Not accepting and not executing tasks
    };

public:
    /**
     * @brief Construct ThreadPool
     *
     * @param[in] name Name of the pool.
     * @param[in] low_watermark Minimal number of working threads.
     * @param[in] high_watermark Maximal number of working threads.
     * @param[in] batch_size Number of tasks to extract by a single worker.
     * @param[in] max_queue_size Maximal number of tasks in the queue.
     * @param[in] idle_time Time (ms) for extra working threads to wait for a new task, and die if no such.
     */
    ThreadPool(const std::string &name, std::size_t low_watermark = 4, std::size_t high_watermark = 8,
               std::size_t batch_size = 1, std::size_t max_queue_size = 0, std::size_t idle_time = 250);

    ~ThreadPool();

    /**
     * @brief Signalize the pool to stop.
     *
     * @details No tasks are accepted for execution. All existing tasks will be finished.
     *
     * @param[in] await If true then the call will block untill all enqueued tasks will be finished.
     */
    void stop(bool await = false);

    /**
     * @brief Start accepting and executing tasks.
     */
    void start();

    std::size_t batch_size() const;

    std::size_t batch_size(std::size_t new_size);


    friend void perform(ThreadPool *pool);


    /**
     * @brief Schedule prioritized task to be executed by the pool.
     *
     * @details ThreadPool does not take care of function's return value.
     * If needed, caller should manually set it by a reference or pointer passed in arguments.
     *
     * @tparam F Function type.
     * @tparam Types Arguments type.
     * @param[in] func Function to be set as a task.
     * @param[in] args Arguments for the function to use.
     * @return Task addition success.
     */
    template<typename F, typename... Types>
    bool execute(F &&func, Types... args) {
        // Prepare "task"
        auto exec = std::bind(std::forward<F>(func), std::forward<Types>(args)...);

        std::unique_lock<std::mutex> _lock(this->mutex);
        if ((_max_queue_size > 0 && tasks.size() >= _max_queue_size) || state != State::running) {
            return false;
        }

        // Enqueue new task
        tasks.push_back(exec);
        if (free_threads > 0) { // Wake up a worker of there is a sleeping one
            new_tasks.notify_one();
        } else if (all_threads < _high_watermark) { // Otherwise add new worker
            ++all_threads;
            std::thread new_thread(perform, this);
            new_thread.detach();
        }
        return true;
    }

    /**
     * @brief Maximal possible number of worker threads.
     *
     * @return Number of workers.
     */
    std::size_t max_threads() const;

    /**
     * @brief Minimal possible number of alive worker threads.
     *
     * @return Number of workers.
     */
    std::size_t min_threads() const;

    ThreadPool(const ThreadPool &) = delete;
    ThreadPool(ThreadPool &&) = delete;
    ThreadPool &operator=(const ThreadPool &) = delete;
    ThreadPool &operator=(ThreadPool &&) = delete;

private:
    mutable std::mutex mutex; ///< State and queue lock

    std::condition_variable new_tasks; ///< Condvar to wait new tasks on
    std::condition_variable stop_cond; ///< Condvar to wait while stopping the pool

    const std::string _name; ///< Pool name
    const std::size_t _low_watermark, _high_watermark, _max_queue_size; ///< Workers and tasks bounds
    std::size_t _batch_size{1}; ///< Number of tasks to extract by a single worker
    const std::chrono::milliseconds _idle_time; ///< Time to wait for a new task, and die if no task

    std::deque<std::function<void()>> tasks; ///< Queue

    std::size_t all_threads, free_threads; ///< Counters of threads

    State state; ///< State flag for workers and pool to communicate
};

}
} // namespace MSk

#endif // MOSAIC_SKELETON_THREAD_POOL_THREAD_POOL_HPP
