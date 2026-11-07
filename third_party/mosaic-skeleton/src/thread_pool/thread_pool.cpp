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


#include "msk/thread_pool/thread_pool.hpp"


namespace MSk {

namespace exec {
ThreadPool::ThreadPool(const std::string &name, std::size_t low_watermark, std::size_t high_watermark,
                       std::size_t batch_size, std::size_t max_queue_size, std::size_t idle_time) :
                       _name{name}, _low_watermark{low_watermark}, _high_watermark{high_watermark},
                       _max_queue_size{max_queue_size}, _batch_size{batch_size}, _idle_time{idle_time},
                       all_threads{0}, free_threads{0}, state{State::stopped} {}

ThreadPool::~ThreadPool() {
    stop(true);
}

void ThreadPool::start() {
    std::unique_lock<std::mutex> _lock(mutex);
    if (state == State::running) {
        return;
    }
    // Finalize previous work in case ThreadPool was running before
    while (state != State::stopped) {
        stop_cond.wait(_lock);
    }
    // Launch workers
    for (std::size_t i = 0; i < _low_watermark; ++i) {
        std::thread new_thread(perform, this);
        new_thread.detach();
    }
    all_threads = _low_watermark;
    free_threads = 0;
    state = State::running;
}

void ThreadPool::stop(bool await) {
    std::unique_lock<std::mutex> _lock(mutex);
    if (state == State::stopped) {
        return;
    }
    if (state == State::running) {
        if (all_threads == 0) {
            state = State::stopped;
        } else {
            state = State::stopping;
            new_tasks.notify_all(); // To skip idle_time
        }
    }
    if (await) {
        while (state != State::stopped) {
            stop_cond.wait(_lock);
        }
    }
}

std::size_t ThreadPool::batch_size() const {
    std::lock_guard<std::mutex> lock{mutex};
    return _batch_size;
}

std::size_t ThreadPool::batch_size(std::size_t new_size) {
    std::lock_guard<std::mutex> lock{mutex};
    std::size_t old_size = _batch_size;
    if (new_size > 0) {
        _batch_size = new_size;
    }
    return old_size;
}

void perform(ThreadPool *pool) {
    using State = ThreadPool::State;
    std::unique_lock<std::mutex> _lock(pool->mutex);
    omp_set_num_threads(1);
    ++pool->free_threads;
    bool exit_flag{false};
    while (!pool->tasks.empty() || pool->state == State::running) {
        auto to_wait = std::chrono::system_clock::now() + pool->_idle_time;
        while (pool->tasks.empty() && pool->state == State::running) {
            if (pool->new_tasks.wait_until(_lock, to_wait) == std::cv_status::timeout &&
                    pool->all_threads > pool->_low_watermark) { // Current thread is not needed
                exit_flag = true;
                break;
            }
        }
        if (exit_flag || (pool->tasks.empty() && pool->state == State::stopping)) {
            break;
        }
        
        // Take a bunch of tasks to work on
        const int n_take = std::min(pool->tasks.size(), pool->_batch_size);
        std::vector<std::function<void()>> tasks(n_take);
        for (int i = 0; i < n_take; ++i) {
            tasks[i] = std::move(pool->tasks.begin()[i]);
        }
        pool->tasks.erase(pool->tasks.begin(), pool->tasks.begin() + n_take);
        // The work
        pool->free_threads -= 1;
        _lock.unlock();

        for (auto &task : tasks) {
            try {
                task();
            } catch (std::exception &ex) {
                std::cerr << "Exception caught while preforming task in ThreadPool: " << ex.what() << std::endl;
            }
        }

        _lock.lock();
        pool->free_threads += 1;
    } // while

    // Thread is dying
    pool->free_threads -= 1;
    if (--pool->all_threads == 0 && pool->state == State::stopping) { // This is last thread
        pool->state = State::stopped;
        pool->stop_cond.notify_all();
    }
}

std::size_t ThreadPool::max_threads() const {
    return _high_watermark;
}

std::size_t ThreadPool::min_threads() const {
    return _low_watermark;
}

}
} // namespace MSk

