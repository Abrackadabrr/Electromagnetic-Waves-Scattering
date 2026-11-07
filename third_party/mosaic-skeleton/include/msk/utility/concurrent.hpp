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
// Created by b.valiakhmetov on 2/6/24.
//


#ifndef MOSAIC_SKELETON_UTILITY_CONCURRENT_HPP
#define MOSAIC_SKELETON_UTILITY_CONCURRENT_HPP


#include <complex>
#include <functional>
#include <vector>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <set>

#include "smart_array/smart_array.hpp"


namespace MSk {

namespace exec {


template<class Item>
class ConcurrentExecutor {
public:
    using id_t = uint64_t;
private:
    struct node_t {
        Item *item_ptr{nullptr};
        std::vector<node_t *> concur{};
        int n_lock{0};
        bool complete{false};
    };
    std::unordered_map<id_t, node_t> nodes{};
    uint64_t n_nodes{0}, n_complete{}, n_active{}, n_locked{};

    mutable std::mutex _m;
    mutable std::condition_variable _cv;

public:
    void set_item(id_t id, Item *item) {
        std::lock_guard<std::mutex> lock(_m);
        node_t node;
        node.item_ptr = item;
        nodes[id] = node;
        ++n_nodes;
    }

    void reset() {
        nodes.clear();
        n_nodes = n_complete = n_active = n_locked = 0;
    }

    void mark_concurrent(id_t id, std::vector<id_t> &&concurrent) {
        std::lock_guard<std::mutex> lock(_m);
        auto it = nodes.find(id);
        if (it == nodes.end()) {
            return;
        }
        auto &concur = it->second.concur;
        int n = concur.size();
        concur.resize(n + concurrent.size());
        for (uint64_t i = 0; i < concurrent.size(); ++i) {
            concur[n + i] = &nodes.at(concurrent[i]);
        }
    }

    void restart_execution() {
        std::lock_guard<std::mutex> lock(_m);
        n_complete = n_active = n_locked = 0;
        n_nodes = nodes.size();
        for (auto &node : nodes) {
            node.second.n_lock = 0;
            node.second.complete = false;
        }
    }

    bool finished() const {
        std::lock_guard<std::mutex> lock(_m);
        return n_complete == n_nodes;
    }

    std::pair<id_t, Item *> get_next() {
        std::unique_lock<std::mutex> lock(_m);
        while (n_complete + n_active + n_locked == n_nodes) { // Check if there is available work to do
            if (n_locked == 0) { // All work is done
                return {-1, nullptr};
            }
            _cv.wait(lock); // Wait if no free items to work on
        }
        // Look for an item to execute
        for (auto &id_node : nodes) {
            auto &node = id_node.second;
            if (!node.complete && node.n_lock == 0) {
                for (auto conc_ptr : node.concur) {
                    if (!conc_ptr->complete) {
                        if (conc_ptr->n_lock == 0) {
                            ++n_locked;
                        }
                        conc_ptr->n_lock += 1;
                    }
                }
                ++n_active;
                node.n_lock = -1;
                return {id_node.first, node.item_ptr};
            }
        }
        return {-1, nullptr};
    }

    bool set_complete(id_t id) {
        std::lock_guard<std::mutex> lock(_m);
        auto &node = nodes.at(id);
        node.complete = true;
        node.n_lock = 0;
        --n_active;
        ++n_complete;
        for (auto conc_ptr : node.concur) {
            if (!conc_ptr->complete) {
                conc_ptr->n_lock -= 1;
                if (conc_ptr->n_lock == 0) {
                    --n_locked;
                    _cv.notify_one();
                }
            }
        }
        return n_complete == n_nodes;
    }
};


template<typename DataType>
class MemoryBank {
public:
    class Cell {
        uint64_t _size{0};
        DataType *_data{nullptr};
    public:
        explicit Cell() = default;
        Cell(uint64_t size, DataType *data) : _size{size}, _data{data} {}
        Cell(const Cell &) = default;
        Cell(Cell &&other) noexcept : _size{other._size}, _data{other._data} {
            other._size = 0;
            other._data = nullptr;
        }
        Cell &operator=(const Cell &) = default;
        Cell &operator=(Cell &&other) noexcept {
            if (this == &other) {
                return *this;
            }
            if (_size > 0) {
                std::free(_data);
            }
            _size = other._size;
            _data = other._data;
            other._size = 0;
            other._data = nullptr;
            return *this;
        }
        DataType *data() const {
            return _data;
        }
        uint64_t size() const {
            return _size;
        }
    };

private:
    std::function<bool(const Cell &, const Cell &)> comp = [](const Cell &n1, const Cell &n2) {
            return std::make_pair(n1.size(), n1.data()) < std::make_pair(n2.size(), n2.data());
    };
    std::set<Cell, decltype(comp)> bank{comp};
    mutable std::mutex m;

public:
    ~MemoryBank() {
        for (auto & node : bank) {
            std::free(node.data());
        }
    }

    Cell get(uint64_t size) {
        std::lock_guard<std::mutex> lock{m};
        Cell required{size, nullptr};
        auto it = bank.lower_bound(required);
        if (it == bank.end()) {
            if (bank.empty()) {;
                return Cell{size, static_cast<DataType *>(std::malloc(size * sizeof(DataType)))};
            }
            auto largest = bank.end();
            --largest;
            DataType *old_data{largest->data()};
            bank.erase(largest);
            return Cell{size, static_cast<DataType *>(std::realloc(old_data, size * sizeof(DataType)))};
        }
        Cell node(*it);
        bank.erase(it);
        return node;
    }

    void ret(Cell &&node) {
        if (node.size() > 0) {
            std::lock_guard<std::mutex> lock{m};
            bank.insert(std::move(node));
        }
    }
};


extern template class MemoryBank<float>;

extern template class MemoryBank<double>;

extern template class MemoryBank<std::complex<float>>;

extern template class MemoryBank<std::complex<double>>;


}
}

#endif //MOSAIC_SKELETON_UTILITY_CONCURRENT_HPP
