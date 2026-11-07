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
// Created by b.valiakhmetov on 2/20/24.
//


/**
 * @file
 *
 * @brief Fourier transform routines.
 */


#ifndef MOSAIC_SKELETON_UTILITY_FOURIER_HPP
#define MOSAIC_SKELETON_UTILITY_FOURIER_HPP


#include <fftw3.h>

#include <mutex>
#include <condition_variable>
#include <map>
#include <complex>
#include <functional>
#include <unordered_map>


namespace MSk {



namespace FFT {

constexpr decltype(FFTW_ESTIMATE) FFTW_PLAN_POLICY = FFTW_MEASURE;

template<typename DataType>
class Problem;

template<>
class Problem<float> {
public:
    using data_type = float;
    using spec_type = std::complex<float>;
    using plan_type = fftwf_plan;
};

template<>
class Problem<std::complex<float>> {
public:
    using data_type = std::complex<float>;
    using spec_type = std::complex<float>;
    using plan_type = fftwf_plan;
};

template<>
class Problem<double> {
public:
    using data_type = double;
    using spec_type = std::complex<double>;
    using plan_type = fftw_plan;
};

template<>
class Problem<std::complex<double>> {
public:
    using data_type = std::complex<double>;
    using spec_type = std::complex<double>;
    using plan_type = fftw_plan;
};


class Daemon {
private:
    template<typename DataType>
    struct MemoryBlock {
        typename Problem<DataType>::data_type *data{nullptr};
        typename Problem<DataType>::spec_type *spec{nullptr};
        bool in_use{false};
    };

    template<typename DataType>
    struct MemoryBank {
        uint64_t free_blocks{0};
        uint64_t block_rows{0}, block_columns{0};
        std::unordered_map<int, MemoryBlock<DataType>> blocks;
        std::map<std::tuple<int, int, bool>, typename Problem<DataType>::plan_type> plans{};
    };

    mutable MemoryBank<float> _f_bank;
    mutable MemoryBank<std::complex<float>> _cf_bank;
    mutable MemoryBank<double> _d_bank;
    mutable MemoryBank<std::complex<double>> _cd_bank;

    mutable std::mutex _m;
    mutable std::condition_variable _cv;

private:
    Daemon() = default;

    ~Daemon();

    template<typename DataType>
    void _init_memory(MemoryBank<DataType> &bank, const std::function<void *(size_t)> &alloc,
                      uint64_t n_blocks, uint64_t data_size, uint64_t spec_size, uint64_t rank) {
        using SpecType = typename Problem<DataType>::spec_type;
        for (auto &block : bank.blocks) {
            fftw_free(block.second.data);
            fftw_free(block.second.spec);
        }
        bank.free_blocks = n_blocks;
        bank.block_rows = data_size;
        bank.block_columns = rank;
        bank.blocks.clear();
        for (int id = 0; id < static_cast<int>(n_blocks); ++id) {
            auto &block = bank.blocks[id];
            block.in_use = false;
            block.data = static_cast<DataType *>(alloc(sizeof(DataType) * data_size * rank * 2));
            block.spec = static_cast<SpecType *>(alloc(sizeof(SpecType) * spec_size * rank * 2));
        }
    }

    template<typename DataType>
    void _free_memory(MemoryBank<DataType> &bank) {
        for (auto &block : bank.blocks) {
            fftw_free(block.second.data);
            fftw_free(block.second.spec);
        }
        bank.free_blocks = 0;
        bank.block_rows = 0;
        bank.block_columns = 0;
        bank.blocks.clear();
    }

    template<typename DataType>
    int _lock_memory(MemoryBank<DataType> &bank, DataType *&data, typename Problem<DataType>::spec_type *&spec) const {
        data = nullptr;
        spec = nullptr;
        if (bank.free_blocks == 0) {
            return -1;
        }
        for (auto &id_block : bank.blocks) {
            int id = id_block.first;
            auto &block = id_block.second;
            if (!block.in_use) {
                data = block.data;
                spec = block.spec;
                block.in_use = true;
                --bank.free_blocks;
                return id;
            }
        }
        return -1;
    }

    template<typename DataType>
    void _unlock_memory(MemoryBank<DataType> &bank, int id) const {
        auto &block = bank.blocks.at(id);
        if (block.in_use) {
            ++bank.free_blocks;
        }
        block.in_use = false;
    }

    template<typename DataType>
    void _destroy_plans(MemoryBank<DataType> &bank,
                        const std::function<void(typename Problem<DataType>::plan_type)> &destroyer) {
        for (auto &plan : bank.plans) {
            destroyer(plan.second);
        }
    }

public:
    static Daemon& get() {
        static Daemon instance;
        return instance;
    }

    template<typename DataType>
    void init_memory(uint64_t n_blocks, uint64_t m, uint64_t r);

    template<typename DataType>
    void free_memory();

    template<typename DataType>
    int lock_memory(DataType *&data, std::complex<decltype(std::abs(DataType{}))> *&spec) const;

    template<typename DataType>
    void unlock_memory(int id) const;

    template<typename DataType>
    typename Problem<DataType>::plan_type get_plan(int size, int num, bool fwd) const;

    template<typename DataType>
    void transform(int size, int num, DataType *signal, typename Problem<DataType>::spec_type *spectrum, bool fwd) const;

    template<typename DataType>
    void destroy_plans();
};



}
}


#endif //MOSAIC_SKELETON_UTILITY_FOURIER_HPP
