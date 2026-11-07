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



#include "msk/utility/fourier.hpp"

namespace MSk {

namespace FFT {


template<>
void Daemon::init_memory<float>(uint64_t n_blocks, uint64_t size, uint64_t r) {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _f_bank;
    const uint64_t spec_size{(size / 2 + 1)};
    const uint64_t bank_spec_size_total{(bank.block_rows / 2 + 1) * bank.block_columns};
    if (bank.blocks.size() >= n_blocks && bank_spec_size_total >= spec_size * r) {
        return;
    }
    while (bank.free_blocks < bank.blocks.size()) {
        _cv.wait(lock);
    }
    _init_memory<float>(bank, fftwf_malloc, n_blocks, size, spec_size, r);
}

template<>
void Daemon::init_memory<std::complex<float>>(uint64_t n_blocks, uint64_t size, uint64_t r) {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _cf_bank;
    const uint64_t spec_size{size};
    const uint64_t bank_spec_size_total{bank.block_rows * bank.block_columns};
    if (bank.blocks.size() >= n_blocks && bank_spec_size_total >= spec_size * r) {
        return;
    }
    while (bank.free_blocks < bank.blocks.size()) {
        _cv.wait(lock);
    }
    _init_memory<std::complex<float>>(bank, fftwf_malloc, n_blocks, size, spec_size, r);
}

template<>
void Daemon::init_memory<double>(uint64_t n_blocks, uint64_t size, uint64_t r) {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _d_bank;
    const uint64_t spec_size{(size / 2 + 1)};
    const uint64_t bank_spec_size_total{(bank.block_rows / 2 + 1) * bank.block_columns};
    if (bank.blocks.size() >= n_blocks && bank_spec_size_total >= spec_size * r) {
        return;
    }
    while (bank.free_blocks < bank.blocks.size()) {
        _cv.wait(lock);
    }
    _init_memory<double>(bank, fftw_malloc, n_blocks, size, spec_size, r);
}

template<>
void Daemon::init_memory<std::complex<double>>(uint64_t n_blocks, uint64_t size, uint64_t r) {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _cd_bank;
    const uint64_t spec_size{size};
    const uint64_t bank_spec_size_total{bank.block_rows * bank.block_columns};
    if (bank.blocks.size() >= n_blocks && bank_spec_size_total >= spec_size * r) {
        return;
    }
    while (bank.free_blocks < bank.blocks.size()) {
        _cv.wait(lock);
    }
    _init_memory<std::complex<double>>(bank, fftw_malloc, n_blocks, size, spec_size, r);
}


template<>
void Daemon::free_memory<float>() {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _f_bank;
    while (bank.free_blocks < bank.blocks.size()) {
        _cv.wait(lock);
    }
    _free_memory(bank);
}

template<>
void Daemon::free_memory<std::complex<float>>() {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _cf_bank;
    while (bank.free_blocks < bank.blocks.size()) {
        _cv.wait(lock);
    }
    _free_memory(bank);
}

template<>
void Daemon::free_memory<double>() {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _d_bank;
    while (bank.free_blocks < bank.blocks.size()) {
        _cv.wait(lock);
    }
    _free_memory(bank);
}

template<>
void Daemon::free_memory<std::complex<double>>() {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _cd_bank;
    while (bank.free_blocks < bank.blocks.size()) {
        _cv.wait(lock);
    }
    _free_memory(bank);
}


template<>
int Daemon::lock_memory(float *&data, std::complex<float> *&spec) const {
    std::unique_lock<std::mutex> lock{_m};
    while (true) {
        int id = _lock_memory(_f_bank, data, spec);
        if (id != -1) {
            return id;
        }
        _cv.wait(lock);
    }
}

template<>
int Daemon::lock_memory(std::complex<float> *&data, std::complex<float> *&spec) const {
    std::unique_lock<std::mutex> lock{_m};
    while (true) {
        int id = _lock_memory(_cf_bank, data, spec);
        if (id != -1) {
            return id;
        }
        _cv.wait(lock);
    }
}

template<>
int Daemon::lock_memory(double *&data, std::complex<double> *&spec) const {
    std::unique_lock<std::mutex> lock{_m};
    while (true) {
        int id = _lock_memory(_d_bank, data, spec);
        if (id != -1) {
            return id;
        }
        _cv.wait(lock);
    }
}

template<>
int Daemon::lock_memory(std::complex<double> *&data, std::complex<double> *&spec) const {
    std::unique_lock<std::mutex> lock{_m};
    while (true) {
        int id = _lock_memory(_cd_bank, data, spec);
        if (id != -1) {
            return id;
        }
        _cv.wait(lock);
    }
}


template<>
void Daemon::unlock_memory<float>(int id) const {
    std::unique_lock<std::mutex> lock{_m};
    _unlock_memory(_f_bank, id);
    lock.unlock();
    _cv.notify_all();
}

template<>
void Daemon::unlock_memory<std::complex<float>>(int id) const {
    std::unique_lock<std::mutex> lock{_m};
    _unlock_memory(_cf_bank, id);
    lock.unlock();
    _cv.notify_all();
}

template<>
void Daemon::unlock_memory<double>(int id) const {
    std::unique_lock<std::mutex> lock{_m};
    _unlock_memory(_d_bank, id);
    lock.unlock();
    _cv.notify_all();
}

template<>
void Daemon::unlock_memory<std::complex<double>>(int id) const {
    std::unique_lock<std::mutex> lock{_m};
    _unlock_memory(_cd_bank, id);
    lock.unlock();
    _cv.notify_all();
}


template<>
fftwf_plan Daemon::get_plan<float>(int size, int num, bool fwd) const {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _f_bank;
    auto &plans = bank.plans;
    auto plan_it = plans.find({size, num, fwd});
    if (plan_it != plans.end()) {
        return plan_it->second;
    }

    int spec_size{size / 2 + 1};
    fftwf_plan plan;
    float *signal{nullptr};
    std::complex<float> *spectrum{nullptr};
    int id{_lock_memory(bank, signal, spectrum)};
    while (id == -1) {
        _cv.wait(lock);
        id = _lock_memory(bank, signal, spectrum);
    }
    if (fwd) {
        plan = fftwf_plan_many_dft_r2c(1, &size, num,
                signal, NULL, 1, size,
                reinterpret_cast<fftwf_complex *>(spectrum), NULL, 1, spec_size,
                FFTW_PLAN_POLICY);
    } else {
        plan = fftwf_plan_many_dft_c2r(1, &size, num,
                reinterpret_cast<fftwf_complex *>(spectrum), NULL, 1, spec_size,
                signal, NULL, 1, size,
                FFTW_PLAN_POLICY);
    }
    _unlock_memory(bank, id);
    plans[{size, num, fwd}] = plan;
    lock.unlock();
    _cv.notify_all();
    return plan;
}


template<>
fftwf_plan Daemon::get_plan<std::complex<float>>(int size, int num, bool fwd) const {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _cf_bank;
    auto &plans = bank.plans;
    auto plan_it = plans.find({size, num, fwd});
    if (plan_it != plans.end()) {
        return plan_it->second;
    }
    fftwf_plan plan;
    std::complex<float> *signal{nullptr};
    std::complex<float> *spectrum{nullptr};
    int id{_lock_memory(bank, signal, spectrum)};
    while (id == -1) {
        _cv.wait(lock);
        id = _lock_memory(bank, signal, spectrum);
    }
    if (fwd) {
        plan = fftwf_plan_many_dft(1, &size, num,
                reinterpret_cast<fftwf_complex *>(signal), NULL, 1, size,
                reinterpret_cast<fftwf_complex *>(spectrum), NULL, 1, size,
                FFTW_FORWARD, FFTW_PLAN_POLICY);
    } else {
        plan = fftwf_plan_many_dft(1, &size, num,
                reinterpret_cast<fftwf_complex *>(spectrum), NULL, 1, size,
                reinterpret_cast<fftwf_complex *>(signal), NULL, 1, size,
                FFTW_BACKWARD, FFTW_PLAN_POLICY);
    }
    _unlock_memory(bank, id);
    plans[{size, num, fwd}] = plan;
    lock.unlock();
    _cv.notify_all();
    return plan;
}

template<>
fftw_plan Daemon::get_plan<double>(int size, int num, bool fwd) const {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _d_bank;
    auto &plans = bank.plans;
    auto plan_it = plans.find({size, num, fwd});
    if (plan_it != plans.end()) {
        return plan_it->second;
    }

    int spec_size{size / 2 + 1};
    fftw_plan plan;
    double *signal{nullptr};
    std::complex<double> *spectrum{nullptr};
    int id{_lock_memory(bank, signal, spectrum)};
    while (id == -1) {
        _cv.wait(lock);
        id = _lock_memory(bank, signal, spectrum);
    }
    if (fwd) {
        plan = fftw_plan_many_dft_r2c(1, &size, num,
                signal, NULL, 1, size,
                reinterpret_cast<fftw_complex *>(spectrum), NULL, 1, spec_size,
                FFTW_PLAN_POLICY);
    } else {
        plan = fftw_plan_many_dft_c2r(1, &size, num,
                reinterpret_cast<fftw_complex *>(spectrum), NULL, 1, spec_size,
                signal, NULL, 1, size,
                FFTW_PLAN_POLICY);
    }
    _unlock_memory(bank, id);
    plans[{size, num, fwd}] = plan;
    lock.unlock();
    _cv.notify_all();
    return plan;
}

template<>
fftw_plan Daemon::get_plan<std::complex<double>>(int size, int num, bool fwd) const {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _cd_bank;
    auto &plans = bank.plans;
    auto plan_it = plans.find({size, num, fwd});
    if (plan_it != plans.end()) {
        return plan_it->second;
    }
    fftw_plan plan;
    std::complex<double> *signal{nullptr};
    std::complex<double> *spectrum{nullptr};
    int id{_lock_memory(bank, signal, spectrum)};
    while (id == -1) {
        _cv.wait(lock);
        id = _lock_memory(bank, signal, spectrum);
    }
    if (fwd) {
        plan = fftw_plan_many_dft(1, &size, num,
                reinterpret_cast<fftw_complex *>(signal), NULL, 1, size,
                reinterpret_cast<fftw_complex *>(spectrum), NULL, 1, size,
                FFTW_FORWARD, FFTW_PLAN_POLICY);
    } else {
        plan = fftw_plan_many_dft(1, &size, num,
                reinterpret_cast<fftw_complex *>(spectrum), NULL, 1, size,
                reinterpret_cast<fftw_complex *>(signal), NULL, 1, size,
                FFTW_BACKWARD, FFTW_PLAN_POLICY);
    }
    _unlock_memory(bank, id);
    plans[{size, num, fwd}] = plan;
    lock.unlock();
    _cv.notify_all();
    return plan;
}


template<>
void Daemon::transform<float>(int size, int num, float *signal, std::complex<float> *spectrum, bool fwd) const {
    fftwf_plan plan = get_plan<float>(size, num, fwd);
    if (fwd) {
        fftwf_execute_dft_r2c(plan, signal, reinterpret_cast<fftwf_complex *>(spectrum));
    } else {
        fftwf_execute_dft_c2r(plan, reinterpret_cast<fftwf_complex *>(spectrum), signal);
    }
}

template<>
void Daemon::transform<std::complex<float>>(int size, int num, std::complex<float> *signal, std::complex<float> *spectrum, bool fwd) const {
    fftwf_plan plan = get_plan<std::complex<float>>(size, num, fwd);
    if (fwd) {
        fftwf_execute_dft(plan, reinterpret_cast<fftwf_complex *>(signal), reinterpret_cast<fftwf_complex *>(spectrum));
    } else {
        fftwf_execute_dft(plan, reinterpret_cast<fftwf_complex *>(spectrum), reinterpret_cast<fftwf_complex *>(signal));
    }
}

template<>
void Daemon::transform<double>(int size, int num, double *signal, std::complex<double> *spectrum, bool fwd) const {
    fftw_plan plan = get_plan<double>(size, num, fwd);
    if (fwd) {
        fftw_execute_dft_r2c(plan, signal, reinterpret_cast<fftw_complex *>(spectrum));
    } else {
        fftw_execute_dft_c2r(plan, reinterpret_cast<fftw_complex *>(spectrum), signal);
    }
}

template<>
void Daemon::transform<std::complex<double>>(int size, int num, std::complex<double> *signal, std::complex<double> *spectrum, bool fwd) const {
    fftw_plan plan = get_plan<std::complex<double>>(size, num, fwd);
    if (fwd) {
        fftw_execute_dft(plan, reinterpret_cast<fftw_complex *>(signal), reinterpret_cast<fftw_complex *>(spectrum));
    } else {
        fftw_execute_dft(plan, reinterpret_cast<fftw_complex *>(spectrum), reinterpret_cast<fftw_complex *>(signal));
    }
}


template<>
void Daemon::destroy_plans<float>() {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _f_bank;
    _destroy_plans(bank, fftwf_destroy_plan);
}

template<>
void Daemon::destroy_plans<std::complex<float>>() {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _cf_bank;
    _destroy_plans(bank, fftwf_destroy_plan);
}

template<>
void Daemon::destroy_plans<double>() {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _d_bank;
    _destroy_plans(bank, fftw_destroy_plan);
}

template<>
void Daemon::destroy_plans<std::complex<double>>() {
    std::unique_lock<std::mutex> lock{_m};
    auto &bank = _cd_bank;
    _destroy_plans(bank, fftw_destroy_plan);
}


Daemon::~Daemon() {
    free_memory<float>();
    free_memory<double>();
    free_memory<std::complex<float>>();
    free_memory<std::complex<double>>();

    destroy_plans<float>();
    destroy_plans<double>();
    destroy_plans<std::complex<float>>();
    destroy_plans<std::complex<double>>();
}


}

}