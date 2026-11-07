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
// Created by bulat_v on 17.01.24.
//


/**
 * @file
 *
 * @brief Containg logging interface utilities to spdlog library.
 */

#ifndef MOSAIC_SKELETON_UTILITY_CONVOLUTION_HPP
#define MOSAIC_SKELETON_UTILITY_CONVOLUTION_HPP


#include "msk_config.h"

#include <omp.h>

#include <algorithm>
#include <complex>
#include <type_traits>
#include <iostream>

#include "msk/utility/common.hpp"

#ifdef MSK_USE_FFT
#include "msk/utility/fourier.hpp"
#endif

namespace MSk {


/**
 * @brief Convolve two vectors using dense matrix as the kernel and add to a vector.
 *
 * @tparam DataType Elements type.
 *
 * @param[in] m Matrix row number.
 * @param[in] n Matrix column number.
 * @param[in] A Matrix array.
 * @param[in] y Left-hand side vector.
 * @param[in] x Right-hand side vector.
 * @param[in, out] z Output vector to add to.
 * @param[in] size Convolution size to compute.
 */
template<typename DataType>
void convolve_add_dense(int m, int n, const DataType *A, const DataType *y, const DataType *x,
                        DataType *z, int size = 0) {
    if (size == 0) {
        size = m + n - 1;
    }
    size = std::min(size, m + n - 1);
#pragma omp parallel for default(shared)
    for (int k = 0; k < size; ++k) {
        // 0 <= i < m
        // 0 <= k-i < n --> k-n+1 <= i < k+1
        const int i_start{std::max(0, k - n + 1)}, i_end{std::min(k + 1, m)};
        for (int i{i_start}; i < i_end; ++i) {
            z[k] += y[i] * A[i + m * (k - i)] * x[k - i];
        }
    }
}

/**
 * @brief Convolve two vectors using dense matrix with hermitian conjugation
 * as the kernel and add to a vector.
 *
 * @tparam DataType Elements type.
 *
 * @param[in] m Matrix row number.
 * @param[in] n Matrix column number.
 * @param[in] A Matrix array.
 * @param[in] y Left-hand side vector.
 * @param[in] x Right-hand side vector.
 * @param[in, out] z Output vector to add to.
 * @param[in] size Convolution size to compute.
 */
template<typename DataType>
void conj_convolve_add_dense(int m, int n, const DataType *A, const DataType *y, const DataType *x,
                             DataType *z, int size = 0) {
    if (size == 0) {
        size = m + n - 1;
    }
    size = std::min(size, m + n - 1);
#pragma omp parallel for default(shared)
    for (int k = 0; k < size; ++k) {
        // 0 <= i < n
        // 0 <= k-i < m --> k-m+1 <= i < k+1
        const int i_start{std::max(0, k - m + 1)}, i_end{std::min(k + 1, n)};
        for (int i{i_start}; i < i_end; ++i) {
            z[k] += y[i] * conjugate(A[(k - i) + m * i]) * x[k - i];
        }
    }
}

/**
 * @brief Convolve two vectors using dense matrix as the kernel and add to a vector.
 *
 * @tparam DataType Elements type.
 *
 * @param[in] m Matrix row number.
 * @param[in] n Matrix column number.
 * @param[in] U Left-hand side factor of the matrix.
 * @param[in] V Right-hand side factor of the matrix.
 * @param[in] y Left-hand side vector.
 * @param[in] x Right-hand side vector.
 * @param[in, out] z Output vector to add to.
 * @param[in] size Convolution size to compute.
 */
#ifdef MSK_USE_FFT
template<typename DataType>
void convolve_add_lowrank(int m, int n, int r, const DataType *U, const DataType *V,
                          const DataType *y, const DataType *x, DataType *z, int size = 0) {
    using SpecType = std::complex<decltype(std::abs(DataType{}))>;
    if (size == 0) {
        size = m + n - 1;
    }
    const int conv_size = m + n - 1;
    int spec_size;
    if (sizeof(DataType) == sizeof(SpecType)) {
        spec_size = conv_size;
    } else {
        spec_size = conv_size / 2 + 1;
    }

    const FFT::Daemon &fft_daemon = FFT::Daemon::get();
    // Pre-allocate plans
    fft_daemon.get_plan<DataType>(conv_size, 2 * r, true);
    fft_daemon.get_plan<DataType>(conv_size, r, false);
    // Columns to be Fourier-transformed
    DataType *UyVx_in{nullptr};
    SpecType *UyVx_out{nullptr};
    int mem_id = fft_daemon.lock_memory(UyVx_in, UyVx_out);
    if (UyVx_in == nullptr || UyVx_out == nullptr) {

    }

    // Prepare data for FFTW: U (*) y
    for (int j = 0; j < r; ++j) {
        for (int i = 0; i < m; ++i) {
            UyVx_in[i + conv_size * j] = U[i + m * j] * y[i];
        }
        for (int i = m; i < conv_size; ++i) {
            UyVx_in[i + conv_size * j] = 0;
        }
    }
    // Prepare data for FFTW: conj(V) (*) x
    for (int j = 0; j < r; ++j) {
        for (int i = 0; i < n; ++i) {
            UyVx_in[i + conv_size * (r + j)] = conjugate(V[i + n * j]) * x[i];
        }
        for (int i = n; i < conv_size; ++i) {
            UyVx_in[i + conv_size * (r + j)] = 0;
        }
    }

    // Perform forward transform
    fft_daemon.transform(conv_size, 2 * r, UyVx_in, UyVx_out, true);

    // Hadamard product in spectral space
    for (int j = 0; j < r; ++j) {
        for (int i = 0; i < spec_size; ++i) {
            UyVx_out[i + spec_size * j] *= UyVx_out[i + spec_size * (j + r)];
        }
    }

    // Perform backward transform
    fft_daemon.transform(conv_size, r, UyVx_in, UyVx_out, false);

    // Accumulate the result
    const DataType inv_factor{DataType{1.0} / static_cast<DataType>(conv_size)};
    for (int j = 0; j < r; ++j) {
        for (int i = 0; i < std::min(conv_size, size); ++i) {
            z[i] += UyVx_in[i + conv_size * j] * inv_factor;
        }
    }

    // Free memory
    fft_daemon.unlock_memory<DataType>(mem_id);
}
#else
template<typename DataType>
void convolve_add_lowrank(int m, int n, int r, const DataType *U, const DataType *V,
                          const DataType *y, const DataType *x, DataType *z, int size = 0) {
    if (size == 0) {
        size = m + n - 1;
    }
    size = std::min(size, m + n - 1);
#pragma omp parallel for default(shared)
    for (int k = 0; k < size; ++k) {
        // 0 <= i < m
        // 0 <= k-i < n --> k-n+1 <= i < k+1
        const int i_start{std::max(0, k - n + 1)}, i_end{std::min(k + 1, m)};
        for (int rank{0}; rank < r; ++rank) {
            for (int i{i_start}; i < i_end; ++i) {
                z[k] += y[i] * U[i + m * rank] * conjugate(V[(k - i) + n * rank]) * x[k - i];
            }
        }
    }
}
#endif



extern template float conjugate(const float &x);

extern template double conjugate(const double &x);

extern template std::complex<float> conjugate(const std::complex<float> &x);

extern template std::complex<double> conjugate(const std::complex<double> &x);


extern template
void convolve_add_dense<float>(int m, int n, const float *A,
const float *y, const float *x, float *z, int size);

extern template
void convolve_add_dense<double>(int m, int n, const double *A,
        const double *y, const double *x, double *z, int size);

extern template
void convolve_add_dense<std::complex<float>>(int m, int n, const std::complex<float> *A,
        const std::complex<float> *y, const std::complex<float> *x, std::complex<float> *z, int size);

extern template
void convolve_add_dense<std::complex<double>>(int m, int n, const std::complex<double> *A,
        const std::complex<double> *y, const std::complex<double> *x, std::complex<double> *z, int size);


extern template
void convolve_add_lowrank<float>(int m, int n, int r, const float *U, const float *V,
        const float *y, const float *x, float *z, int size);

extern template
void convolve_add_lowrank<double>(int m, int n, int r, const double *U, const double *V,
        const double *y, const double *x, double *z, int size);

extern template
void convolve_add_lowrank<std::complex<float>>(int m, int n, int r,
        const std::complex<float> *U, const std::complex<float> *V,
        const std::complex<float> *y, const std::complex<float> *x, std::complex<float> *z, int size);

extern template
void convolve_add_lowrank<std::complex<double>>(int m, int n, int r,
        const std::complex<double> *U, const std::complex<double> *V,
        const std::complex<double> *y, const std::complex<double> *x, std::complex<double> *z, int size);


}

#endif //MOSAIC_SKELETON_UTILITY_CONVOLUTION_HPP
