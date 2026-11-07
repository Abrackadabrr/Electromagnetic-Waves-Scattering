/*
 * Copyright 2024-2026 Bulat Valiakhmetov
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
// Created by b.valiakhmetov on 4/10/24.
//


#include "msk/utility/svd.hpp"


template uint64_t MSk::truncate_singular_values(uint64_t n, const float *S, const comp::Parameters &params, const bool frob);

template uint64_t MSk::truncate_singular_values(uint64_t n, const double *S, const comp::Parameters &params, const bool frob);


template MSk::comp::Result MSk::perform_svd(uint64_t m, uint64_t n, SmartArray<float> &data,
         SmartArray<float> &U, SmartArray<float> &S, SmartArray<float> &V,
         comp::Parameters params, bool frob_truncate);

template MSk::comp::Result MSk::perform_svd(uint64_t m, uint64_t n, SmartArray<double> &data,
         SmartArray<double> &U, SmartArray<double> &S, SmartArray<double> &V,
         comp::Parameters params, bool frob_truncate);

template MSk::comp::Result MSk::perform_svd(uint64_t m, uint64_t n, SmartArray<std::complex<float>> &data,
         SmartArray<std::complex<float>> &U, SmartArray<float> &S, SmartArray<std::complex<float>> &V,
         comp::Parameters params, bool frob_truncate);

template MSk::comp::Result MSk::perform_svd(uint64_t m, uint64_t n, SmartArray<std::complex<double>> &data,
         SmartArray<std::complex<double>> &U, SmartArray<double> &S, SmartArray<std::complex<double>> &V,
         comp::Parameters params, bool frob_truncate);


template MSk::comp::Result MSk::perform_svd(uint64_t m, uint64_t n, SmartArray<float> &data,
         SmartArray<float> &col_factor, SmartArray<float> &row_factor,
         comp::Parameters params);

template MSk::comp::Result MSk::perform_svd(uint64_t m, uint64_t n, SmartArray<double> &data,
         SmartArray<double> &col_factor, SmartArray<double> &row_factor,
         comp::Parameters params);

template MSk::comp::Result MSk::perform_svd(uint64_t m, uint64_t n, SmartArray<std::complex<float>> &data,
         SmartArray<std::complex<float>> &col_factor, SmartArray<std::complex<float>> &row_factor,
         comp::Parameters params);

template MSk::comp::Result MSk::perform_svd(uint64_t m, uint64_t n, SmartArray<std::complex<double>> &data,
         SmartArray<std::complex<double>> &col_factor, SmartArray<std::complex<double>> &row_factor,
         comp::Parameters params);


template MSk::comp::Result MSk::recompress_lowrank(uint64_t m, uint64_t n, uint64_t rank,
        SmartArray<float> &col_factor, SmartArray<float> &row_factor,
        const comp::Parameters &params, const bool orth_cols, const bool orth_rows);

template MSk::comp::Result MSk::recompress_lowrank(uint64_t m, uint64_t n, uint64_t rank,
        SmartArray<double> &col_factor, SmartArray<double> &row_factor,
        const comp::Parameters &params, const bool orth_cols, const bool orth_rows);

template MSk::comp::Result MSk::recompress_lowrank(uint64_t m, uint64_t n, uint64_t rank,
        SmartArray<std::complex<float>> &col_factor, SmartArray<std::complex<float>> &row_factor,
        const comp::Parameters &params, const bool orth_cols, const bool orth_rows);

template MSk::comp::Result MSk::recompress_lowrank(uint64_t m, uint64_t n, uint64_t rank,
        SmartArray<std::complex<double>> &col_factor, SmartArray<std::complex<double>> &row_factor,
        const comp::Parameters &params, const bool orth_cols, const bool orth_rows);
