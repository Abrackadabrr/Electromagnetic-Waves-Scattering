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
// Created by b.valiakhmetov on 4/9/24.
//


#include "msk/utility/convolution.hpp"


template
float MSk::conjugate(const float &x);

template
double MSk::conjugate(const double &x);

template
std::complex<float> MSk::conjugate(const std::complex<float> &x);

template
std::complex<double> MSk::conjugate(const std::complex<double> &x);


template
void MSk::convolve_add_dense(int m, int n, const float *A,
        const float *y, const float *x, float *z, int size);

template
void MSk::convolve_add_dense(int m, int n, const double *A,
        const double *y, const double *x, double *z, int size);

template
void MSk::convolve_add_dense(int m, int n, const std::complex<float> *A,
        const std::complex<float> *y, const std::complex<float> *x, std::complex<float> *z, int size);

template
void MSk::convolve_add_dense(int m, int n, const std::complex<double> *A,
        const std::complex<double> *y, const std::complex<double> *x, std::complex<double> *z, int size);


template
void MSk::convolve_add_lowrank(int m, int n, int r, const float *U, const float *V,
        const float *y, const float *x, float *z, int size);

template
void MSk::convolve_add_lowrank(int m, int n, int r, const double *U, const double *V,
        const double *y, const double *x, double *z, int size);

template
void MSk::convolve_add_lowrank(int m, int n, int r,
        const std::complex<float> *U, const std::complex<float> *V,
        const std::complex<float> *y, const std::complex<float> *x, std::complex<float> *z, int size);

template
void MSk::convolve_add_lowrank(int m, int n, int r,
        const std::complex<double> *U, const std::complex<double> *V,
        const std::complex<double> *y, const std::complex<double> *x, std::complex<double> *z, int size);
