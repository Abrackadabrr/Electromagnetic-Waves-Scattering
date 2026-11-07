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


#include "msk/utility/testing.hpp"



template void MSk::run_gmres(const MSk::MosaicSkeleton<float> &msk_matrix,
       uint64_t n_rhs, float *rhs, uint64_t ld_rhs, float *sol, uint64_t ld_sol,
       MSk::solvers::linsys::Parameters<float> params, std::shared_ptr<MSk::log::Logger> log);

template void MSk::run_gmres(const MSk::MosaicSkeleton<double> &msk_matrix,
       uint64_t n_rhs, double *rhs, uint64_t ld_rhs, double *sol, uint64_t ld_sol,
       MSk::solvers::linsys::Parameters<double> params, std::shared_ptr<MSk::log::Logger> log);

template void MSk::run_gmres(const MSk::MosaicSkeleton<std::complex<float>> &msk_matrix,
       uint64_t n_rhs, std::complex<float> *rhs, uint64_t ld_rhs, std::complex<float> *sol, uint64_t ld_sol,
       MSk::solvers::linsys::Parameters<float> params, std::shared_ptr<MSk::log::Logger> log);

template void MSk::run_gmres(const MSk::MosaicSkeleton<std::complex<double>> &msk_matrix,
       uint64_t n_rhs, std::complex<double> *rhs, uint64_t ld_rhs, std::complex<double> *sol, uint64_t ld_sol,
       MSk::solvers::linsys::Parameters<double> params, std::shared_ptr<MSk::log::Logger> log);


template uint64_t MSk::block_svd(uint64_t n_blocks, const std::vector<const MSk::block::LowRank<float> *> &blocks,
       SmartArray<float> &U, SmartArray<float> &S, SmartArray<float> &V,
       float tol, bool frob);

template uint64_t MSk::block_svd(uint64_t n_blocks, const std::vector<const MSk::block::LowRank<double> *> &blocks,
       SmartArray<double> &U, SmartArray<double> &S, SmartArray<double> &V,
       double tol, bool frob);

template uint64_t MSk::block_svd(uint64_t n_blocks, const std::vector<const MSk::block::LowRank<std::complex<float>> *> &blocks,
       SmartArray<std::complex<float>> &U, SmartArray<float> &S, SmartArray<std::complex<float>> &V,
       float tol, bool frob);

template uint64_t MSk::block_svd(uint64_t n_blocks, const std::vector<const MSk::block::LowRank<std::complex<double>> *> &blocks,
       SmartArray<std::complex<double>> &U, SmartArray<double> &S, SmartArray<std::complex<double>> &V,
       double tol, bool frob);
