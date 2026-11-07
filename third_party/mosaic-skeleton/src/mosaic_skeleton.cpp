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


#include "msk/mosaic_skeleton.hpp"


template class MSk::MosaicSkeleton<float>;

template class MSk::MosaicSkeleton<double>;

template class MSk::MosaicSkeleton<std::complex<float>>;

template class MSk::MosaicSkeleton<std::complex<double>>;


template MSk::MosaicSkeleton<float> MSk::MosaicSkeleton<float>::approximate(matrix::Elementary<float> &matrix,
        oracle::Elementary<float> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<float> MSk::MosaicSkeleton<float>::approximate(matrix::Dense<float, 2> &matrix,
        oracle::Distant<float, matrix::Dense<float, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<float> MSk::MosaicSkeleton<float>::approximate(matrix::Dense<float, 3> &matrix,
        oracle::Distant<float, matrix::Dense<float, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<float> MSk::MosaicSkeleton<float>::approximate(matrix::BiPoint<float, float, 2> &matrix,
        oracle::Distant<float, matrix::BiPoint<float, float, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<float> MSk::MosaicSkeleton<float>::approximate(matrix::BiPoint<float, float, 3> &matrix,
        oracle::Distant<float, matrix::BiPoint<float, float, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);


template MSk::MosaicSkeleton<double> MSk::MosaicSkeleton<double>::approximate(matrix::Elementary<double> &matrix,
        oracle::Elementary<double> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<double> MSk::MosaicSkeleton<double>::approximate(matrix::Dense<double, 2> &matrix,
        oracle::Distant<double, matrix::Dense<double, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<double> MSk::MosaicSkeleton<double>::approximate(matrix::Dense<double, 3> &matrix,
        oracle::Distant<double, matrix::Dense<double, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<double> MSk::MosaicSkeleton<double>::approximate(matrix::BiPoint<double, double, 2> &matrix,
        oracle::Distant<double, matrix::BiPoint<double, double, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<double> MSk::MosaicSkeleton<double>::approximate(matrix::BiPoint<double, double, 3> &matrix,
        oracle::Distant<double, matrix::BiPoint<double, double, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);
    

template MSk::MosaicSkeleton<std::complex<float>> MSk::MosaicSkeleton<std::complex<float>>::approximate(matrix::Elementary<std::complex<float>> &matrix,
        oracle::Elementary<std::complex<float>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<std::complex<float>> MSk::MosaicSkeleton<std::complex<float>>::approximate(matrix::Dense<std::complex<float>, 2> &matrix,
        oracle::Distant<std::complex<float>, matrix::Dense<std::complex<float>, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<std::complex<float>> MSk::MosaicSkeleton<std::complex<float>>::approximate(matrix::Dense<std::complex<float>, 3> &matrix,
        oracle::Distant<std::complex<float>, matrix::Dense<std::complex<float>, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<std::complex<float>> MSk::MosaicSkeleton<std::complex<float>>::approximate(matrix::BiPoint<std::complex<float>, float, 2> &matrix,
        oracle::Distant<std::complex<float>, matrix::BiPoint<std::complex<float>, float, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<std::complex<float>> MSk::MosaicSkeleton<std::complex<float>>::approximate(matrix::BiPoint<std::complex<float>, float, 3> &matrix,
        oracle::Distant<std::complex<float>, matrix::BiPoint<std::complex<float>, float, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

    
template MSk::MosaicSkeleton<std::complex<double>> MSk::MosaicSkeleton<std::complex<double>>::approximate(matrix::Elementary<std::complex<double>> &matrix,
        oracle::Elementary<std::complex<double>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<std::complex<double>> MSk::MosaicSkeleton<std::complex<double>>::approximate(matrix::Dense<std::complex<double>, 2> &matrix,
        oracle::Distant<std::complex<double>, matrix::Dense<std::complex<double>, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<std::complex<double>> MSk::MosaicSkeleton<std::complex<double>>::approximate(matrix::Dense<std::complex<double>, 3> &matrix,
        oracle::Distant<std::complex<double>, matrix::Dense<std::complex<double>, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<std::complex<double>> MSk::MosaicSkeleton<std::complex<double>>::approximate(matrix::BiPoint<std::complex<double>, double, 2> &matrix,
        oracle::Distant<std::complex<double>, matrix::BiPoint<std::complex<double>, double, 2>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);

template MSk::MosaicSkeleton<std::complex<double>> MSk::MosaicSkeleton<std::complex<double>>::approximate(matrix::BiPoint<std::complex<double>, double, 3> &matrix,
        oracle::Distant<std::complex<double>, matrix::BiPoint<std::complex<double>, double, 3>> &oracle, uint64_t n_jobs,
#ifdef MSK_USE_MPI
        MPI_Comm comm,
#endif
        std::shared_ptr<MSk::log::Logger> logger_ptr);


template void MSk::MosaicSkeleton<float>::build_mosaic(oracle::Elementary<float> &oracle);
template void MSk::MosaicSkeleton<float>::build_mosaic(oracle::Distant<float, matrix::Dense<float, 2>> &oracle);
template void MSk::MosaicSkeleton<float>::build_mosaic(oracle::Distant<float, matrix::Dense<float, 3>> &oracle);
template void MSk::MosaicSkeleton<float>::build_mosaic(oracle::Distant<float, matrix::BiPoint<float, float, 2>> &oracle);
template void MSk::MosaicSkeleton<float>::build_mosaic(oracle::Distant<float, matrix::BiPoint<float, float, 3>> &oracle);

template void MSk::MosaicSkeleton<float>::approximate(matrix::Elementary<float> &matrix);
template void MSk::MosaicSkeleton<float>::approximate(matrix::Dense<float, 2> &matrix);
template void MSk::MosaicSkeleton<float>::approximate(matrix::Dense<float, 3> &matrix);
template void MSk::MosaicSkeleton<float>::approximate(matrix::BiPoint<float, float, 2> &matrix);
template void MSk::MosaicSkeleton<float>::approximate(matrix::BiPoint<float, float, 3> &matrix);


template void MSk::MosaicSkeleton<double>::build_mosaic(oracle::Elementary<double> &oracle);
template void MSk::MosaicSkeleton<double>::build_mosaic(oracle::Distant<double, matrix::Dense<double, 2>> &oracle);
template void MSk::MosaicSkeleton<double>::build_mosaic(oracle::Distant<double, matrix::Dense<double, 3>> &oracle);
template void MSk::MosaicSkeleton<double>::build_mosaic(oracle::Distant<double, matrix::BiPoint<double, double, 2>> &oracle);
template void MSk::MosaicSkeleton<double>::build_mosaic(oracle::Distant<double, matrix::BiPoint<double, double, 3>> &oracle);

template void MSk::MosaicSkeleton<double>::approximate(matrix::Elementary<double> &matrix);
template void MSk::MosaicSkeleton<double>::approximate(matrix::Dense<double, 2> &matrix);
template void MSk::MosaicSkeleton<double>::approximate(matrix::Dense<double, 3> &matrix);
template void MSk::MosaicSkeleton<double>::approximate(matrix::BiPoint<double, double, 2> &matrix);
template void MSk::MosaicSkeleton<double>::approximate(matrix::BiPoint<double, double, 3> &matrix);


template void MSk::MosaicSkeleton<std::complex<float>>::build_mosaic(oracle::Elementary<std::complex<float>> &oracle);
template void MSk::MosaicSkeleton<std::complex<float>>::build_mosaic(oracle::Distant<std::complex<float>, matrix::Dense<std::complex<float>, 2>> &oracle);
template void MSk::MosaicSkeleton<std::complex<float>>::build_mosaic(oracle::Distant<std::complex<float>, matrix::Dense<std::complex<float>, 3>> &oracle);
template void MSk::MosaicSkeleton<std::complex<float>>::build_mosaic(oracle::Distant<std::complex<float>, matrix::BiPoint<std::complex<float>, float, 2>> &oracle);
template void MSk::MosaicSkeleton<std::complex<float>>::build_mosaic(oracle::Distant<std::complex<float>, matrix::BiPoint<std::complex<float>, float, 3>> &oracle);

template void MSk::MosaicSkeleton<std::complex<float>>::approximate(matrix::Elementary<std::complex<float>> &matrix);
template void MSk::MosaicSkeleton<std::complex<float>>::approximate(matrix::Dense<std::complex<float>, 2> &matrix);
template void MSk::MosaicSkeleton<std::complex<float>>::approximate(matrix::Dense<std::complex<float>, 3> &matrix);
template void MSk::MosaicSkeleton<std::complex<float>>::approximate(matrix::BiPoint<std::complex<float>, float, 2> &matrix);
template void MSk::MosaicSkeleton<std::complex<float>>::approximate(matrix::BiPoint<std::complex<float>, float, 3> &matrix);


template void MSk::MosaicSkeleton<std::complex<double>>::build_mosaic(oracle::Elementary<std::complex<double>> &oracle);
template void MSk::MosaicSkeleton<std::complex<double>>::build_mosaic(oracle::Distant<std::complex<double>, matrix::Dense<std::complex<double>, 2>> &oracle);
template void MSk::MosaicSkeleton<std::complex<double>>::build_mosaic(oracle::Distant<std::complex<double>, matrix::Dense<std::complex<double>, 3>> &oracle);
template void MSk::MosaicSkeleton<std::complex<double>>::build_mosaic(oracle::Distant<std::complex<double>, matrix::BiPoint<std::complex<double>, double, 2>> &oracle);
template void MSk::MosaicSkeleton<std::complex<double>>::build_mosaic(oracle::Distant<std::complex<double>, matrix::BiPoint<std::complex<double>, double, 3>> &oracle);

template void MSk::MosaicSkeleton<std::complex<double>>::approximate(matrix::Elementary<std::complex<double>> &matrix);
template void MSk::MosaicSkeleton<std::complex<double>>::approximate(matrix::Dense<std::complex<double>, 2> &matrix);
template void MSk::MosaicSkeleton<std::complex<double>>::approximate(matrix::Dense<std::complex<double>, 3> &matrix);
template void MSk::MosaicSkeleton<std::complex<double>>::approximate(matrix::BiPoint<std::complex<double>, double, 2> &matrix);
template void MSk::MosaicSkeleton<std::complex<double>>::approximate(matrix::BiPoint<std::complex<double>, double, 3> &matrix);
