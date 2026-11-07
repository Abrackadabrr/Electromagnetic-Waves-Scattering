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
 *
 * @brief Contains C language side of C/C++ interface to class MosaicSkeleton.
 */


#ifndef MOSAIC_SKELETON_MSK_INTERFACE_FORTRAN_WRAP_BUILD_H
#define MOSAIC_SKELETON_MSK_INTERFACE_FORTRAN_WRAP_BUILD_H


#include "msk/interface/C/msk.h"


int msk_build_fort_s3(
        mskMosaicMatrix *msk_matrix,
        mskParameters *parameters,
        uint64_t m, uint64_t n,
        float (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const float *row_pts, const float *col_pts,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);

int msk_build_fort_d3(
        mskMosaicMatrix *msk_matrix,
        mskParameters *parameters,
        uint64_t m, uint64_t n,
        double (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const double *row_pts, const double *col_pts,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);

int msk_build_fort_c3(
        mskMosaicMatrix *msk_matrix,
        mskParameters *parameters,
        uint64_t m, uint64_t n,
        float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const float *row_pts, const float *col_pts,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);

int msk_build_fort_z3(
        mskMosaicMatrix *msk_matrix,
        mskParameters *parameters,
        uint64_t m, uint64_t n,
        double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const double *row_pts, const double *col_pts,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);


int msk_build_block_fort_s3(
        mskMosaicMatrix *msk_matrix,
        mskParameters *parameters,
        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
        float (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const float *row_pts, const float *col_pts,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);

int msk_build_block_fort_d3(
        mskMosaicMatrix *msk_matrix,
        mskParameters *parameters,
        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
        double (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const double *row_pts, const double *col_pts,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);

int msk_build_block_fort_c3(
        mskMosaicMatrix *msk_matrix,
        mskParameters *parameters,
        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
        float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const float *row_pts, const float *col_pts,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);

int msk_build_block_fort_z3(
        mskMosaicMatrix *msk_matrix,
        mskParameters *parameters,
        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
        double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const double *row_pts, const double *col_pts,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);

int msk_approx_block_solve_fort_s3(
        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
        float (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const float *row_pts, const float *col_pts,
        mskParameters *parameters,
        uint64_t n_rhs, float *B, float *X,
        mskSolveParameters *solve_params,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);


int msk_approx_block_solve_fort_d3(
        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
        double (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const double *row_pts, const double *col_pts,
        mskParameters *parameters,
        uint64_t n_rhs, double *B, double *X,
        mskSolveParameters *solve_params,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);

int msk_approx_block_solve_fort_c3(
        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
        float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const float *row_pts, const float *col_pts,
        mskParameters *parameters,
        uint64_t n_rhs, float _Complex *B, float _Complex *X,
        mskSolveParameters *solve_params,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);

int msk_approx_block_solve_fort_z3(
        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
        double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
        const double *row_pts, const double *col_pts,
        mskParameters *parameters,
        uint64_t n_rhs, double _Complex *B, double _Complex *X,
        mskSolveParameters *solve_params,
        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Fint f_comm
#endif
);

#endif // MOSAIC_SKELETON_MSK_INTERFACE_FORTRAN_WRAP_BUILD_H
