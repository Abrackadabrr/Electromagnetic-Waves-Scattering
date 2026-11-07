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
// Created by b.valiakhmetov on 9/6/23.
//


/**
 * @file
 *
 * @brief Contains C language side of C/C++ interface to class MosaicSkeleton.
 */


#ifndef MOSAIC_SKELETON_INTERFACE_C_MSK_H
#define MOSAIC_SKELETON_INTERFACE_C_MSK_H


#include "msk_config.h"

#ifdef MSK_USE_MPI
#include <mpi.h>
#endif
#include <omp.h>
#include <stdbool.h>
#include <stdint.h>
#include <complex.h>

#undef complex
#undef I
#undef J


#ifdef __cplusplus
extern "C"
{
#endif

typedef struct mskMosaicMatrix {
    void *_matrix;

    uint64_t row_size, col_size;
    double frob_norm;
    uint64_t memory;

#ifdef MSK_USE_MPI
    uint64_t local_row_size, local_col_size;
    uint64_t local_row_off, local_col_off;

    double local_frob_norm;

    uint64_t local_memory;
#endif
} mskMosaicMatrix;


typedef struct mskParameters {
    uint64_t min_rows, min_cols;
    double rho;
    uint64_t sphere_norm;
    bool force_lowrank;

    double abs_tol, rel_tol;
    uint64_t max_rank;

    const char *log_file;
    const char *log_level;
} mskParameters;


typedef struct mskSolveParameters {
    double abs_tol, rel_tol;
    uint64_t iter_limit;
} mskSolveParameters;


int msk_build_basic_s3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t n,
                        float (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const float *row_pts, const float *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
                        , MPI_Comm comm
#endif
);

int msk_build_basic_d3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t n,
                        double (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const double *row_pts, const double *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
);

int msk_build_basic_c3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t n,
                        float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const float *row_pts, const float *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
);

int msk_build_basic_z3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t n,
                        double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const double *row_pts, const double *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
);


int msk_build_block_s3(mskMosaicMatrix *msk_matrix,
                    mskParameters *parameters,
                    uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                    uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                    float (*f_elem)(uint64_t, uint64_t, void *), void *data,
                    const float *row_pts, const float *col_pts,
                    uint64_t n_jobs
#ifdef MSK_USE_MPI
                    , MPI_Comm comm
#endif
);

int msk_build_block_d3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                        double (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const double *row_pts, const double *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
);

int msk_build_block_c3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                        float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const float *row_pts, const float *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
);

int msk_build_block_z3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                        double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const double *row_pts, const double *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
);


int msk_recompress_s3(mskMosaicMatrix *msk_matrix, double rel_tol);

int msk_recompress_d3(mskMosaicMatrix *msk_matrix, double rel_tol);

int msk_recompress_c3(mskMosaicMatrix *msk_matrix, double rel_tol);

int msk_recompress_z3(mskMosaicMatrix *msk_matrix, double rel_tol);


int msk_destroy_s3(mskMosaicMatrix *msk_matrix);

int msk_destroy_d3(mskMosaicMatrix *msk_matrix);

int msk_destroy_c3(mskMosaicMatrix *msk_matrix);

int msk_destroy_z3(mskMosaicMatrix *msk_matrix);


int msk_matvec_s3(const mskMosaicMatrix *msk_matrix, const char trans, const float *x, float *y);

int msk_matvec_d3(const mskMosaicMatrix *msk_matrix, const char trans, const double *x, double *y);

int msk_matvec_c3(const mskMosaicMatrix *msk_matrix, const char trans, const float _Complex *x, float _Complex *y);

int msk_matvec_z3(const mskMosaicMatrix *msk_matrix, const char trans, const double _Complex *x, double _Complex *y);


#ifdef MSK_USE_MPI
int msk_to_dense_s3(const mskMosaicMatrix *msk_matrix, float *ptr, uint64_t ldim, int root);

int msk_to_dense_d3(const mskMosaicMatrix *msk_matrix, double *ptr, uint64_t ldim, int root);

int msk_to_dense_c3(const mskMosaicMatrix *msk_matrix, float _Complex *ptr, uint64_t ldim, int root);

int msk_to_dense_z3(const mskMosaicMatrix *msk_matrix, double _Complex *ptr, uint64_t ldim, int root);
#else
int msk_to_dense_s3(const mskMosaicMatrix *msk_matrix, float *ptr, uint64_t ldim);

int msk_to_dense_d3(const mskMosaicMatrix *msk_matrix, double *ptr, uint64_t ldim);

int msk_to_dense_c3(const mskMosaicMatrix *msk_matrix, float _Complex *ptr, uint64_t ldim);

int msk_to_dense_z3(const mskMosaicMatrix *msk_matrix, double _Complex *ptr, uint64_t ldim);
#endif


#ifdef MSK_USE_MPI
int msk_permute_col_s3(const mskMosaicMatrix *msk_matrix, const float *x_in, float *x_out, bool global);

int msk_permute_col_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out, bool global);

int msk_permute_col_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out, bool global);

int msk_permute_col_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out, bool global);


int msk_permute_col_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x, bool global);

int msk_permute_col_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x, bool global);

int msk_permute_col_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x, bool global);

int msk_permute_col_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x, bool global);


int msk_unpermute_col_s3(const mskMosaicMatrix *msk_matrix, const float *x_in, float *x_out, bool global);

int msk_unpermute_col_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out, bool global);

int msk_unpermute_col_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out, bool global);

int msk_unpermute_col_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out, bool global);


int msk_unpermute_col_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x, bool global);

int msk_unpermute_col_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x, bool global);

int msk_unpermute_col_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x, bool global);

int msk_unpermute_col_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x, bool global);


int msk_permute_row_s3(const mskMosaicMatrix *msk_matrix, const float *x_in, float *x_out, bool global);

int msk_permute_row_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out, bool global);

int msk_permute_row_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out, bool global);

int msk_permute_row_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out, bool global);


int msk_permute_row_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x, bool global);

int msk_permute_row_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x, bool global);

int msk_permute_row_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x, bool global);

int msk_permute_row_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x, bool global);


int msk_unpermute_row_s3(const mskMosaicMatrix *msk_matrix, const float *x_in, float *x_out, bool global);

int msk_unpermute_row_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out, bool global);

int msk_unpermute_row_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out, bool global);

int msk_unpermute_row_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out, bool global);


int msk_unpermute_row_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x, bool global);

int msk_unpermute_row_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x, bool global);

int msk_unpermute_row_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x, bool global);

int msk_unpermute_row_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x, bool global);
#else
int msk_permute_col_s3(const mskMosaicMatrix *msk_matrix, const float *x_in, float *x_out);

int msk_permute_col_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out);

int msk_permute_col_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out);

int msk_permute_col_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out);

int msk_permute_col_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x);

int msk_permute_col_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x);

int msk_permute_col_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x);

int msk_permute_col_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x);


int msk_unpermute_col_s3(const mskMosaicMatrix *msk_matrix, const float *x_in, float *x_out);

int msk_unpermute_col_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out);

int msk_unpermute_col_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out);

int msk_unpermute_col_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out);


int msk_unpermute_col_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x);

int msk_unpermute_col_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x);

int msk_unpermute_col_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x);

int msk_unpermute_col_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x);


int msk_permute_row_s3(const mskMosaicMatrix *msk_matrix, const float *x_in, float *x_out);

int msk_permute_row_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out);

int msk_permute_row_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out);

int msk_permute_row_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out);


int msk_permute_row_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x);

int msk_permute_row_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x);

int msk_permute_row_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x);

int msk_permute_row_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x);


int msk_unpermute_row_s3(const mskMosaicMatrix *msk_matrix, const float *x_in, float *x_out);

int msk_unpermute_row_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out);

int msk_unpermute_row_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out);

int msk_unpermute_row_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out);


int msk_unpermute_row_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x);

int msk_unpermute_row_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x);

int msk_unpermute_row_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x);

int msk_unpermute_row_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x);
#endif


#ifdef MSK_USE_MPI
int msk_check_error_s3(mskMosaicMatrix *msk_matrix, float (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, bool global, double *error);

int msk_check_error_d3(mskMosaicMatrix *msk_matrix, double (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, bool global, double *error);

int msk_check_error_c3(mskMosaicMatrix *msk_matrix, float _Complex (*f_elem)(uint64_t, uint64_t, void *),void *data,
                        bool relative, bool global, double *error);

int msk_check_error_z3(mskMosaicMatrix *msk_matrix, double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, bool global, double *error);
#else
int msk_check_error_s3(mskMosaicMatrix *msk_matrix, float (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, double *error);

int msk_check_error_d3(mskMosaicMatrix *msk_matrix, double (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, double *error);

int msk_check_error_c3(mskMosaicMatrix *msk_matrix, float _Complex (*f_elem)(uint64_t, uint64_t, void *),void *data,
                        bool relative, double *error);

int msk_check_error_z3(mskMosaicMatrix *msk_matrix, double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, double *error);
#endif


int msk_approx_block_solve_s3(uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                               uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                               float (*f_elem)(uint64_t, uint64_t, void *), void *data,
                               const float *row_pts, const float *col_pts,
                               mskParameters *parameters,
                               uint64_t n_rhs, float *B, float *X,
                               mskSolveParameters *solve_params,
                               uint64_t n_jobs
#ifdef MSK_USE_MPI
                               , MPI_Comm comm
#endif
);

int msk_approx_block_solve_d3(uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                               uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                               double (*f_elem)(uint64_t, uint64_t, void *), void *data,
                               const double *row_pts, const double *col_pts,
                               mskParameters *parameters,
                               uint64_t n_rhs, double *B, double *X,
                               mskSolveParameters *solve_params,
                               uint64_t n_jobs
#ifdef MSK_USE_MPI
                               , MPI_Comm comm
#endif
);

int msk_approx_block_solve_c3(uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                               uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                               float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                               const float *row_pts, const float *col_pts,
                               mskParameters *parameters,
                               uint64_t n_rhs, float _Complex *B, float _Complex *X,
                               mskSolveParameters *solve_params,
                               uint64_t n_jobs
#ifdef MSK_USE_MPI
                               , MPI_Comm comm
#endif
);

int msk_approx_block_solve_z3(uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                               uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                               double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                               const double *row_pts, const double *col_pts,
                               mskParameters *parameters,
                               uint64_t n_rhs, double _Complex *B, double _Complex *X,
                               mskSolveParameters *solve_params,
                               uint64_t n_jobs
#ifdef MSK_USE_MPI
                               , MPI_Comm comm
#endif
);


#ifdef __cplusplus
}
#endif




#endif // MOSAIC_SKELETON_INTERFACE_C_MSK_H
