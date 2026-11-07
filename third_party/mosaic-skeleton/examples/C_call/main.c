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
// Created by bulat_v on 04.02.2022.
//


#ifdef MSK_EXAMPLES_USE_MPI
#include <mpi.h>
#endif
#include <stdint.h>
#include <complex.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>

#undef complex
#undef I
#undef J


#include "msk_examples_config.h"
#include "msk/interface/C/msk.h"


typedef double _Complex DataType;
typedef double RealType;

const DataType PROBLEM_WAVE_NUMBER = 1.0;
const uint64_t PROBLEM_DIM = 3;
const uint64_t PROBLEM_SIZE = 4096;


typedef struct Point3d {
    RealType x;
    RealType y;
    RealType z;
} Point3d;


Point3d *make_line(uint64_t n) {
    if (n <= 0) {
        return NULL;
    }
    Point3d *ptr = malloc(n * sizeof(ptr[0]));
    for (int i = 0; i < n; ++i) {
        ptr[i].x = ((RealType)i) / ((RealType)n);
        ptr[i].y = 0;
        ptr[i].z = 0;
    }
    return ptr;
}

DataType helmgoltz(uint64_t i, uint64_t j, void *data) {
    if (i == j) {
        return 0;
    }
    Point3d *points = data;
    Point3d p1 = points[i], p2 = points[j];
    p1.x -= p2.x;
    p1.y -= p2.y;
    p1.z -= p2.z;
    RealType r = sqrt(p1.x * p1.x + p1.y * p1.y + p1.z * p1.z);
    return cexp(_Complex_I * PROBLEM_WAVE_NUMBER * r) / (4 * M_PI * r);
}

void densify(uint64_t m, uint64_t n, DataType *matrix,
             DataType (*f_elem)(uint64_t, uint64_t, void *), void *data) {
    for (uint64_t j = 0; j < n; ++j) {
        for (uint64_t i = 0; i < m; ++i) {
            matrix[i + m * j] = f_elem(i, j, data);
        }
    }
}

void dense_matvec(uint64_t m, uint64_t n, const DataType *matrix, const char trans, const DataType *x, DataType *y) {
    if (trans == 'N' || trans == 'n') {
        for (uint64_t i = 0; i < m; ++i) {
            y[i] = 0;
            for (uint64_t j = 0; j < n; ++j) {
                y[i] += matrix[i + m * j] * x[j];
            }
        }
    } else if (trans == 'T' || trans == 't') {
        for (uint64_t j = 0; j < n; ++j) {
            y[j] = 0;
            for (uint64_t i = 0; i < m; ++i) {
                y[j] += matrix[i + m * j] * x[i];
            }
        }
    } else if (trans == 'C' || trans == 'c') {
        for (uint64_t j = 0; j < n; ++j) {
            y[j] = 0;
            for (uint64_t i = 0; i < m; ++i) {
                y[j] += conj(matrix[i + m * j]) * x[i];
            }
        }
    }
}

RealType norm(uint64_t n, const DataType *x) {
    RealType nrm = 0;
    for (uint64_t i = 0; i < n; ++i) {
        nrm += cabs(x[i]) * cabs(x[i]);
    }
    return sqrt(nrm);
}

void print_matrix(uint64_t m, uint64_t n, const DataType *matrix) {
    for (uint64_t i = 0; i < m; ++i) {
        for (uint64_t j = 0; j < n; ++j) {
            printf("(%f, %f) ", creal(matrix[i + m * j]), cimag(matrix[i + m * j]));
        }
        printf("\n");
    }
}
void print_vector(uint64_t n, const DataType *vector) {
    for (uint64_t i = 0; i < n; ++i) {
        printf("(%f, %f) ", creal(vector[i]), cimag(vector[i]));
    }
    printf("\n");
}


int main(int argc, char *argv[]) {
#ifdef MSK_EXAMPLES_USE_MPI
    int mpi_thread_req = MPI_THREAD_MULTIPLE, mpi_thread_prov = 0;
    MPI_Init_thread(&argc, &argv, mpi_thread_req, &mpi_thread_prov);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
#endif

    int proc_id = 0;
#ifdef MSK_EXAMPLES_USE_MPI
    MPI_Comm_rank(comm, &proc_id);
#endif

    // Build geometry
    uint64_t n_pts = PROBLEM_SIZE;
    Point3d *points = make_line(n_pts);

    // Gain parameters
    uint64_t min_rows = 128, min_cols = 128;
    double rho = 1.01;
    uint64_t sphere_norm = 0;
    bool force_lowrank = false;
    double abs_tol = 0, rel_tol = 1e-4;
    uint64_t max_rank = 0;
    const char *log_file = "stdout";
    const char *log_level = "info"; // "off"
    mskParameters params = {min_rows, min_cols, rho, sphere_norm, force_lowrank,
                            abs_tol, rel_tol, max_rank, log_file, log_level};

    // Approximate
    mskMosaicMatrix msk_matrix;
    /*
    msk_build_basic_z3(&msk_matrix, &params, n_pts, n_pts, helmgoltz, points,
                       (RealType *)points, (RealType *)points, omp_get_max_threads()
#ifdef MSK_EXAMPLES_USE_MPI
                       , MPI_COMM_WORLD
#endif
    );
    */
    uint64_t blocks[2];
    blocks[0] = n_pts / 3;
    blocks[1] = n_pts - n_pts / 3;

    // Solve linear system using mosaic approximation
    mskSolveParameters solve_params;
    solve_params.abs_tol = 0;
    solve_params.rel_tol = 1e-6;
    solve_params.iter_limit = 0;


    DataType *rhs = malloc(n_pts * sizeof(rhs[0]));
    DataType *solution = malloc(n_pts * sizeof(solution[0]));
    for (int i = 0; i < n_pts; ++i) {
        rhs[i] = 1;
        solution[i] = 0;
    }

    msk_approx_block_solve_z3(n_pts, 2, blocks, n_pts, 2, blocks,
                              helmgoltz, points, (RealType *)points, (RealType *)points, &params,
                              1, rhs, solution,
                              &solve_params,
                              omp_get_max_threads()
#ifdef MSK_EXAMPLES_USE_MPI
                              , MPI_COMM_WORLD
#endif
    );

    // Build new mosaic matrix
    msk_build_block_z3(&msk_matrix, &params,
            n_pts, 2, blocks,
            n_pts, 2, blocks, helmgoltz, points,
            (RealType *)points, (RealType *)points, omp_get_max_threads()
#ifdef MSK_EXAMPLES_USE_MPI
            , MPI_COMM_WORLD
#endif
    );

    // Recompress
    RealType recomp_tol = 1e-4;
    msk_recompress_z3(&msk_matrix, recomp_tol);

    // True dense matrix
    DataType *dense_matrix = NULL;
    if (proc_id == 0) {
        dense_matrix = malloc(n_pts * n_pts * sizeof(dense_matrix[0]));
        densify(n_pts, n_pts, dense_matrix, helmgoltz, points);
    }

    // Check approximation error
    DataType *msk_dense = NULL;
    if (proc_id == 0) {
        msk_dense = malloc(n_pts * n_pts * sizeof(msk_dense[0]));
    }
#ifdef MSK_EXAMPLES_USE_MPI
    msk_to_dense_z3(&msk_matrix, msk_dense, n_pts, 0);
#else
    msk_to_dense_z3(&msk_matrix, msk_dense, n_pts);
#endif
    if (proc_id == 0) {
        for (uint64_t i = 0; i < n_pts * n_pts; ++i) {
            msk_dense[i] -= dense_matrix[i];
        }
        RealType true_norm, err_norm;
        true_norm = norm(n_pts * n_pts, dense_matrix);
        err_norm = norm(n_pts * n_pts, msk_dense);
        printf("Approximation error = %.3e\n", err_norm / true_norm);
    }
    double error;
#ifdef MSK_EXAMPLES_USE_MPI
    msk_check_error_z3(&msk_matrix, helmgoltz, points, true, true, &error);
#else
    msk_check_error_z3(&msk_matrix, helmgoltz, points, true, &error);
#endif
    if (proc_id == 0) {
        printf("Checked error = %.3e\n", error);
    }

    // Check matvec error
    DataType *x = malloc(n_pts * sizeof(x[0]));
    DataType *z = malloc(n_pts * sizeof(z[0]));

    // Prepare x and z = Ax and distribute them
    if (proc_id == 0) {
        for (int i = 0; i < n_pts; ++i) {
            x[i] = ((DataType) random()) / ((DataType) RAND_MAX) * 2 - 1;
            z[i] = 0;
        }
        RealType x_norm = norm(n_pts, x);
        for (int i = 0; i < n_pts; ++i) {
            x[i] /= x_norm;
        }
        dense_matvec(n_pts, n_pts, dense_matrix, 'N', x, z);
#ifdef MSK_EXAMPLES_USE_MPI
        MPI_Bcast(&x[0], n_pts * sizeof(x[0]), MPI_UINT8_T, 0, comm);
        MPI_Bcast(&z[0], n_pts * sizeof(x[0]), MPI_UINT8_T, 0, comm);
#endif
    } else {
#ifdef MSK_EXAMPLES_USE_MPI
        MPI_Bcast(&x[0], n_pts * sizeof(x[0]), MPI_UINT8_T, 0, comm);
        MPI_Bcast(&z[0], n_pts * sizeof(x[0]), MPI_UINT8_T, 0, comm);
#endif
    }

    DataType *y = malloc(n_pts * sizeof(y[0]));
    for (int i = 0; i < n_pts; ++i) {
        y[i] = 0;
    }

#ifdef MSK_EXAMPLES_USE_MPI
    msk_permute_col_inplace_z3(&msk_matrix, x, false);
    msk_matvec_z3(&msk_matrix, 'N', x + msk_matrix.local_col_off, y + msk_matrix.local_row_off);
    msk_unpermute_row_inplace_z3(&msk_matrix, y, false);
#else
    msk_permute_col_inplace_z3(&msk_matrix, x);
    msk_matvec_z3(&msk_matrix, 'N', x, y);
    msk_unpermute_row_inplace_z3(&msk_matrix, y);
#endif

    if (proc_id == 0) {
#ifdef MSK_EXAMPLES_USE_MPI
        MPI_Reduce(MPI_IN_PLACE, &y[0], n_pts * sizeof(y[0]), MPI_UINT8_T, MPI_SUM, 0, comm);
#endif
        RealType true_norm, err_norm;
        true_norm = norm(n_pts, z);
        for (uint64_t i = 0; i < n_pts; ++i) {
            y[i] -= z[i];
        }
        err_norm = norm(n_pts, y);
        printf("Matvec error = %.3e\n", err_norm / true_norm);
    } else {
#ifdef MSK_EXAMPLES_USE_MPI
        MPI_Reduce(&y[0], NULL, n_pts * sizeof(y[0]), MPI_UINT8_T, MPI_SUM, 0, comm);
#endif
    }

    msk_destroy_z3(&msk_matrix);

    free(points);
    free(dense_matrix);
    free(msk_dense);
    free(x);
    free(y);
    free(z);
    free(rhs);
    free(solution);

#ifdef MSK_EXAMPLES_USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
