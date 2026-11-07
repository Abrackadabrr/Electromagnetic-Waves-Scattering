! Copyright 2023-2026 Bulat Valiakhmetov
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
! Created by b.valiakhmetov on 9/19/23.
!


module msk_c_wrap

#include "msk_config.h"

use, intrinsic :: iso_c_binding
#ifdef MSK_USE_MPI
use :: mpi_f08
#endif
implicit none


type, bind(C) :: MSkWrapCParameters
    integer(c_int64_t) :: min_rows, min_cols
    real(c_double) :: rho
    integer(c_int64_t) :: sphere_norm
    logical(c_bool) :: force_lowrank

    real(c_double) :: abs_tol, rel_tol
    integer(c_int64_t) :: max_rank

    type(c_ptr) :: log_file
    type(c_ptr) :: log_level
end type MSkWrapCParameters


type, bind(C) :: MSkWrapCSolveParameters
    real(c_double) :: abs_tol, rel_tol;
    integer(c_int64_t) :: iter_limit
end type MSkWrapCSolveParameters


type, bind(C) :: MSkWrapCMosaicMatrix
    type(C_PTR) :: matrix

    integer(c_int64_t) :: row_size, col_size
    real(c_double) :: frob_norm
    integer(c_int64_t) :: memory

#ifdef MSK_USE_MPI
    integer(c_int64_t) :: local_row_size, local_col_size
    integer(c_int64_t) :: local_row_off, local_col_off

    real(c_double) :: local_frob_norm

    integer(c_int64_t) :: local_memory
#endif
end type MSkWrapCMosaicMatrix


! Approximation
interface
#ifdef MSK_USE_MPI
function msk_s3_build(msk_matrix, parameters, &
        m, n, f_element, data, &
        row_pts, col_pts, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_build_fort_s3") &
        result(info)
    use :: mpi_f08
#else
function msk_s3_build(msk_matrix, parameters, &
        m, n, f_element, data, &
        row_pts, col_pts, &
        n_jobs) &
        bind(C, name = "msk_build_fort_s3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: parameters
    integer(c_int64_t), value :: m, n

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: info
end function msk_s3_build
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_build(msk_matrix, parameters, &
        m, n, f_element, data, &
        row_pts, col_pts, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_build_fort_d3") &
        result(info)
        use :: mpi_f08
#else
function msk_d3_build(msk_matrix, parameters, &
        m, n, f_element, data, &
        row_pts, col_pts, &
        n_jobs) &
        bind(C, name = "msk_build_fort_d3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: parameters
    integer(c_int64_t), value :: m, n

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: info
end function msk_d3_build
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_build(msk_matrix, parameters, &
        m, n, f_element, data, &
        row_pts, col_pts, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_build_fort_c3") &
        result(info)
    use :: mpi_f08
#else
function msk_c3_build(msk_matrix, parameters, &
        m, n, f_element, data, &
        row_pts, col_pts, &
        n_jobs) &
        bind(C, name = "msk_build_fort_c3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: parameters
    integer(c_int64_t), value :: m, n

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: result
end function msk_c3_build
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_build(msk_matrix, parameters, &
        m, n, f_element, data, &
        row_pts, col_pts, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_build_fort_z3") &
        result(info)
    use :: mpi_f08
#else
function msk_z3_build(msk_matrix, parameters, &
        m, n, f_element, data, &
        row_pts, col_pts, &
        n_jobs) &
        bind(C, name = "msk_build_fort_z3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: parameters
    integer(c_int64_t), value :: m, n

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: result
end function msk_z3_build
end interface


! Blockwise approximation
interface
#ifdef MSK_USE_MPI
function msk_s3_build_block(msk_matrix, parameters, &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_build_block_fort_s3") &
        result(info)
    use :: mpi_f08
#else
function msk_s3_build_block(msk_matrix, parameters, &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        n_jobs) &
        bind(C, name = "msk_build_block_fort_s3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: parameters
    integer(c_int64_t), value :: m, m_blocks
    type(c_ptr), value :: m_block_sizes
    integer(c_int64_t), value :: n, n_blocks
    type(c_ptr), value :: n_block_sizes

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: result
end function msk_s3_build_block
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_build_block(msk_matrix, parameters, &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_build_block_fort_d3") &
        result(info)
    use :: mpi_f08
#else
function msk_d3_build_block(msk_matrix, parameters, &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        n_jobs) &
        bind(C, name = "msk_build_block_fort_d3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: parameters
    integer(c_int64_t), value :: m, m_blocks
    type(c_ptr), value :: m_block_sizes
    integer(c_int64_t), value :: n, n_blocks
    type(c_ptr), value :: n_block_sizes

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: result
end function msk_d3_build_block
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_build_block(msk_matrix, parameters, &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_build_block_fort_c3") &
        result(info)
    use :: mpi_f08
#else
function msk_c3_build_block(msk_matrix, parameters, &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        n_jobs) &
        bind(C, name = "msk_build_block_fort_c3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: parameters
    integer(c_int64_t), value :: m, m_blocks
    type(c_ptr), value :: m_block_sizes
    integer(c_int64_t), value :: n, n_blocks
    type(c_ptr), value :: n_block_sizes

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: result
end function msk_c3_build_block
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_build_block(msk_matrix, parameters, &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_build_block_fort_z3") &
        result(info)
    use :: mpi_f08
#else
function msk_z3_build_block(msk_matrix, parameters, &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        n_jobs) &
        bind(C, name = "msk_build_block_fort_z3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: parameters
    integer(c_int64_t), value :: m, m_blocks
    type(c_ptr), value :: m_block_sizes
    integer(c_int64_t), value :: n, n_blocks
    type(c_ptr), value :: n_block_sizes

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: result
end function msk_z3_build_block
end interface


! Destruction
interface
function msk_s3_destroy(msk_matrix) bind(C, name = "msk_destroy_s3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    integer(c_int) :: result
end function msk_s3_destroy
end interface

interface
function msk_d3_destroy(msk_matrix) bind(C, name = "msk_destroy_d3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    integer(c_int) :: result
end function msk_d3_destroy
end interface

interface
function msk_c3_destroy(msk_matrix) bind(C, name = "msk_destroy_c3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    integer(c_int) :: result
end function msk_c3_destroy
end interface

interface
function msk_z3_destroy(msk_matrix) bind(C, name = "msk_destroy_z3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    integer(c_int) :: result
end function msk_z3_destroy
end interface


! Recompression
interface
function msk_s3_recompress(msk_matrix, rel_tol) bind(C, name = "msk_recompress_s3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    real(c_double), value :: rel_tol
    integer(c_int) :: result
end function msk_s3_recompress
end interface

interface
function msk_d3_recompress(msk_matrix, rel_tol) bind(C, name = "msk_recompress_d3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    real(c_double), value :: rel_tol
    integer(c_int) :: result
end function msk_d3_recompress
end interface

interface
function msk_c3_recompress(msk_matrix, rel_tol) bind(C, name = "msk_recompress_c3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    real(c_double), value :: rel_tol
    integer(c_int) :: result
end function msk_c3_recompress
end interface

interface
function msk_z3_recompress(msk_matrix, rel_tol) bind(C, name = "msk_recompress_z3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    real(c_double), value :: rel_tol
    integer(c_int) :: result
end function msk_z3_recompress
end interface


! Matvec
interface
function msk_s3_matvec(msk_matrix, trans, x, y) bind(C, name = "msk_matvec_s3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    character(c_char), value :: trans
    type(c_ptr), value :: x, y
    integer(c_int) :: result
end function msk_s3_matvec
end interface

interface
function msk_d3_matvec(msk_matrix, trans, x, y) bind(C, name = "msk_matvec_d3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    character(c_char), value :: trans
    type(c_ptr), value :: x, y
    integer(c_int) :: result
end function msk_d3_matvec
end interface

interface
function msk_c3_matvec(msk_matrix, trans, x, y) bind(C, name = "msk_matvec_c3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    character(c_char), value :: trans
    type(c_ptr), value :: x, y
    integer(c_int) :: result
end function msk_c3_matvec
end interface

interface
function msk_z3_matvec(msk_matrix, trans, x, y) bind(C, name = "msk_matvec_z3") result(info)
#ifdef MSK_USE_MPI
    use :: mpi_f08
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    character(c_char), value :: trans
    type(c_ptr), value :: x, y
    integer(c_int) :: result
end function msk_z3_matvec
end interface


! Saving as dense
interface
#ifdef MSK_USE_MPI
function msk_s3_to_dense(msk_matrix, ptr, ldim, root) bind(C, name = "msk_to_dense_s3") result(info)
    use :: mpi_f08
#else
function msk_s3_to_dense(msk_matrix, ptr, ldim) bind(C, name = "msk_to_dense_s3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: ptr
    integer(c_int64_t), value :: ldim
#ifdef MSK_USE_MPI
    integer(c_int), value :: root
#endif
    integer(c_int) :: result
end function msk_s3_to_dense
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_to_dense(msk_matrix, ptr, ldim, root) bind(C, name = "msk_to_dense_d3") result(info)
    use :: mpi_f08
#else
function msk_d3_to_dense(msk_matrix, ptr, ldim) bind(C, name = "msk_to_dense_d3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: ptr
    integer(c_int64_t), value :: ldim
#ifdef MSK_USE_MPI
    integer(c_int), value :: root
#endif
    integer(c_int) :: result
end function msk_d3_to_dense
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_to_dense(msk_matrix, ptr, ldim, root) bind(C, name = "msk_to_dense_c3") result(info)
    use :: mpi_f08
#else
function msk_c3_to_dense(msk_matrix, ptr, ldim) bind(C, name = "msk_to_dense_c3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: ptr
    integer(c_int64_t), value :: ldim
#ifdef MSK_USE_MPI
    integer(c_int), value :: root
#endif
    integer(c_int) :: result
end function msk_c3_to_dense
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_to_dense(msk_matrix, ptr, ldim, root) bind(C, name = "msk_to_dense_z3") result(info)
    use :: mpi_f08
#else
function msk_z3_to_dense(msk_matrix, ptr, ldim) bind(C, name = "msk_to_dense_z3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: ptr
    integer(c_int64_t), value :: ldim
#ifdef MSK_USE_MPI
    integer(c_int), value :: root
#endif
    integer(c_int) :: result
end function msk_z3_to_dense
end interface


! Direct (by columns) permutation
interface
#ifdef MSK_USE_MPI
function msk_s3_permute_col(msk_matrix, x_in, x_out, global) bind(C, name = "msk_permute_col_s3") result(info)
    use :: mpi_f08
#else
function msk_s3_permute_col(msk_matrix, x_in, x_out) bind(C, name = "msk_permute_col_s3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_s3_permute_col
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_permute_col(msk_matrix, x_in, x_out, global) bind(C, name = "msk_permute_col_d3") result(info)
    use :: mpi_f08
#else
function msk_d3_permute_col(msk_matrix, x_in, x_out) bind(C, name = "msk_permute_col_d3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_d3_permute_col
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_permute_col(msk_matrix, x_in, x_out, global) bind(C, name = "msk_permute_col_c3") result(info)
    use :: mpi_f08
#else
function msk_c3_permute_col(msk_matrix, x_in, x_out) bind(C, name = "msk_permute_col_c3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_c3_permute_col
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_permute_col(msk_matrix, x_in, x_out, global) bind(C, name = "msk_permute_col_z3") result(info)
    use :: mpi_f08
#else
function msk_z3_permute_col(msk_matrix, x_in, x_out) bind(C, name = "msk_permute_col_z3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_z3_permute_col
end interface


interface
#ifdef MSK_USE_MPI
function msk_s3_permute_col_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_permute_col_inplace_s3") result(info)
use :: mpi_f08
#else
function msk_s3_permute_col_inplace(msk_matrix, x_in_out) bind(C, name = "msk_permute_col_inplace_s3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_s3_permute_col_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_permute_col_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_permute_col_inplace_d3") result(info)
use :: mpi_f08
#else
function msk_d3_permute_col_inplace(msk_matrix, x_in_out) bind(C, name = "msk_permute_col_inplace_d3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_d3_permute_col_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_permute_col_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_permute_col_inplace_c3") result(info)
use :: mpi_f08
#else
function msk_c3_permute_col_inplace(msk_matrix, x_in_out) bind(C, name = "msk_permute_col_inplace_c3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_c3_permute_col_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_permute_col_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_permute_col_inplace_z3") result(info)
use :: mpi_f08
#else
function msk_z3_permute_col_inplace(msk_matrix, x_in_out) bind(C, name = "msk_permute_col_inplace_z3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_z3_permute_col_inplace
end interface


! Direct (by columns) inverse permutation
interface
#ifdef MSK_USE_MPI
function msk_s3_unpermute_col(msk_matrix, x_in, x_out, global) bind(C, name = "msk_unpermute_col_s3") result(info)
    use :: mpi_f08
#else
function msk_s3_unpermute_col(msk_matrix, x_in, x_out) bind(C, name = "msk_unpermute_col_s3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_s3_unpermute_col
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_unpermute_col(msk_matrix, x_in, x_out, global) bind(C, name = "msk_unpermute_col_d3") result(info)
    use :: mpi_f08
#else
function msk_d3_unpermute_col(msk_matrix, x_in, x_out) bind(C, name = "msk_unpermute_col_d3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_d3_unpermute_col
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_unpermute_col(msk_matrix, x_in, x_out, global) bind(C, name = "msk_unpermute_col_c3") result(info)
    use :: mpi_f08
#else
function msk_c3_unpermute_col(msk_matrix, x_in, x_out) bind(C, name = "msk_unpermute_col_c3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_c3_unpermute_col
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_unpermute_col(msk_matrix, x_in, x_out, global) bind(C, name = "msk_unpermute_col_z3") result(info)
    use :: mpi_f08
#else
function msk_z3_unpermute_col(msk_matrix, x_in, x_out) bind(C, name = "msk_unpermute_col_z3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_z3_unpermute_col
end interface


interface
#ifdef MSK_USE_MPI
function msk_s3_unpermute_col_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_unpermute_col_inplace_s3") result(info)
use :: mpi_f08
#else
function msk_s3_unpermute_col_inplace(msk_matrix, x_in_out) bind(C, name = "msk_unpermute_col_inplace_s3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_s3_unpermute_col_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_unpermute_col_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_unpermute_col_inplace_d3") result(info)
use :: mpi_f08
#else
function msk_d3_unpermute_col_inplace(msk_matrix, x_in_out) bind(C, name = "msk_unpermute_col_inplace_d3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_d3_unpermute_col_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_unpermute_col_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_unpermute_col_inplace_c3") result(info)
use :: mpi_f08
#else
function msk_c3_unpermute_col_inplace(msk_matrix, x_in_out) bind(C, name = "msk_unpermute_col_inplace_c3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_c3_unpermute_col_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_unpermute_col_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_unpermute_col_inplace_z3") result(info)
use :: mpi_f08
#else
function msk_z3_unpermute_col_inplace(msk_matrix, x_in_out) bind(C, name = "msk_unpermute_col_inplace_z3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_z3_unpermute_col_inplace
end interface


! Transposed (by rows) permutation
interface
#ifdef MSK_USE_MPI
function msk_s3_permute_row(msk_matrix, x_in, x_out, global) bind(C, name = "msk_permute_row_s3") result(info)
    use :: mpi_f08
#else
function msk_s3_permute_row(msk_matrix, x_in, x_out) bind(C, name = "msk_permute_row_s3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_s3_permute_row
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_permute_row(msk_matrix, x_in, x_out, global) bind(C, name = "msk_permute_row_d3") result(info)
    use :: mpi_f08
#else
function msk_d3_permute_row(msk_matrix, x_in, x_out) bind(C, name = "msk_permute_row_d3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_d3_permute_row
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_permute_row(msk_matrix, x_in, x_out, global) bind(C, name = "msk_permute_row_c3") result(info)
    use :: mpi_f08
#else
function msk_c3_permute_row(msk_matrix, x_in, x_out) bind(C, name = "msk_permute_row_c3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_c3_permute_row
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_permute_row(msk_matrix, x_in, x_out, global) bind(C, name = "msk_permute_row_z3") result(info)
    use :: mpi_f08
#else
function msk_z3_permute_row(msk_matrix, x_in, x_out) bind(C, name = "msk_permute_row_z3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_z3_permute_row
end interface


interface
#ifdef MSK_USE_MPI
function msk_s3_permute_row_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_permute_row_inplace_s3") result(info)
use :: mpi_f08
#else
function msk_s3_permute_row_inplace(msk_matrix, x_in_out) bind(C, name = "msk_permute_row_inplace_s3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_s3_permute_row_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_permute_row_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_permute_row_inplace_d3") result(info)
use :: mpi_f08
#else
function msk_d3_permute_row_inplace(msk_matrix, x_in_out) bind(C, name = "msk_permute_row_inplace_d3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_d3_permute_row_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_permute_row_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_permute_row_inplace_c3") result(info)
use :: mpi_f08
#else
function msk_c3_permute_row_inplace(msk_matrix, x_in_out) bind(C, name = "msk_permute_row_inplace_c3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_c3_permute_row_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_permute_row_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_permute_row_inplace_z3") result(info)
use :: mpi_f08
#else
function msk_z3_permute_row_inplace(msk_matrix, x_in_out) bind(C, name = "msk_permute_row_inplace_z3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_z3_permute_row_inplace
end interface


! Transposed (by rows) inverse permutation
interface
#ifdef MSK_USE_MPI
function msk_s3_unpermute_row(msk_matrix, x_in, x_out, global) bind(C, name = "msk_unpermute_row_s3") result(info)
    use :: mpi_f08
#else
function msk_s3_unpermute_row(msk_matrix, x_in, x_out) bind(C, name = "msk_unpermute_row_s3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_s3_unpermute_row
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_unpermute_row(msk_matrix, x_in, x_out, global) bind(C, name = "msk_unpermute_row_d3") result(info)
    use :: mpi_f08
#else
function msk_d3_unpermute_row(msk_matrix, x_in, x_out) bind(C, name = "msk_unpermute_row_d3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_d3_unpermute_row
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_unpermute_row(msk_matrix, x_in, x_out, global) bind(C, name = "msk_unpermute_row_c3") result(info)
    use :: mpi_f08
#else
function msk_c3_unpermute_row(msk_matrix, x_in, x_out) bind(C, name = "msk_unpermute_row_c3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_c3_unpermute_row
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_unpermute_row(msk_matrix, x_in, x_out, global) bind(C, name = "msk_unpermute_row_z3") result(info)
    use :: mpi_f08
#else
function msk_z3_unpermute_row(msk_matrix, x_in, x_out) bind(C, name = "msk_unpermute_row_z3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in, x_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_z3_unpermute_row
end interface


interface
#ifdef MSK_USE_MPI
function msk_s3_unpermute_row_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_unpermute_row_inplace_s3") result(info)
use :: mpi_f08
#else
function msk_s3_unpermute_row_inplace(msk_matrix, x_in_out) bind(C, name = "msk_unpermute_row_inplace_s3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_s3_unpermute_row_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_unpermute_row_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_unpermute_row_inplace_d3") result(info)
use :: mpi_f08
#else
function msk_d3_unpermute_row_inplace(msk_matrix, x_in_out) bind(C, name = "msk_unpermute_row_inplace_d3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_d3_unpermute_row_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_unpermute_row_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_unpermute_row_inplace_c3") result(info)
use :: mpi_f08
#else
function msk_c3_unpermute_row_inplace(msk_matrix, x_in_out) bind(C, name = "msk_unpermute_row_inplace_c3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_c3_unpermute_row_inplace
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_unpermute_row_inplace(msk_matrix, x_in_out, global) bind(C, name = "msk_unpermute_row_inplace_z3") result(info)
use :: mpi_f08
#else
function msk_z3_unpermute_row_inplace(msk_matrix, x_in_out) bind(C, name = "msk_unpermute_row_inplace_z3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    type(c_ptr), value :: x_in_out
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    integer(c_int) :: result
end function msk_z3_unpermute_row_inplace
end interface


! Approximation error check
interface
#ifdef MSK_USE_MPI
function msk_s3_check_error(msk_matrix, f_element, data, relative, global, error) bind(C, name = "msk_check_error_s3") result(info)
    use :: mpi_f08
#else
function msk_s3_check_error(msk_matrix, f_element, data, relative, error) bind(C, name = "msk_check_error_s3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    logical(c_bool), value :: relative
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    real(c_double) :: error
    integer(c_int) :: result
end function msk_s3_check_error
end interface

interface
#ifdef MSK_USE_MPI
function msk_d3_check_error(msk_matrix, f_element, data, relative, global, error) bind(C, name = "msk_check_error_d3") result(info)
    use :: mpi_f08
#else
function msk_d3_check_error(msk_matrix, f_element, data, relative, error) bind(C, name = "msk_check_error_d3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix
    logical(c_bool), value :: relative

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    real(c_double) :: error
    integer(c_int) :: result
end function msk_d3_check_error
end interface

interface
#ifdef MSK_USE_MPI
function msk_c3_check_error(msk_matrix, f_element, data, relative, global, error) bind(C, name = "msk_check_error_c3") result(info)
    use :: mpi_f08
#else
function msk_c3_check_error(msk_matrix, f_element, data, relative, error) bind(C, name = "msk_check_error_c3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    logical(c_bool), value :: relative
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    real(c_double) :: error
    integer(c_int) :: result
end function msk_c3_check_error
end interface

interface
#ifdef MSK_USE_MPI
function msk_z3_check_error(msk_matrix, f_element, data, relative, global, error) bind(C, name = "msk_check_error_z3") result(info)
    use :: mpi_f08
#else
function msk_z3_check_error(msk_matrix, f_element, data, relative, error) bind(C, name = "msk_check_error_z3") result(info)
#endif
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: msk_matrix

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    logical(c_bool), value :: relative
#ifdef MSK_USE_MPI
    logical(c_bool), value :: global
#endif
    real(c_double) :: error
    integer(c_int) :: result
end function msk_z3_check_error
end interface


interface
#ifdef MSK_USE_MPI
function fmsk_s3_approx_block_solve( &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        parameters, &
        n_rhs, rhs, solution, sol_params, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_approx_block_solve_fort_s3") &
        result(info)
        use :: mpi_f08
#else
function fmsk_s3_approx_block_solve( &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        parameters, &
        n_rhs, rhs, solution, sol_params, &
        n_jobs) &
        bind(C, name = "msk_approx_block_solve_fort_s3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    integer(c_int64_t), value :: m, m_blocks
    type(c_ptr), value :: m_block_sizes
    integer(c_int64_t), value :: n, n_blocks
    type(c_ptr), value :: n_block_sizes

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    type(c_ptr), value :: parameters

    integer(c_int64_t), value :: n_rhs
    type(c_ptr), value :: rhs, solution

    type(c_ptr), value :: sol_params

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: result
end function fmsk_s3_approx_block_solve
end interface

interface
#ifdef MSK_USE_MPI
function fmsk_d3_approx_block_solve( &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        parameters, &
        n_rhs, rhs, solution, sol_params, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_approx_block_solve_fort_d3") &
        result(info)
        use :: mpi_f08
#else
function fmsk_d3_approx_block_solve( &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        parameters, &
        n_rhs, rhs, solution, sol_params, &
        n_jobs) &
        bind(C, name = "msk_approx_block_solve_fort_d3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    integer(c_int64_t), value :: m, m_blocks
    type(c_ptr), value :: m_block_sizes
    integer(c_int64_t), value :: n, n_blocks
    type(c_ptr), value :: n_block_sizes

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    type(c_ptr), value :: parameters

    integer(c_int64_t), value :: n_rhs
    type(c_ptr), value :: rhs, solution

    type(c_ptr), value :: sol_params

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: result
end function fmsk_d3_approx_block_solve
end interface

interface
#ifdef MSK_USE_MPI
function fmsk_c3_approx_block_solve( &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        parameters, &
        n_rhs, rhs, solution, sol_params, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_approx_block_solve_fort_c3") &
        result(info)
        use :: mpi_f08
#else
function fmsk_c3_approx_block_solve( &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        parameters, &
        n_rhs, rhs, solution, sol_params, &
        n_jobs) &
        bind(C, name = "msk_approx_block_solve_fort_c3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    integer(c_int64_t), value :: m, m_blocks
    type(c_ptr), value :: m_block_sizes
    integer(c_int64_t), value :: n, n_blocks
    type(c_ptr), value :: n_block_sizes

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    type(c_ptr), value :: parameters

    integer(c_int64_t), value :: n_rhs
    type(c_ptr), value :: rhs, solution

    type(c_ptr), value :: sol_params

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: result
end function fmsk_c3_approx_block_solve
end interface

interface
#ifdef MSK_USE_MPI
function fmsk_z3_approx_block_solve( &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        parameters, &
        n_rhs, rhs, solution, sol_params, &
        n_jobs, f_comm_id) &
        bind(C, name = "msk_approx_block_solve_fort_z3") &
        result(info)
        use :: mpi_f08
#else
function fmsk_z3_approx_block_solve( &
        m, m_blocks, m_block_sizes, &
        n, n_blocks, n_block_sizes, &
        f_element, data, &
        row_pts, col_pts, &
        parameters, &
        n_rhs, rhs, solution, sol_params, &
        n_jobs) &
        bind(C, name = "msk_approx_block_solve_fort_z3") &
        result(info)
#endif
    use, intrinsic :: iso_c_binding
    integer(c_int64_t), value :: m, m_blocks
    type(c_ptr), value :: m_block_sizes
    integer(c_int64_t), value :: n, n_blocks
    type(c_ptr), value :: n_block_sizes

    type(c_funptr), value :: f_element
    type(c_ptr), value :: data

    type(c_ptr), value :: row_pts, col_pts

    type(c_ptr), value :: parameters

    integer(c_int64_t), value :: n_rhs
    type(c_ptr), value :: rhs, solution

    type(c_ptr), value :: sol_params

    integer(c_int64_t), value :: n_jobs
#ifdef MSK_USE_MPI
    integer, value :: f_comm_id
#endif
    integer(c_int) :: result
end function fmsk_z3_approx_block_solve
end interface


end module msk_c_wrap