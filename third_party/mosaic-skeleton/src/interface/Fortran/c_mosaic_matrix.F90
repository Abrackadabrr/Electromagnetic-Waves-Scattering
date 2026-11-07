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


module msk_c

#include "msk_config.h"

use, intrinsic :: iso_c_binding
#ifdef MSK_USE_MPI
use :: mpi_f08
#endif

use :: mosaic_matrix
use :: msk_c_wrap
use :: msk_fort_wrap

implicit none

type, extends(MSk_MosaicMatrix) :: MSkcMosaicMatrix

contains
    ! approximation
    procedure, public :: init => fmsk_c3_init
    procedure, public :: init_block => fmsk_c3_init_block
    procedure, public :: destroy => fmsk_c3_destroy
    procedure, public :: recompress => fmsk_c3_recompress
    ! matrix-vector product
    procedure, public :: matvec => fmsk_c3_matvec
    ! permutation
    procedure, public :: permute_col => fmsk_c3_permute_col
    procedure, public :: permute_col_inplace => fmsk_c3_permute_col_inplace
    procedure, public :: unpermute_col => fmsk_c3_unpermute_col
    procedure, public :: unpermute_col_inplace => fmsk_c3_unpermute_col_inplace
    procedure, public :: permute_row => fmsk_c3_permute_row
    procedure, public :: permute_row_inplace => fmsk_c3_permute_row_inplace
    procedure, public :: unpermute_row => fmsk_c3_unpermute_row
    procedure, public :: unpermute_row_inplace => fmsk_c3_unpermute_row_inplace
    ! checks
    procedure, public :: to_dense => fmsk_c3_to_dense
    procedure, public :: check_error => fmsk_c3_check_error
end type MSkcMosaicMatrix


contains

function fmsk_c3_init(this, matrix, rel_tol, &
        min_block_size, rho, n_jobs, &
        log_file, log_level, &
        abs_tol, min_rows, min_cols, max_rank) &
        result(info)
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    type(MSkcElementMatrix), target :: matrix
    real(c_double) :: rel_tol

    integer, optional :: min_block_size
    real(c_double), optional :: rho
    character(kind=c_char, len=*), target, optional :: log_file, log_level
    integer, optional :: n_jobs

    real(c_double), optional :: abs_tol
    integer, optional :: min_rows, min_cols, max_rank

    integer :: info

    ! Local variables


    type(MSkWrapCParameters), target :: c_params
#ifdef MSK_USE_MPI
    type(MPI_Comm) :: comm = MPI_COMM_WORLD
#endif

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_params_ptr
    type(c_ptr) :: c_data_ptr
    type(c_funptr) :: c_func_ptr
    type(c_ptr) :: c_row_pts_ptr, c_col_pts_ptr
    integer(c_int64_t) :: c_n_jobs = 4

    ! Set approximation parameters
    call this%set_parameters(rel_tol, abs_tol, max_rank, &
            min_rows, min_cols, min_block_size, rho, &
            log_file, log_level, c_params)
    if (present(n_jobs)) c_n_jobs = n_jobs
    
    ! Allocate matrix and set pointers
    allocate(this%c_mosaic_matrix)
    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_params_ptr = C_LOC(c_params)
    c_data_ptr = C_LOC(matrix)
    c_func_ptr = C_FUNLOC(msk_c_fort_eval_element)
    c_row_pts_ptr = C_LOC(matrix%row_points(1, 1))
    c_col_pts_ptr = C_LOC(matrix%col_points(1, 1))

#ifdef MSK_USE_MPI
    info = msk_c3_build(c_mosaic_ptr, c_params_ptr, &
                      matrix%M, matrix%N, c_func_ptr, c_data_ptr, &
                      c_row_pts_ptr, c_col_pts_ptr, c_n_jobs, comm%MPI_VAL)
#else
    info = msk_c3_build(c_mosaic_ptr, c_params_ptr, &
                      matrix%M, matrix%N, c_func_ptr, c_data_ptr, &
                      c_row_pts_ptr, c_col_pts_ptr, c_n_jobs)
#endif

end function fmsk_c3_init



function fmsk_c3_init_block(this, matrix, &
        row_block_sizes, column_block_sizes, &
        rel_tol, min_block_size, rho, n_jobs, &
        log_file, log_level, &
        abs_tol, min_rows, min_cols, max_rank) &
        result(info)
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    type(MSkcElementMatrix), target :: matrix
    integer, dimension(:) :: row_block_sizes, column_block_sizes

    real(c_double) :: rel_tol

    integer, optional :: min_block_size
    real(c_double), optional :: rho
    character(kind=c_char, len=*), target, optional :: log_file, log_level
    integer, optional :: n_jobs

    real(c_double), optional :: abs_tol
    integer, optional :: min_rows, min_cols, max_rank

    integer :: info

    ! Local variables


    type(MSkWrapCParameters), target :: c_params
#ifdef MSK_USE_MPI
    type(MPI_Comm) :: comm = MPI_COMM_WORLD
#endif

    type(c_ptr) :: c_mosaic_ptr

    integer(c_int64_t) :: c_m_blocks, c_n_blocks
    integer(c_int64_t), dimension(:), pointer :: c_m_block_sizes, c_n_block_sizes
    type(c_ptr) :: c_m_block_sizes_ptr, c_n_block_sizes_ptr

    type(c_ptr) :: c_params_ptr
    type(c_ptr) :: c_data_ptr
    type(c_funptr) :: c_func_ptr
    type(c_ptr) :: c_row_pts_ptr, c_col_pts_ptr
    integer(c_int64_t) :: c_n_jobs = 4

    integer :: i


    ! Set approximation parameters
    call this%set_parameters(rel_tol, abs_tol, max_rank, &
            min_rows, min_cols, min_block_size, rho, &
            log_file, log_level, c_params)
    if (present(n_jobs)) c_n_jobs = n_jobs

    ! Blocks
    c_m_blocks = size(row_block_sizes, 1)
    c_n_blocks = size(column_block_sizes, 1)
    allocate(c_m_block_sizes(c_m_blocks), c_n_block_sizes(c_n_blocks))
    do i = 1, c_m_blocks
        c_m_block_sizes(i) = row_block_sizes(i)
    end do
    do i = 1, c_n_blocks
        c_n_block_sizes(i) = column_block_sizes(i)
    end do
    c_m_block_sizes_ptr = C_LOC(c_m_block_sizes(1))
    c_n_block_sizes_ptr = C_LOC(c_n_block_sizes(1))

    ! Allocate matrix and set pointers
    allocate(this%c_mosaic_matrix)
    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_params_ptr = C_LOC(c_params)
    c_data_ptr = C_LOC(matrix)
    c_func_ptr = C_FUNLOC(msk_c_fort_eval_element)
    c_row_pts_ptr = C_LOC(matrix%row_points(1, 1))
    c_col_pts_ptr = C_LOC(matrix%col_points(1, 1))

#ifdef MSK_USE_MPI
    info = msk_c3_build_block(c_mosaic_ptr, c_params_ptr, &
                      matrix%M, c_m_blocks, c_m_block_sizes_ptr, &
                      matrix%N, c_n_blocks, c_n_block_sizes_ptr, &
                      c_func_ptr, c_data_ptr, &
                      c_row_pts_ptr, c_col_pts_ptr, c_n_jobs, comm%MPI_VAL)
#else
    info = msk_c3_build_block(c_mosaic_ptr, c_params_ptr, &
                      matrix%M, c_m_blocks, c_m_block_sizes_ptr, &
                      matrix%N, c_n_blocks, c_n_block_sizes_ptr, &
                      c_func_ptr, c_data_ptr, &
                      c_row_pts_ptr, c_col_pts_ptr, c_n_jobs)
#endif
    deallocate(c_m_block_sizes, c_n_block_sizes)

end function fmsk_c3_init_block


function fmsk_c3_destroy(this) result(info)
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    info = msk_c3_destroy(c_mosaic_ptr)
    deallocate(this%c_mosaic_matrix)

end function fmsk_c3_destroy


function fmsk_c3_recompress(this, rel_tol) result(info)
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    real(c_double) :: rel_tol
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    info = msk_c3_recompress(c_mosaic_ptr, rel_tol)

end function fmsk_c3_recompress


function fmsk_c3_matvec(this, trans, x, y) result(info)
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    character :: trans
    complex(c_float_complex), dimension(:), target :: x
    complex(c_float_complex), dimension(:), target :: y
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    character(c_char) :: c_trans
    type(c_ptr) :: c_x_ptr, c_y_ptr

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_trans = trans
    c_x_ptr = C_LOC(x(1))
    c_y_ptr = C_LOC(y(1))
    info = msk_c3_matvec(c_mosaic_ptr, c_trans, c_x_ptr, c_y_ptr)

end function fmsk_c3_matvec

#ifdef MSK_USE_MPI
function fmsk_c3_to_dense(this, arr, ldim, root) result(info)
#else
function fmsk_c3_to_dense(this, arr, ldim) result(info)
#endif
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    complex(c_float_complex), dimension(:, :), target :: arr
    integer(c_int64_t) :: ldim
#ifdef MSK_USE_MPI
    integer :: root
#endif
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_arr_ptr
    integer(c_int) :: c_root

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_arr_ptr = C_LOC(arr(1, 1))
#ifdef MSK_USE_MPI
    c_root = root
    info = msk_c3_to_dense(c_mosaic_ptr, c_arr_ptr, ldim, c_root)
#else
    info = msk_c3_to_dense(c_mosaic_ptr, c_arr_ptr, ldim)
#endif

end function fmsk_c3_to_dense


#ifdef MSK_USE_MPI
function fmsk_c3_permute_col(this, x_in, x_out, global) result(info)
#else
function fmsk_c3_permute_col(this, x_in, x_out) result(info)
#endif
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    complex(c_float_complex), dimension(:), target :: x_in, x_out
#ifdef MSK_USE_MPI
    logical :: global
#endif
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_x_in_ptr, c_x_out_ptr
#ifdef MSK_USE_MPI
    logical(c_bool) :: c_global
#endif

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_x_in_ptr = C_LOC(x_in(1))
    c_x_out_ptr = C_LOC(x_out(1))
#ifdef MSK_USE_MPI
    c_global = global
#endif

#ifdef MSK_USE_MPI
    info = msk_c3_permute_col(c_mosaic_ptr, c_x_in_ptr, c_x_out_ptr, c_global)
#else
    info = msk_c3_permute_col(c_mosaic_ptr, c_x_in_ptr, c_x_out_ptr)
#endif
end function fmsk_c3_permute_col


#ifdef MSK_USE_MPI
function fmsk_c3_permute_col_inplace(this, x_in_out, global) result(info)
#else
function fmsk_c3_permute_col_inplace(this, x_in_out) result(info)
#endif
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    complex(c_float_complex), dimension(:), target :: x_in_out
#ifdef MSK_USE_MPI
    logical :: global
#endif
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_x_in_out_ptr
#ifdef MSK_USE_MPI
    logical(c_bool) :: c_global
#endif

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_x_in_out_ptr = C_LOC(x_in_out(1))
#ifdef MSK_USE_MPI
    c_global = global
#endif

#ifdef MSK_USE_MPI
    info = msk_c3_permute_col_inplace(c_mosaic_ptr, c_x_in_out_ptr, c_global)
#else
    info = msk_c3_permute_col_inplace(c_mosaic_ptr, c_x_in_out_ptr)
#endif
end function fmsk_c3_permute_col_inplace


#ifdef MSK_USE_MPI
function fmsk_c3_unpermute_col(this, x_in, x_out, global) result(info)
#else
function fmsk_c3_unpermute_col(this, x_in, x_out) result(info)
#endif
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    complex(c_float_complex), dimension(:), target :: x_in, x_out
#ifdef MSK_USE_MPI
    logical :: global
#endif
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_x_in_ptr, c_x_out_ptr
#ifdef MSK_USE_MPI
    logical(c_bool) :: c_global
#endif

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_x_in_ptr = C_LOC(x_in(1))
    c_x_out_ptr = C_LOC(x_out(1))
#ifdef MSK_USE_MPI
    c_global = global
#endif

#ifdef MSK_USE_MPI
    info = msk_c3_unpermute_col(c_mosaic_ptr, c_x_in_ptr, c_x_out_ptr, c_global)
#else
    info = msk_c3_unpermute_col(c_mosaic_ptr, c_x_in_ptr, c_x_out_ptr)
#endif
end function fmsk_c3_unpermute_col


#ifdef MSK_USE_MPI
function fmsk_c3_unpermute_col_inplace(this, x_in_out, global) result(info)
#else
function fmsk_c3_unpermute_col_inplace(this, x_in_out) result(info)
#endif
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    complex(c_float_complex), dimension(:), target :: x_in_out
#ifdef MSK_USE_MPI
    logical :: global
#endif
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_x_in_out_ptr
#ifdef MSK_USE_MPI
    logical(c_bool) :: c_global
#endif

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_x_in_out_ptr = C_LOC(x_in_out(1))
#ifdef MSK_USE_MPI
    c_global = global
#endif

#ifdef MSK_USE_MPI
    info = msk_c3_unpermute_col_inplace(c_mosaic_ptr, c_x_in_out_ptr, c_global)
#else
    info = msk_c3_unpermute_col_inplace(c_mosaic_ptr, c_x_in_out_ptr)
#endif
end function fmsk_c3_unpermute_col_inplace


#ifdef MSK_USE_MPI
function fmsk_c3_permute_row(this, x_in, x_out, global) result(info)
#else
function fmsk_c3_permute_row(this, x_in, x_out) result(info)
#endif
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    complex(c_float_complex), dimension(:), target :: x_in, x_out
#ifdef MSK_USE_MPI
    logical :: global
#endif
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_x_in_ptr, c_x_out_ptr
#ifdef MSK_USE_MPI
    logical(c_bool) :: c_global
#endif

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_x_in_ptr = C_LOC(x_in(1))
    c_x_out_ptr = C_LOC(x_out(1))
#ifdef MSK_USE_MPI
    c_global = global
#endif

#ifdef MSK_USE_MPI
    info = msk_c3_permute_row(c_mosaic_ptr, c_x_in_ptr, c_x_out_ptr, c_global)
#else
    info = msk_c3_permute_row(c_mosaic_ptr, c_x_in_ptr, c_x_out_ptr)
#endif
end function fmsk_c3_permute_row


#ifdef MSK_USE_MPI
function fmsk_c3_permute_row_inplace(this, x_in_out, global) result(info)
#else
function fmsk_c3_permute_row_inplace(this, x_in_out) result(info)
#endif
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    complex(c_float_complex), dimension(:), target :: x_in_out
#ifdef MSK_USE_MPI
    logical :: global
#endif
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_x_in_out_ptr
#ifdef MSK_USE_MPI
    logical(c_bool) :: c_global
#endif

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_x_in_out_ptr = C_LOC(x_in_out(1))
#ifdef MSK_USE_MPI
    c_global = global
#endif

#ifdef MSK_USE_MPI
    info = msk_c3_permute_row_inplace(c_mosaic_ptr, c_x_in_out_ptr, c_global)
#else
    info = msk_c3_permute_row_inplace(c_mosaic_ptr, c_x_in_out_ptr)
#endif
end function fmsk_c3_permute_row_inplace


#ifdef MSK_USE_MPI
function fmsk_c3_unpermute_row(this, x_in, x_out, global) result(info)
#else
function fmsk_c3_unpermute_row(this, x_in, x_out) result(info)
#endif
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    complex(c_float_complex), dimension(:), target :: x_in, x_out
#ifdef MSK_USE_MPI
    logical :: global
#endif
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_x_in_ptr, c_x_out_ptr
#ifdef MSK_USE_MPI
    logical(c_bool) :: c_global
#endif

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_x_in_ptr = C_LOC(x_in(1))
    c_x_out_ptr = C_LOC(x_out(1))
#ifdef MSK_USE_MPI
    c_global = global
#endif

#ifdef MSK_USE_MPI
    info = msk_c3_unpermute_row(c_mosaic_ptr, c_x_in_ptr, c_x_out_ptr, c_global)
#else
    info = msk_c3_unpermute_row(c_mosaic_ptr, c_x_in_ptr, c_x_out_ptr)
#endif

end function fmsk_c3_unpermute_row


#ifdef MSK_USE_MPI
function fmsk_c3_unpermute_row_inplace(this, x_in_out, global) result(info)
#else
function fmsk_c3_unpermute_row_inplace(this, x_in_out) result(info)
#endif
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    complex(c_float_complex), dimension(:), target :: x_in_out
#ifdef MSK_USE_MPI
    logical :: global
#endif
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_x_in_out_ptr
#ifdef MSK_USE_MPI
    logical(c_bool) :: c_global
#endif

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)
    c_x_in_out_ptr = C_LOC(x_in_out(1))
#ifdef MSK_USE_MPI
    c_global = global
#endif

#ifdef MSK_USE_MPI
    info = msk_c3_unpermute_row_inplace(c_mosaic_ptr, c_x_in_out_ptr, c_global)
#else
    info = msk_c3_unpermute_row_inplace(c_mosaic_ptr, c_x_in_out_ptr)
#endif

end function fmsk_c3_unpermute_row_inplace


#ifdef MSK_USE_MPI
function fmsk_c3_check_error(this, matrix, relative, global, error) result(info)
#else
function fmsk_c3_check_error(this, matrix, relative, error) result(info)
#endif
    use, intrinsic :: iso_c_binding

    class(MSkcMosaicMatrix) :: this
    type(MSkcElementMatrix), target :: matrix

    logical :: relative
#ifdef MSK_USE_MPI
    logical :: global
#endif
    real(8) :: error
    integer :: info

    ! Local variables

    type(c_ptr) :: c_mosaic_ptr
    type(c_ptr) :: c_data_ptr
    type(c_funptr) :: c_func_ptr
    logical(c_bool) :: c_relative
#ifdef MSK_USE_MPI
    logical(c_bool) :: c_global
#endif
    real(c_double) :: c_error

    c_mosaic_ptr = C_LOC(this%c_mosaic_matrix)

    c_func_ptr = C_FUNLOC(msk_c_fort_eval_element)
    c_data_ptr = C_LOC(matrix)

    c_relative = relative
#ifdef MSK_USE_MPI
    c_global = global
#endif
#ifdef MSK_USE_MPI
    info = msk_c3_check_error(c_mosaic_ptr, c_func_ptr, c_data_ptr, c_relative, c_global, c_error)
#else
    info = msk_c3_check_error(c_mosaic_ptr, c_func_ptr, c_data_ptr, c_relative, c_error)
#endif
    error = c_error
end function fmsk_c3_check_error


function msk_c3_approx_block_solve(matrix, &
        row_block_sizes, column_block_sizes, rel_tol, &
        n_rhs, rhs, solution, sol_rel_tol, &
        min_block_size, rho, n_jobs, &
        log_file, log_level, &
        abs_tol, min_rows, min_cols, max_rank) &
        result(info)
    use, intrinsic :: iso_c_binding

    type(MSkcElementMatrix), target :: matrix
    integer, dimension(:) :: row_block_sizes, column_block_sizes
    real(c_double) :: rel_tol

    integer :: n_rhs
    complex(c_float_complex), dimension(:), pointer :: rhs, solution
    real(c_double) :: sol_rel_tol

    integer, optional :: min_block_size
    real(c_double), optional :: rho
    character(kind=c_char, len=*), target, optional :: log_file, log_level
    integer, optional :: n_jobs

    real(c_double), optional :: abs_tol
    integer, optional :: min_rows, min_cols, max_rank

    integer :: info

    ! Local variables

    type(MSkWrapCParameters), target :: c_params
#ifdef MSK_USE_MPI
    type(MPI_Comm) :: comm = MPI_COMM_WORLD
#endif

    integer(c_int64_t) :: c_m_blocks, c_n_blocks
    integer(c_int64_t), dimension(:), pointer :: c_m_block_sizes, c_n_block_sizes
    type(c_ptr) :: c_m_block_sizes_ptr, c_n_block_sizes_ptr

    integer(c_int64_t) :: c_n_rhs
    type(c_ptr) :: c_rhs_ptr, c_solution_ptr
    type(MSkWrapCSolveParameters), target :: c_sol_params
    type(c_ptr) :: c_sol_params_ptr

    type(c_ptr) :: c_params_ptr
    type(c_ptr) :: c_data_ptr
    type(c_funptr) :: c_func_ptr
    type(c_ptr) :: c_row_pts_ptr, c_col_pts_ptr
    integer(c_int64_t) :: c_n_jobs = 4

    integer :: i

    ! Set approximation parameters
    call fmsk_set_parameters(rel_tol, abs_tol, max_rank, &
            min_rows, min_cols, min_block_size, rho, &
            log_file, log_level, c_params)
    if (present(n_jobs)) c_n_jobs = n_jobs

    ! Blocks
    c_m_blocks = size(row_block_sizes, 1)
    c_n_blocks = size(column_block_sizes, 1)
    allocate(c_m_block_sizes(c_m_blocks), c_n_block_sizes(c_n_blocks))
    do i = 1, c_m_blocks
        c_m_block_sizes(i) = row_block_sizes(i)
    end do
    do i = 1, c_n_blocks
        c_n_block_sizes(i) = column_block_sizes(i)
    end do
    c_m_block_sizes_ptr = C_LOC(c_m_block_sizes(1))
    c_n_block_sizes_ptr = C_LOC(c_n_block_sizes(1))

    ! Solution params
    c_sol_params%abs_tol = 0
    c_sol_params%rel_tol = sol_rel_tol
    c_sol_params%iter_limit = 0

    c_n_rhs = n_rhs
    c_rhs_ptr = C_LOC(rhs(1))
    c_solution_ptr = C_LOC(solution(1))
    c_sol_params_ptr = C_LOC(c_sol_params)

    ! Set pointers
    c_params_ptr = C_LOC(c_params)
    c_data_ptr = C_LOC(matrix)
    c_func_ptr = C_FUNLOC(msk_c_fort_eval_element)
    c_row_pts_ptr = C_LOC(matrix%row_points(1, 1))
    c_col_pts_ptr = C_LOC(matrix%col_points(1, 1))

#ifdef MSK_USE_MPI
    info = fmsk_c3_approx_block_solve(matrix%M, c_m_blocks, c_m_block_sizes_ptr, &
            matrix%N, c_n_blocks, c_n_block_sizes_ptr, &
            c_func_ptr, c_data_ptr, &
            c_row_pts_ptr, c_col_pts_ptr, c_params_ptr, &
            c_n_rhs, c_rhs_ptr, c_solution_ptr, c_sol_params_ptr, &
            c_n_jobs, comm%MPI_VAL)
#else
    info = fmsk_c3_approx_block_solve(matrix%M, c_m_blocks, c_m_block_sizes_ptr, &
            matrix%N, c_n_blocks, c_n_block_sizes_ptr, &
            c_func_ptr, c_data_ptr, &
            c_row_pts_ptr, c_col_pts_ptr, c_params_ptr, &
            c_n_rhs, c_rhs_ptr, c_solution_ptr, c_sol_params_ptr, &
            c_n_jobs)
#endif
    deallocate(c_m_block_sizes, c_n_block_sizes)
end function msk_c3_approx_block_solve

end module msk_c