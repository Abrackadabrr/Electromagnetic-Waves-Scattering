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
! Created by b.valiakhmetov on 9/23/23.
!


program main


#include "msk_examples_config.h"

#ifdef MSK_EXAMPLES_USE_MPI
use :: mpi_f08
#endif


use :: msk
use :: inv_dist


implicit none


real(4), external :: snrm2
external :: saxpy, sscal, sgemv, scopy


real(4) :: s_m_one = -1, s_one = 1, s_zero = 0
integer(c_int64_t) :: N = 1000
real(c_float), dimension(:, :), pointer :: points

real(c_float), dimension(:), pointer :: x, y, z, x_full, rhs, solution, solution_full
real(c_float), dimension(:, :), pointer :: msk_dense, full_dense

integer, dimension(2) :: row_blocks, col_blocks

type(sInvDist), target :: matrix
type(MSksMosaicMatrix), target :: msk_matrix
type(MSksElementMatrix), target :: func_matrix

integer(c_int64_t) :: idx, i, j
integer :: req_thr, prov_thr, ierr, proc_id = 0, n_procs = 1, size1, size2
real(4) :: error, norm

integer :: info = 0
character(len=256) :: info_message

integer(c_int64_t) :: row_size, col_size
integer(c_int64_t) :: memory
real(8) :: msk_norm, approx_error
integer(c_int64_t) :: local_row_size, local_col_size, local_row_off, local_col_off
#ifdef MSK_EXAMPLES_USE_MPI
integer(c_int64_t) :: local_memory
real(8) :: local_msk_norm
#endif


type(c_ptr) :: ptr

#ifdef MSK_EXAMPLES_USE_MPI
req_thr = MPI_THREAD_MULTIPLE
call MPI_Init_thread(req_thr, prov_thr, ierr)
if (req_thr /= prov_thr) then
    print *, 'Weak threading'
    call MPI_Finalize()
    error stop 1
end if

call MPI_Comm_rank(MPI_COMM_WORLD, proc_id, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, n_procs, ierr)
#endif

allocate(points(3, N))

do idx = 1, N
    points(1, idx) = real(idx) / real(N)
    points(2, idx) = 0
    points(3, idx) = 0
end do


call matrix%init(N, points, points)

func_matrix%M = N
func_matrix%N = N
func_matrix%row_points => points
func_matrix%col_points => points
func_matrix%matrix_data => matrix
func_matrix%matrix_element => s_inv_dist_wrap
!
! Approximate & check
!
! Perform approximation
row_blocks(1) = N / 3
row_blocks(2) = N - N / 3
col_blocks(1) = N / 3
col_blocks(2) = N - N / 3
info = msk_matrix%init_block(func_matrix, row_blocks, col_blocks, 1e-6_8)
!info = msk_matrix%init(func_matrix, 1e-6_8)
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at matrix initialization"
    error stop trim(info_message)
end if

! Check error
#ifdef MSK_EXAMPLES_USE_MPI
info = msk_matrix%check_error(func_matrix, .true., .true., approx_error)
#else
info = msk_matrix%check_error(func_matrix, .true., approx_error)
#endif
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at approx error calculation"
    error stop trim(info_message)
end if
if (proc_id == 0) then
    print "(A,E10.3,A,E8.2)", "First approx err = ", approx_error, ", required = ", 1e-6_8
end if! Recompress
info = msk_matrix%recompress(1e-4_8)
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at matrix recompression"
    error stop trim(info_message)
end if
! Check error again
#ifdef MSK_EXAMPLES_USE_MPI
info = msk_matrix%check_error(func_matrix, .true., .true., approx_error)
#else
info = msk_matrix%check_error(func_matrix, .true., approx_error)
#endif
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at approx error calculation"
    error stop trim(info_message)
end if
if (proc_id == 0) then
    print "(A,E10.3,A,E8.2)", "Recomp approx err = ", approx_error, ", required = ", 1e-4_8
end if
#ifdef MSK_EXAMPLES_USE_MPI
call MPI_Barrier(MPI_COMM_WORLD)
#endif
! Get and print stats
row_size = msk_matrix%row_size()
col_size = msk_matrix%col_size()
memory = msk_matrix%memory()
msk_norm = msk_matrix%frob_norm()
#ifdef MSK_EXAMPLES_USE_MPI
local_row_size = msk_matrix%local_row_size()
local_col_size = msk_matrix%local_col_size()
local_row_off = msk_matrix%local_row_off()
local_col_off = msk_matrix%local_col_off()
local_memory = msk_matrix%local_memory()
local_msk_norm = msk_matrix%local_frob_norm()
print "(I2,A,I6,I6,A,I6,I6,A,I8)", proc_id, " :: Local: size = ", local_row_size, local_col_size, &
        ", off = ", local_row_off, local_col_off, ", memory = ", local_memory
call MPI_Barrier(MPI_COMM_WORLD)
#else
local_row_size = row_size
local_col_size = col_size
local_row_off = 0
local_col_off = 0
#endif
if (proc_id == 0) then
    print "(A,I6,I6)", "Size = ", row_size, col_size
    print "(A,F12.2,A)", "Compression = ", 100 * real(memory, kind=8) / real(N * N, kind=8), "%"
    print "(A,F12.2)", "Norm (approx) = ", msk_norm
    print "(A,I8)", "Memory = ", memory
end if
!
! Check true error via dense
!
! Get true dense matrix
allocate(full_dense(N, N))
do j = 1, N
    do i = 1, N
        full_dense(i, j) = matrix%element(i, j)
    end do
end do
! Get dense from msk
allocate(msk_dense(N, N))
#ifdef MSK_EXAMPLES_USE_MPI
info = msk_matrix%to_dense(msk_dense, N, 0)
#else
info = msk_matrix%to_dense(msk_dense, N)
#endif
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at matrix densification"
    error stop trim(info_message)
end if
if (proc_id == 0) then
    size1 = N * N
    call saxpy(size1, s_m_one, full_dense, 1, msk_dense, 1)
    norm = snrm2(size1, full_dense, 1)
    error = snrm2(size1, msk_dense, 1)
    print "(A,E10.3)", "Dense approx err = ", error / norm
end if
!
! Check matvec
!
allocate(x(local_col_size), y(local_row_size), z(local_row_size))
! Normalized random x
call random_number(x)
x = 2 * x - 1
size1 = local_col_size
norm = snrm2(size1, x, 1)
#ifdef MSK_EXAMPLES_USE_MPI
norm = norm * norm
call MPI_Allreduce(MPI_IN_PLACE, norm, 1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
norm = sqrt(norm)
#endif
size1 = local_col_size
call sscal(size1, real(1) / norm, x, 1)
! True matvec with true dense matrix
allocate(x_full(N))
x_full = 0
size1 = local_col_size
call scopy(size1, x, 1, x_full(local_col_off + 1 : local_col_off + local_col_size), 1)
#ifdef MSK_EXAMPLES_USE_MPI
size1 = N
call MPI_Allreduce(MPI_IN_PLACE, x_full, size1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
size1 = local_row_size
size2 = N
call sgemv('N', size1, size2, s_one, full_dense(1 + local_row_off, 1), size2, x_full, 1, s_zero, z, 1)
! Pre and postpermutation and msk-matvec
#ifdef MSK_EXAMPLES_USE_MPI
info = msk_matrix%permute_col_inplace(x, .true.)
#else
info = msk_matrix%permute_col_inplace(x)
#endif
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at vector permutation"
    error stop trim(info_message)
end if

info = msk_matrix%matvec('N', x, y)
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at matvec"
    error stop trim(info_message)
end if
#ifdef MSK_EXAMPLES_USE_MPI
info = msk_matrix%unpermute_row_inplace(y, .true.)
#else
info = msk_matrix%unpermute_row_inplace(y)
#endif
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at vector backward permutation"
    error stop trim(info_message)
end if

size1 = local_row_size
call saxpy(size1, s_m_one, z, 1, y, 1)
error = snrm2(size1, y, 1)
norm = snrm2(size1, z, 1)
#ifdef MSK_EXAMPLES_USE_MPI
print "(I2,A,E10.3,A,E10.3)", proc_id, " :: Local matvec abs_err = ", error, ", norm = ", norm
norm = norm * norm
error = error * error
call MPI_Allreduce(MPI_IN_PLACE, norm, 1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_Allreduce(MPI_IN_PLACE, error, 1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
norm = sqrt(norm)
error = sqrt(error)
#endif
if (proc_id == 0) then
    print "(A,E10.3)", "Relative matvec err = ", error / norm
end if

! Solve linear system: Ax = z
allocate(solution(local_col_size), rhs(local_row_size), solution_full(N))
size1 = local_row_size
call scopy(size1, z, 1, rhs, 1)
solution = 0
info = msk_s3_approx_block_solve(func_matrix, row_blocks, col_blocks, 1e-4_8, 1, rhs, solution, 1e-3_8, min_block_size=64, &
        log_level="info")
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at block approximation and solving"
    error stop trim(info_message)
end if

! Get unpermuted x back
#ifdef MSK_EXAMPLES_USE_MPI
info = msk_matrix%unpermute_col_inplace(x, .true.)
#else
info = msk_matrix%unpermute_col_inplace(x)
#endif
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at vector backward permutation"
    error stop trim(info_message)
end if

! Calculate true residual
size1 = local_col_size
solution_full = 0
call scopy(size1, solution, 1, solution_full(local_col_off + 1 : local_col_off + local_col_size), 1)
#ifdef MSK_EXAMPLES_USE_MPI
size1 = N
call MPI_Allreduce(MPI_IN_PLACE, solution_full, size1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
size1 = local_row_size
size2 = N
call sgemv('N', size1, size2, s_one, full_dense(1 + local_row_off, 1), size2, solution_full, 1, s_zero, rhs, 1)

size1 = local_row_size
call saxpy(size1, s_m_one, z, 1, rhs, 1)
error = snrm2(size1, rhs, 1)
norm = snrm2(size1, z, 1)
#ifdef MSK_EXAMPLES_USE_MPI
print "(I2,A,E10.3,A,E10.3)", proc_id, " :: Local matvec abs_err = ", error, ", norm = ", norm
norm = norm * norm
error = error * error
call MPI_Allreduce(MPI_IN_PLACE, norm, 1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_Allreduce(MPI_IN_PLACE, error, 1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
norm = sqrt(norm)
error = sqrt(error)
#endif
if (proc_id == 0) then
    print "(A,E10.3)", "Relative residual = ", error / norm
end if

! Solution error
size1 = local_col_size
call saxpy(size1, s_m_one, x, 1, solution, 1)
error = snrm2(size1, solution, 1)
norm = snrm2(size1, x, 1)
#ifdef MSK_EXAMPLES_USE_MPI
print "(I2,A,E10.3,A,E10.3)", proc_id, " :: Local solution abs_err = ", error, ", norm = ", norm
norm = norm * norm
error = error * error
call MPI_Allreduce(MPI_IN_PLACE, norm, 1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_Allreduce(MPI_IN_PLACE, error, 1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
norm = sqrt(norm)
error = sqrt(error)
#endif
if (proc_id == 0) then
    print "(A,E10.3)", "Relative solution err = ", error / norm
end if

! Destruction
info = msk_matrix%destroy()
if (info /= 0) then
    write (info_message, "(I1, A)") info, ": Failed at matrix destruction"
    error stop trim(info_message)
end if

deallocate(points)
deallocate(msk_dense, full_dense)
deallocate(x, y, z)
deallocate(x_full)
deallocate(solution, rhs, solution_full)
#ifdef MSK_EXAMPLES_USE_MPI
call MPI_Finalize()
#endif


end program main
