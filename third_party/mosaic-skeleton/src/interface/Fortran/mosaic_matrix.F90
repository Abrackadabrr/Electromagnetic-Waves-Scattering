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
! Created by b.valiakhmetov on 06.04.26.
!


module mosaic_matrix



#include "msk_config.h"

use, intrinsic :: iso_c_binding
#ifdef MSK_USE_MPI
use :: mpi_f08
#endif

use :: msk_c_wrap
use :: msk_fort_wrap

implicit none

type MSk_MosaicMatrix
    type(MSkWrapCMosaicMatrix), pointer, public :: c_mosaic_matrix

contains
    ! sizes and offsets
    procedure, public :: row_size => fmsk_row_size
    procedure, public :: col_size => fmsk_col_size

    ! norm
    procedure, public :: frob_norm => fmsk_frob_norm

    ! memory
    procedure, public :: memory => fmsk_memory

    ! MPI related offsets and sizes
#ifdef MSK_USE_MPI
    procedure, public :: local_row_size => fmsk_local_row_size
    procedure, public :: local_col_size => fmsk_local_col_size
    procedure, public :: local_row_off => fmsk_local_row_off
    procedure, public :: local_col_off => fmsk_local_col_off

    procedure, public :: local_frob_norm => fmsk_local_frob_norm

    procedure, public :: local_memory => fmsk_local_memory
#endif

    ! helpers
    procedure, public, nopass :: set_parameters => fmsk_set_parameters


end type MSk_MosaicMatrix


contains


integer(c_int64_t) function fmsk_row_size(this) result(size)
    use, intrinsic :: iso_c_binding

    class(MSk_MosaicMatrix) :: this

    size = this%c_mosaic_matrix%row_size

end function fmsk_row_size


integer(c_int64_t) function fmsk_col_size(this) result(size)
    use, intrinsic :: iso_c_binding

    class(MSk_MosaicMatrix) :: this

    size = this%c_mosaic_matrix%col_size

end function fmsk_col_size


#ifdef MSK_USE_MPI
integer(c_int64_t) function fmsk_local_row_size(this) result(size)
    use, intrinsic :: iso_c_binding

    class(MSk_MosaicMatrix) :: this

    size = this%c_mosaic_matrix%local_row_size

end function fmsk_local_row_size


integer(c_int64_t) function fmsk_local_col_size(this) result(size)
    use, intrinsic :: iso_c_binding

    class(MSk_MosaicMatrix) :: this

    size = this%c_mosaic_matrix%local_col_size

end function fmsk_local_col_size


integer(c_int64_t) function fmsk_local_row_off(this) result(off)
    use, intrinsic :: iso_c_binding

    class(MSk_MosaicMatrix) :: this

    off = this%c_mosaic_matrix%local_row_off

end function fmsk_local_row_off


integer(c_int64_t) function fmsk_local_col_off(this) result(off)
    use, intrinsic :: iso_c_binding

    class(MSk_MosaicMatrix) :: this

    off = this%c_mosaic_matrix%local_col_off

end function fmsk_local_col_off
#endif


real(c_double) function fmsk_frob_norm(this) result(norm)
    use, intrinsic :: iso_c_binding

    class(MSk_MosaicMatrix) :: this

    norm = this%c_mosaic_matrix%frob_norm

end function fmsk_frob_norm


#ifdef MSK_USE_MPI
real(c_double) function fmsk_local_frob_norm(this) result(norm)
    use, intrinsic :: iso_c_binding

    class(MSk_MosaicMatrix) :: this

    norm = this%c_mosaic_matrix%local_frob_norm

end function fmsk_local_frob_norm
#endif


integer(c_int64_t) function fmsk_memory(this) result(mem)
    use, intrinsic :: iso_c_binding

    class(MSk_MosaicMatrix) :: this

    mem = this%c_mosaic_matrix%memory

end function fmsk_memory


#ifdef MSK_USE_MPI
integer(c_int64_t) function fmsk_local_memory(this) result(mem)
    use, intrinsic :: iso_c_binding

    class(MSk_MosaicMatrix) :: this

    mem = this%c_mosaic_matrix%local_memory

end function fmsk_local_memory
#endif



subroutine fmsk_set_parameters(&
        rel_tol, abs_tol, max_rank, &
        min_rows, min_cols, min_block_size, rho, &
        log_file, log_level, &
        c_params)
    use, intrinsic :: iso_c_binding
    
    real(c_double) :: rel_tol
    real(c_double), optional :: abs_tol
    integer, optional :: max_rank

    integer, optional :: min_rows, min_cols, min_block_size
    real(c_double), optional :: rho

    character(kind=c_char, len=*), target, optional :: log_file, log_level

    type(MSkWrapCParameters), intent(out) :: c_params

    ! Local variables

    character(kind=c_char, len=256), target :: c_log_file = "stdout"//c_null_char
    character(kind=c_char, len=256), target :: c_log_level = "off"//c_null_char


    ! Fill struct of parameters

    c_params%rel_tol = rel_tol

    ! Defaults
    c_params%abs_tol = 0
    c_params%max_rank = 0
    c_params%min_rows = 64
    c_params%min_cols = 64
    c_params%rho = 1.01
    c_params%log_file = C_LOC(c_log_file(1:1))
    c_params%log_level = C_LOC(c_log_level(1:1))

    ! Use present parameters
    if (present(abs_tol)) c_params%abs_tol = abs_tol
    if (present(max_rank)) c_params%max_rank = max_rank

    if (present(min_block_size)) c_params%min_rows = min_block_size
    if (present(min_block_size)) c_params%min_cols = min_block_size
    if (present(min_rows)) c_params%min_rows = min_rows
    if (present(min_cols)) c_params%min_cols = min_cols
    if (present(rho)) c_params%rho =rho

    if (present(log_file)) then
        c_log_file = trim(log_file(1:255))//c_null_char
        c_params%log_file = C_LOC(c_log_file(1:1))
    end if
    if (present(log_level)) then
        c_log_level = trim(log_level(1:255))//c_null_char
        c_params%log_level = C_LOC(c_log_level(1:1))
    end if
    
    ! Hard sets
    c_params%sphere_norm = 0
    c_params%force_lowrank = .false.

end subroutine fmsk_set_parameters



end module mosaic_matrix
