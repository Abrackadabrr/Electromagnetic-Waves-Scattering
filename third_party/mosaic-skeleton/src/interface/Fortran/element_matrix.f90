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


module msk_element_matrix


use,intrinsic :: iso_c_binding
implicit none


! Fortran matrix types for user to be passed to MosaicMatrix approximation
type MSksElementMatrix
    integer(c_int64_t)                             :: M, N
    class(*), pointer                              :: matrix_data => NULL()
    procedure(msk_s_eval_element), pointer, nopass :: matrix_element
    real(c_float), dimension(:, :), pointer        :: row_points
    real(c_float), dimension(:, :), pointer        :: col_points
end type MSksElementMatrix

type MSkdElementMatrix
    integer(c_int64_t)                             :: M, N
    class(*), pointer                              :: matrix_data => NULL()
    procedure(msk_d_eval_element), pointer, nopass :: matrix_element
    real(c_double), dimension(:, :), pointer       :: row_points(:, :)
    real(c_double), dimension(:, :), pointer       :: col_points(:, :)
end type MSkdElementMatrix

type MSkcElementMatrix
    integer(c_int64_t)                             :: M, N
    class(*), pointer                              :: matrix_data => NULL()
    procedure(msk_c_eval_element), pointer, nopass :: matrix_element
    real(c_float), dimension(:, :), pointer        :: row_points(:, :)
    real(c_float), dimension(:, :), pointer        :: col_points(:, :)
end type MSkcElementMatrix

type MSkzElementMatrix
    integer(c_int64_t)                             :: M, N
    class(*), pointer                              :: matrix_data => NULL()
    procedure(msk_z_eval_element), pointer, nopass :: matrix_element
    real(c_double), dimension(:, :), pointer       :: row_points
    real(c_double), dimension(:, :), pointer       :: col_points
end type MSkzElementMatrix


abstract interface
real(c_float) function msk_s_eval_element(i, j, data)
    use, intrinsic     :: iso_c_binding
    integer(c_int64_t) :: i, j
    class(*), pointer  :: data
end function msk_s_eval_element
end interface

abstract interface
real(c_double) function msk_d_eval_element(i, j, data)
    use, intrinsic     :: iso_c_binding
    integer(c_int64_t) :: i, j
    class(*), pointer  :: data
end function msk_d_eval_element
end interface

abstract interface
complex(c_float_complex) function msk_c_eval_element(i, j, data)
    use, intrinsic     :: iso_c_binding
    integer(c_int64_t) :: i, j
    class(*), pointer  :: data
end function msk_c_eval_element
end interface

abstract interface
complex(c_double_complex) function msk_z_eval_element(i, j, data)
    use, intrinsic     :: iso_c_binding
    integer(c_int64_t) :: i, j
    class(*), pointer  :: data
end function msk_z_eval_element
end interface




















end module msk_element_matrix