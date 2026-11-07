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


module msk_fort_wrap


use,intrinsic :: iso_c_binding
use :: msk_element_matrix
implicit none


contains

real(c_float) function msk_s_fort_eval_element(i, j, data) bind(C)
    integer(c_int64_t), value :: i, j
    type(c_ptr), value :: data

    type(MSksElementMatrix), pointer :: matrix

    call C_F_POINTER(data, matrix)
    msk_s_fort_eval_element = matrix%matrix_element(i + 1, j + 1, matrix%matrix_data)
end function msk_s_fort_eval_element


real(c_double) function msk_d_fort_eval_element(i, j, data) bind(C)
    integer(c_int64_t), value :: i, j
    type(c_ptr), value :: data

    type(MSkdElementMatrix), pointer :: matrix

    call C_F_POINTER(data, matrix)
    msk_d_fort_eval_element = matrix%matrix_element(i + 1, j + 1, matrix%matrix_data)
end function msk_d_fort_eval_element


complex(c_float_complex) function msk_c_fort_eval_element(i, j, data) bind(C)
    integer(c_int64_t), value :: i, j
    type(c_ptr), value :: data

    type(MSkcElementMatrix), pointer :: matrix

    call C_F_POINTER(data, matrix)
    msk_c_fort_eval_element = matrix%matrix_element(i + 1, j + 1, matrix%matrix_data)
end function msk_c_fort_eval_element


complex(c_double_complex) function msk_z_fort_eval_element(i, j, data) bind(C)
    integer(c_int64_t), value :: i, j
    type(c_ptr), value :: data

    type(MSkzElementMatrix), pointer :: matrix

    call C_F_POINTER(data, matrix)
    msk_z_fort_eval_element = matrix%matrix_element(i + 1, j + 1, matrix%matrix_data)
end function msk_z_fort_eval_element













end module msk_fort_wrap