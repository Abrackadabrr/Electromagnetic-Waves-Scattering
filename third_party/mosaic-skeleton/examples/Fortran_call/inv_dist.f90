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


module inv_dist

use, intrinsic :: iso_c_binding

implicit none


type sInvDist
    integer(c_int64_t) :: N
    real(4), dimension(:, :), pointer :: points_x, points_y

contains
    procedure, public :: init => s_inv_dist_init
    procedure, public :: element => s_inv_dist_elem

end type sInvDist


type dInvDist
    integer(c_int64_t) :: N
    real(8), dimension(:, :), pointer :: points_x, points_y

contains
    procedure, public :: init => d_inv_dist_init
    procedure, public :: element => d_inv_dist_elem

end type dInvDist

contains

! Single procision routines
subroutine s_inv_dist_init(this, N, pts_x, pts_y)
    class(sInvDist) :: this

    integer(c_int64_t), value :: N
    real(4), dimension(3, N), target :: pts_x, pts_y

    this%N = N
    this%points_x => pts_x
    this%points_y => pts_y
end subroutine s_inv_dist_init

real(4) function s_inv_dist_elem(this, i, j) result(elem)
    class(sInvDist) :: this
    integer(c_int64_t) :: i, j

    real(4), dimension(3) :: diff
    real(4) :: dist
    integer :: idx

    if (i == j) then
        elem = 0
    else
        do idx = 1, 3
            diff(idx) = this%points_x(idx, i) - this%points_y(idx, j)
        end do
        dist = 0
        do idx = 1, 3
            dist = dist + diff(idx) * diff(idx)
        end do
        dist = sqrt(dist)

        elem = 1.0 / dist
    end if
end function s_inv_dist_elem


real(4) function s_inv_dist_wrap(i, j, data) result(elem)
    integer(c_int64_t) :: i, j
    class(*) :: data

    select type(data)
    type is (sInvDist)
        elem = data%element(i, j)
    class default
        print *, 'Wrong matrix passed'
    end select
end function s_inv_dist_wrap


! Double procision routines
subroutine d_inv_dist_init(this, N, pts_x, pts_y)
    class(dInvDist) :: this

    integer(c_int64_t), value :: N
    real(8), dimension(3, N), target :: pts_x, pts_y

    this%N = N
    this%points_x => pts_x
    this%points_y => pts_y
end subroutine d_inv_dist_init

real(8) function d_inv_dist_elem(this, i, j) result(elem)
    class(dInvDist) :: this
    integer(c_int64_t) :: i, j

    real(8), dimension(3) :: diff
    real(8) :: dist
    integer :: idx

    if (i == j) then
        elem = 0
    else
        do idx = 1, 3
            diff(idx) = this%points_x(idx, i) - this%points_y(idx, j)
        end do
        dist = 0
        do idx = 1, 3
            dist = dist + diff(idx) * diff(idx)
        end do
        dist = sqrt(dist)

        elem = 1.0 / dist
    end if
end function d_inv_dist_elem

real(8) function d_inv_dist_wrap(i, j, data) result(elem)
    integer(c_int64_t) :: i, j
    class(*) :: data

    select type(data)
    type is (dInvDist)
        elem = data%element(i, j)
    class default
        print *, 'Wrong matrix passed'
    end select
end function d_inv_dist_wrap


end module inv_dist