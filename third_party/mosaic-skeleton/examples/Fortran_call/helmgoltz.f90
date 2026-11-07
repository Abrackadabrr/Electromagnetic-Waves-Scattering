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


module helmgoltz

use, intrinsic :: iso_c_binding

implicit none


type cHelmgoltz
    complex(4) :: k
    integer(c_int64_t) :: N
    real(4), dimension(:, :), pointer :: points_x, points_y

contains
    procedure, public :: init => c_helmgoltz_init
    procedure, public :: element => c_helm_elem

end type cHelmgoltz

type zHelmgoltz
    complex(8) :: k
    integer(c_int64_t) :: N
    real(8), dimension(:, :), pointer :: points_x, points_y

contains
    procedure, public :: init => z_helmgoltz_init
    procedure, public :: element => z_helm_elem

end type zHelmgoltz


contains

subroutine c_helmgoltz_init(this, k, N, pts_x, pts_y)
    class(cHelmgoltz) :: this

    integer(c_int64_t), value :: N
    real(4), dimension(3, N), target :: pts_x, pts_y
    complex(4) :: k

    this%k = k
    this%N = N
    this%points_x => pts_x
    this%points_y => pts_y
end subroutine c_helmgoltz_init

complex(4) function c_helm_elem(this, i, j) result(elem)
    class(cHelmgoltz) :: this
    integer(c_int64_t) :: i, j

    real(4), parameter :: PI = 4 * atan(1.0_8)
    real(4), dimension(3) :: diff
    real(4) :: rad
    integer :: idx

    if (i == j) then
        elem = 0
    else
        do idx = 1, 3
            diff(idx) = this%points_x(idx, i) - this%points_y(idx, j)
        end do
        rad = 0
        do idx = 1, 3
            rad = rad + diff(idx) * diff(idx)
        end do
        rad = sqrt(rad)

        elem = 1.0 / (4 * PI) * exp(cmplx(0, 1) * this%k * rad) / rad
    end if
end function c_helm_elem

complex(4) function c_helm_wrap(i, j, data) result(elem)
    integer(c_int64_t) :: i, j
    class(*) :: data

    select type(data)
    type is (cHelmgoltz)
        elem = data%element(i, j)
    class default
        print *, 'Wrong matrix passed'
    end select
end function c_helm_wrap


subroutine z_helmgoltz_init(this, k, N, pts_x, pts_y)
    class(zHelmgoltz) :: this

    integer(c_int64_t), value :: N
    real(8), dimension(3, N), target :: pts_x, pts_y
    complex(8) :: k

    this%k = k
    this%N = N
    this%points_x => pts_x
    this%points_y => pts_y
end subroutine z_helmgoltz_init

complex(8) function z_helm_elem(this, i, j) result(elem)
    class(zHelmgoltz) :: this
    integer(c_int64_t) :: i, j

    real(8), parameter :: PI = 4 * atan(1.0_8)
    real(8), dimension(3) :: diff
    real(8) :: rad
    integer :: idx

    if (i == j) then
        elem = 0
    else
        do idx = 1, 3
            diff(idx) = this%points_x(idx, i) - this%points_y(idx, j)
        end do
        rad = 0
        do idx = 1, 3
            rad = rad + diff(idx) * diff(idx)
        end do
        rad = sqrt(rad)

        elem = 1.0 / (4 * PI) * exp(cmplx(0, 1) * this%k * rad) / rad
    end if
end function z_helm_elem

complex(8) function z_helm_wrap(i, j, data) result(elem)
    integer(c_int64_t) :: i, j
    class(*) :: data

    select type(data)
    type is (zHelmgoltz)
        elem = data%element(i, j)
    class default
        print *, 'Wrong matrix passed'
    end select
end function z_helm_wrap



end module helmgoltz