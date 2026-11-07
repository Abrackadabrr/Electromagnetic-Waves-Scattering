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


      subroutine cdotc_stackreturn(r, n, x, incx, y, incy)
        implicit none
        integer n, incx, incy
        complex r, x(*), y(*), cdotc
        external cdotc
        r = cdotc(n, x, incx, y, incy)
      end

      subroutine cdotu_stackreturn(r, n, x, incx, y, incy)
        implicit none
        integer n, incx, incy
        complex r, x(*), y(*), cdotu
        external cdotu
        r = cdotu(n, x, incx, y, incy)
      end

      subroutine zdotc_stackreturn(r, n, x, incx, y, incy)
        implicit none
        integer n, incx, incy
        complex(8) r, x(*), y(*), zdotc
        external zdotc
        r = zdotc(n, x, incx, y, incy)
      end

      subroutine zdotu_stackreturn(r, n, x, incx, y, incy)
        implicit none
        integer n, incx, incy
        complex(8) r, x(*), y(*), zdotu
        external zdotu
        r = zdotu(n, x, incx, y, incy)
      end
