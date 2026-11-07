# Copyright 2023-2026 Bulat Valiakhmetov
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


find_package(BLAS REQUIRED)

set(EXAMPLE_FORTRAN_HELMGOLTZ_SOURCE
        helmgoltz.f90)
add_library(example_fortran_helmgoltz ${EXAMPLE_FORTRAN_HELMGOLTZ_SOURCE})

set(EXAMPLE_FORTRAN_INV_DIST_SOURCE
        inv_dist.f90)
add_library(example_fortran_inv_dist ${EXAMPLE_FORTRAN_INV_DIST_SOURCE})

set(EXAMPLE_FORTRAN_S_CALL_SOURCE
        main_s.f90)
set_source_files_properties(${EXAMPLE_FORTRAN_S_CALL_SOURCE} PROPERTIES Fortran_PREPROCESS ON)
add_executable(example_fortran_s_call ${EXAMPLE_FORTRAN_S_CALL_SOURCE})
target_link_libraries(example_fortran_s_call PRIVATE example_fortran_inv_dist zaimsk::zaimsk ${BLAS_LIBRARIES})

set(EXAMPLE_FORTRAN_D_CALL_SOURCE
        main_d.f90)
set_source_files_properties(${EXAMPLE_FORTRAN_D_CALL_SOURCE} PROPERTIES Fortran_PREPROCESS ON)
add_executable(example_fortran_d_call ${EXAMPLE_FORTRAN_D_CALL_SOURCE})
target_link_libraries(example_fortran_d_call PRIVATE example_fortran_inv_dist zaimsk::zaimsk ${BLAS_LIBRARIES})

set(EXAMPLE_FORTRAN_C_CALL_SOURCE
        main_c.f90)
set_source_files_properties(${EXAMPLE_FORTRAN_C_CALL_SOURCE} PROPERTIES Fortran_PREPROCESS ON)
add_executable(example_fortran_c_call ${EXAMPLE_FORTRAN_C_CALL_SOURCE})
target_link_libraries(example_fortran_c_call PRIVATE example_fortran_helmgoltz zaimsk::zaimsk ${BLAS_LIBRARIES})

set(EXAMPLE_FORTRAN_Z_CALL_SOURCE
        main_z.f90)
set_source_files_properties(${EXAMPLE_FORTRAN_Z_CALL_SOURCE} PROPERTIES Fortran_PREPROCESS ON)
add_executable(example_fortran_z_call ${EXAMPLE_FORTRAN_Z_CALL_SOURCE})
target_link_libraries(example_fortran_z_call PRIVATE example_fortran_helmgoltz zaimsk::zaimsk ${BLAS_LIBRARIES})
