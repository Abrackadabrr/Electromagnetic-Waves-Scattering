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


set(SOURCE_FILES_FORTRAN
        c_wrap.F90
        fort_wrap.f90
        element_matrix.f90
        mosaic_matrix.F90
        msk.f90
        s_mosaic_matrix.F90
        d_mosaic_matrix.F90
        c_mosaic_matrix.F90
        z_mosaic_matrix.F90
)
set(SOURCE_FILES_C
        wrap_build.c
)

set_source_files_properties(${SOURCE_FILES_FORTRAN} PROPERTIES Fortran_PREPROCESS ON)
set(SOURCE_FILES ${SOURCE_FILES_C} ${SOURCE_FILES_FORTRAN})

target_sources(zaimsk PRIVATE ${SOURCE_FILES})
