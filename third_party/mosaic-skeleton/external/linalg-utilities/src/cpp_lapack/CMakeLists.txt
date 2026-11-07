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


find_package(LAPACK REQUIRED)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER(mangleinclude/LAPACKFortranMangle.h MACRO_NAMESPACE "LAPACK_")

set(SOURCE_FILES
        LAPACK.cpp
        LAPACK.h)

add_library(msk_cpp_lapack ${SOURCE_FILES})
target_compile_options(msk_cpp_lapack PUBLIC ${LAPACK_LINKER_FLAGS})
target_link_libraries(msk_cpp_lapack PUBLIC ${LAPACK_LIBRARIES})
target_include_directories(msk_cpp_lapack PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_include_directories(msk_cpp_lapack PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/mangleinclude)

target_link_libraries(linalg_utilities INTERFACE msk_cpp_lapack)
