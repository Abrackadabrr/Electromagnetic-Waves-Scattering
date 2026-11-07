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


cmake_minimum_required(VERSION 3.20)
project(mosaic_skeleton_examples
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
add_library(${PROJECT_NAME} INTERFACE)

if(NOT DEFINED MSK_MASTER_PROJECT OR NOT ${MSK_MASTER_PROJECT})
        option(MSK_EXAMPLES_C "Build C interface" OFF)
        option(MSK_EXAMPLES_FORTRAN "Build Fortran interface" OFF)
        option(MSK_EXAMPLES_USE_MPI "Enable MPI usage" OFF)
else()
        option(MSK_EXAMPLES_C "Build C interface" ${MSK_BUILD_C_INTERFACE})
        option(MSK_EXAMPLES_FORTRAN "Build Fortran interface" ${MSK_BUILD_FORTRAN_INTERFACE})
        option(MSK_EXAMPLES_USE_MPI "Enable MPI usage" ${MSK_USE_MPI})
endif()

if(MSK_EXAMPLES_C)
        enable_language(C)
endif()
if(MSK_EXAMPLES_FORTRAN)
        enable_language(Fortran)
endif()


if(NOT TARGET zaimsk)
        find_package(zaimsk REQUIRED)
endif()

if(DEFINED MSK_INSTALL)
        set(TOML11_INSTALL_FROM_SOURCE ${MSK_INSTALL})
endif()
find_package(toml11 REQUIRED)


set(MSK_EXAMPLES_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
configure_file(msk_examples_config.h.in msk_examples_config.h @ONLY)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(common)

add_subdirectory(cross)
add_subdirectory(cgnr)
add_subdirectory(svd)
add_subdirectory(magnetic)
add_subdirectory(helmgoltz)
add_subdirectory(gmres)
add_subdirectory(dense)
if (MSK_EXAMPLES_C)
        add_subdirectory(C_call)
endif()
if (MSK_EXAMPLES_FORTRAN)
        add_subdirectory(Fortran_call)
endif()
add_subdirectory(particles)
add_subdirectory(prandtl)
add_subdirectory(thread_pool)
