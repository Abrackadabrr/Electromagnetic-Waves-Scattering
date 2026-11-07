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
project(zaiMSk
        VERSION 1.2.0
        DESCRIPTION "Mosaic-skeleton approximations library"
        LANGUAGES CXX
)


if(NOT DEFINED MSK_MASTER_PROJECT)
        if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
                set(MSK_MASTER_PROJECT ON)
        else()
                set(MSK_MASTER_PROJECT OFF)
        endif()
endif()

option(MSK_BUILD_EXAMPLES "Build usage examples" ${MSK_MASTER_PROJECT})
option(MSK_INSTALL "Generate the MSk library installation rules" ${MSK_MASTER_PROJECT})
option(MSK_BUILD_DOCS "Build documentation" OFF)
option(MSK_BUILD_C_INTERFACE "Build C interface" OFF)
option(MSK_BUILD_FORTRAN_INTERFACE "Build Fortran interface" OFF)
option(MSK_USE_MPI "Enable MPI usage" OFF)
option(MSK_USE_FFT "Enable FFT-based convolution usage" OFF)
option(MSK_USE_SPDLOG "Enable spdlog logging library" OFF)

if(MSK_BUILD_FORTRAN_INTERFACE)
        set(MSK_BUILD_C_INTERFACE ON)
endif()

if(MSK_BUILD_C_INTERFACE)
        enable_language(C)
endif()

if(MSK_BUILD_FORTRAN_INTERFACE)
        enable_language(Fortran)
        include(FortranCInterface)
        FortranCInterface_VERIFY(CXX)
endif()

set(MSK_FORTRAN_MODULES ${MSK_BUILD_FORTRAN_INTERFACE})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_Debug "-Wall -Wextra -pedantic -fstack-protector-all -fsanitize=hwaddress -fsanitize=shadow-call-stack -fsanitize=undefined -fsanitize=address -Og -g -DDEBUG")

include(GNUInstallDirs)

add_library(zaimsk)
add_library(zaimsk::zaimsk ALIAS zaimsk)
# Keep these aliases for backward compatibility
add_library(mosaic_skeleton ALIAS zaimsk)
add_library(mosaic_skeleton::mosaic_skeleton ALIAS zaimsk)


if(MSK_FORTRAN_MODULES)
        target_include_directories(zaimsk PUBLIC
                $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fmod>
                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/msk/fmod>
        )
else()
        target_include_directories(zaimsk PUBLIC
                $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )
endif()

set_target_properties(zaimsk PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fmod)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/msk_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/msk_config.h)


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set(LINALG_BUILD_MPI ${MSK_USE_MPI})
if(MSK_USE_MPI)
        if(MSK_BUILD_FORTRAN_INTERFACE)
                find_package(MPI REQUIRED COMPONENTS CXX C Fortran)
                target_link_libraries(zaimsk PUBLIC MPI::MPI_CXX MPI::MPI_C MPI::MPI_Fortran)
        elseif(MSK_BUILD_C_INTERFACE)
                find_package(MPI REQUIRED COMPONENTS CXX C)
                target_link_libraries(zaimsk PUBLIC MPI::MPI_CXX MPI::MPI_C)
        else()
                find_package(MPI REQUIRED COMPONENTS CXX)
                target_link_libraries(zaimsk PUBLIC MPI::MPI_CXX)
        endif()
endif()

if(MSK_USE_FFT)
        find_package(FFTW REQUIRED COMPONENTS Float Double)
        target_link_libraries(zaimsk PUBLIC FFTW::Float FFTW::Double)
endif()

if (MSK_BUILD_DOCS)
        find_package(Doxygen REQUIRED)
        set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs)
        set(DOXYGEN_EXCLUDE ${PROJECT_SOURCE_DIR}/external)

        set(GENERATE_TREEVIEW YES)
        set(DOXYGEN_DISABLE_INDEX NO)
        set(DOXYGEN_FULL_SIDEBAR NO)

#        set(DOXYGEN_ENABLE_PREPROCESSING YES)
        set(DOXYGEN_ENABLED_SECTIONS MSK_USE_MPI)

        doxygen_add_docs(docs ${PROJECT_SOURCE_DIR})
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)
target_link_libraries(zaimsk PUBLIC OpenMP::OpenMP_CXX)

find_package(Threads REQUIRED)
target_link_libraries(zaimsk PUBLIC Threads::Threads)

if (MSK_USE_SPDLOG)
        set(SPDLOG_INSTALL_FROM_SOURCE ${MSK_INSTALL})
        find_package(spdlog REQUIRED)
        target_link_libraries(zaimsk PUBLIC spdlog::spdlog)
endif()


set(LINALG_UTILITIES_INSTALL ${MSK_INSTALL})
add_subdirectory(external/linalg-utilities)
target_link_libraries(zaimsk PUBLIC linalg_utilities)


add_subdirectory(src)

if(MSK_BUILD_EXAMPLES)
        add_subdirectory(examples)
endif()

if(MSK_INSTALL)
        include(cmake/install.cmake)
endif()
