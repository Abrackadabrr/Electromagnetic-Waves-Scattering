project(electromagnetic_waves_scattering)
cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 23)

# Introduce library name
set(EMW_LIBRARY_NAME emw)
set(EMW_LIBRARIES PARENT_SCOPE)

list(APPEND EMW_LIBRARIES ${EMW_LIBRARY_NAME})

# Find packages to use in core lib build
find_package(OpenMP)
if (OPENMP_FOUND)
    message("OpenMP is found")
else()
    message(FATAL_ERROR "src/CMakeLists.txt: You use EMW core lib, but CMake cannot find OpenMP library")
endif ()

find_package(LAPACK REQUIRED)
if (LAPACK_FOUND)
    message("Lapack is found")
else()
    message(FATAL_ERROR "src/CMakeLists.txt: You use EMW core lib, but CMake cannot find LAPACK library")
endif ()

# Find specific BLAS lib
set(BLA_VENDOR "OpenBLAS")
find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    message("OpenBlas is found")
else()
    message(FATAL_ERROR "src/CMakeLists.txt: You use EMW core lib, but CMake cannot find BLAS library")
endif ()

add_subdirectory(third_party)
add_subdirectory(src)

# Do we need VTK visualization
if (EMW_USE_VTK_VISUALIZATION)
    set(EMW_VISUAL_LIB emw_visual_vtk)
    list(APPEND EMW_LIBRARIES ${EMW_VISUAL_LIB})
    add_subdirectory(visualisation)
endif()

if (EMW_WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

if (EMW_WITH_TASKS)
    add_subdirectory(examples)
endif ()

if (EMW_WITH_BENCHMARKS)
    add_subdirectory(benchmarks)
endif ()

# Prints EMW-related libraries that is build by cmake and that should be linked to the tests and programs
message("EMW related libraries compiled: ${EMW_LIBRARIES}")
