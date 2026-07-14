# Cmake for a calculation library compilation

message("Building a library: " ${EMW_LIBRARY_NAME})

file(GLOB_RECURSE SOURCES *.cpp)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto -march=native ") #-ftemplate-depth=10000")

add_library(${EMW_LIBRARY_NAME} STATIC ${SOURCES})

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

target_link_libraries(${EMW_LIBRARY_NAME} PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX BLAS::BLAS) # LAPACK::LAPACK -llapacke
target_compile_definitions(${EMW_LIBRARY_NAME} PUBLIC EIGEN_USE_BLAS)
target_include_directories(${EMW_LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include/"
        "${CMAKE_CURRENT_SOURCE_DIR}/../include/EMW/" "${CMAKE_CURRENT_SOURCE_DIR}/..")
