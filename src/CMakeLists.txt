# Cmake for a calculation library compilation

message("Building a library: " ${EMW_LIBRARY_NAME})

file(GLOB_RECURSE SOURCES *.cpp)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto") # -march=native -ftemplate-depth=10000")

add_library(${EMW_LIBRARY_NAME} STATIC ${SOURCES})

# Find packages to use in core lib build
find_package(OpenMP)
if (OPENMP_FOUND)
    message("OpenMP Found")
else()
    message(FATAL_ERROR "src/CMakeLists.txt: You use EMW core lib, but CMake cannot find OpenMP library")
endif ()

find_package(LAPACK REQUIRED)
if (LAPACK_FOUND)
    message("Lapack found")
else()
    message(FATAL_ERROR "src/CMakeLists.txt: You use EMW core lib, but CMake cannot find LAPACK library")
endif ()

find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    message("Blas found")
else()
    message(FATAL_ERROR "src/CMakeLists.txt: You use EMW core lib, but CMake cannot find BLAS library")
endif ()

target_link_libraries(${EMW_LIBRARY_NAME} PUBLIC Eigen3::Eigen LAPACK::LAPACK -llapacke OpenMP::OpenMP_CXX ${BLAS_LIBRARIES})
target_compile_definitions(${EMW_LIBRARY_NAME} PUBLIC EIGEN_USE_BLAS)
target_include_directories(${EMW_LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include/"
        "${CMAKE_CURRENT_SOURCE_DIR}/../include/EMW/" "${CMAKE_CURRENT_SOURCE_DIR}/.."
        "${CMAKE_CURRENT_SOURCE_DIR}/../visualization/include/")

# target_compile_options(${PROJECT_NAME} PUBLIC)

# ENABLING PGO
# target_compile_options(${PROJECT_NAME} PUBLIC -fprofile-generate)
# target_link_options(${PROJECT_NAME} PUBLIC -fprofile-generate)
