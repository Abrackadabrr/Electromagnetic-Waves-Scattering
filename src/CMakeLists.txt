# Cmake for a calculation library compilation

message("Building a library: " ${EMW_LIBRARY_NAME})

file(GLOB_RECURSE SOURCES *.cpp)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}") # -flto -march=native ftemplate-depth=10000")

add_library(${EMW_LIBRARY_NAME} STATIC ${SOURCES})

target_link_libraries(${EMW_LIBRARY_NAME} PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX BLAS::BLAS) # LAPACK::LAPACK -llapacke
target_compile_definitions(${EMW_LIBRARY_NAME} PUBLIC EIGEN_USE_BLAS)
target_include_directories(${EMW_LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include/"
        "${CMAKE_CURRENT_SOURCE_DIR}/../include/EMW/" "${CMAKE_CURRENT_SOURCE_DIR}/..")

