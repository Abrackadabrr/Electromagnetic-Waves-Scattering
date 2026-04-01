# Cmake for a calculation library compilation

message("Building a library: " ${PROJECT_NAME})

file(GLOB_RECURSE SOURCES *.cpp)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto") # -march=native -ftemplate-depth=10000")

add_library(${PROJECT_NAME} STATIC ${SOURCES})

find_package(OpenMP)
if (OPENMP_FOUND)
    message("OpenMP Found")
endif()

find_package(LAPACK REQUIRED)
if (LAPACK_FOUND)
    message("Lapack found")
endif()

find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    message("Blas found")
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen LAPACK::LAPACK -llapacke OpenMP::OpenMP_CXX ${BLAS_LIBRARIES})
target_compile_definitions(${PROJECT_NAME} PUBLIC EIGEN_USE_BLAS)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/..")

# target_compile_options(${PROJECT_NAME} PUBLIC)

# ENABLING PGO
# target_compile_options(${PROJECT_NAME} PUBLIC -fprofile-generate)
# target_link_options(${PROJECT_NAME} PUBLIC -fprofile-generate)
