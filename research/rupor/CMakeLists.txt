project(rupor)
cmake_minimum_required(VERSION 3.22)
message("RUPOR STUDY BUILDING\n")
set(CMAKE_CXX_STANDARD 23)

add_subdirectory(../../src src)
add_subdirectory(../../third_party third_party)
add_subdirectory(../../visualisation visualisation)

set(EXE_NAME rupor_study)
add_executable(${EXE_NAME} first_calculation.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    message("OpenMP Found")
endif ()

find_package (Eigen3 REQUIRED NO_MODULE)
if (EIGEN_FOUND)
    message("Eigen Found")
endif ()

target_link_libraries (${EXE_NAME} eigen)
target_link_libraries(${EXE_NAME} OpenMP::OpenMP_CXX)
target_link_libraries(${EXE_NAME} ${PROJECT_NAME} visual)
target_include_directories(${EXE_NAME} PUBLIC )
target_include_directories(${EXE_NAME} PUBLIC "../../src")