project(cylinder_parallel)
cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(../../src src)
add_subdirectory(../../third_party trird_party)
add_subdirectory(../../visualisation visualisation)

find_package(OpenMP)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(EXE_NAME cylinder_esa_parallel)
add_executable(${EXE_NAME} task_cylinder_parallel.cpp)
target_link_libraries(${EXE_NAME} ${PROJECT_NAME} visual)
target_include_directories(${EXE_NAME} PUBLIC "../../src")
#target_compile_options(${EXE_NAME} PRIVATE -Werror)
#target_compile_options(${EXE_NAME} PUBLIC -fopenmp)
