project(lattice)
cmake_minimum_required(VERSION 3.22)
message("RUPOR STUDY BUILDING\n")
set(CMAKE_CXX_STANDARD 23)

add_subdirectory(../../src src)
add_subdirectory(../../third_party third_party)
add_subdirectory(../../visualisation visualisation)

find_package(OpenMP)
if (OPENMP_FOUND)
    message("OpenMP Found")
endif ()

find_package(Eigen3 REQUIRED NO_MODULE)
if (EIGEN3_FOUND)
    message("Eigen Found")
endif ()

set(TARGET_LIBS eigen visual OpenMP::OpenMP_CXX)  #  + OpenMP::OpenMP_CXX

list(APPEND TASKS tests/matrix_comparison.cpp
        hexogonal/first_attempt.cpp
        hexogonal/hax_coarse_method.cpp
        shifted_structure.cpp
        lattice_formation_simple.cpp
        lattice_formation_real.cpp
        lattice_matvec_speed_up.cpp
        lattice_toeplitz_matrix_equivalence.cpp
        lattice_formation_scientific.cpp
        lattice_formation_analysis_of_blocks.cpp
        lattice_formation_compressed_matrix.cpp
        lattice_formation.cpp)


foreach (TASK IN LISTS TASKS)
    get_filename_component(FILE_NAME ${TASK} NAME_WE)
    set(EXE_NAME ${PROJECT_NAME}_${FILE_NAME})
    add_executable(${EXE_NAME} ${TASK})
    message("TASK: ${EXE_NAME}")
    target_link_libraries(${EXE_NAME} ${PROJECT_NAME} ${TARGET_LIBS})
    target_compile_options(${EXE_NAME} PRIVATE -Wconversion)
    target_include_directories(${EXE_NAME} PUBLIC "../../src")
endforeach ()
