message("TESTS BUILDING STARTS\n")

file(GLOB_RECURSE BENCHS bench_*.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    message("OpenMP Found")
endif()

find_package(benchmark REQUIRED)
if (BENCHMARK_FOUND)
    message("benchmark Found")
endif()

find_package(BLAS REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

foreach (TEST IN LISTS BENCHS)
    get_filename_component(FILE_NAME ${TEST} NAME_WE)
    set(EXE_NAME ${FILE_NAME})
    add_executable(${EXE_NAME} ${TEST})
    message("BENCH: ${EXE_NAME}")
    target_link_libraries(${EXE_NAME} PRIVATE ${PROJECT_NAME} ${BLAS_LIBRARIES} benchmark::benchmark benchmark::benchmark_main)
#    target_compile_options(${EXE_NAME} PUBLIC -flto)
#    target_link_libraries(${EXE_NAME} OpenMP::OpenMP_CXX)
#    target_compile_options(${EXE_NAME} PRIVATE -Werror)
    # target_link_options(${EXE_NAME} PRIVATE -fsanitize=address)
#    target_compile_options(${EXE_NAME} PUBLIC -fsanitize=undefined)
#    target_link_options(${EXE_NAME} PUBLIC -fsanitize=undefined)
endforeach ()
