message("TESTS BUILDING STARTS\n")

file(GLOB_RECURSE TESTS test_*.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    message("OpenMP Found")
endif()

foreach (TEST IN LISTS TESTS)
    get_filename_component(FILE_NAME ${TEST} NAME_WE)
    set(EXE_NAME run_${FILE_NAME})
    add_executable(${EXE_NAME} ${TEST})
    message("TEST: ${EXE_NAME}")
    target_link_libraries(${EXE_NAME} PRIVATE ${EMW_LIBRARIES} gtest gtest_main)
    add_test(${EXE_NAME} ${EXE_NAME})
endforeach ()
