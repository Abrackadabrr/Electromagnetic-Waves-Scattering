message("TESTS BUILDING STARTS\n")

file(GLOB_RECURSE TESTS test_*.cpp)

foreach (TEST IN LISTS TESTS)
    get_filename_component(FILE_NAME ${TEST} NAME_WE)
    set(EXE_NAME run_${FILE_NAME})
    add_executable(${EXE_NAME} ${TEST}
            research/test_comparison_of_automatic_assembling_with_manual.cpp
            research/lattice_tests.hpp
            ../research/lattice/lattice_formation_real.cpp)
    message("TEST: ${EXE_NAME}")
    target_link_libraries(${EXE_NAME} ${PROJECT_NAME} gtest gtest_main visual)
    target_compile_options(${EXE_NAME} PRIVATE -Werror)
    # target_link_options(${EXE_NAME} PRIVATE -fsanitize=address)
    # target_compile_options(${EXE_NAME} PUBLIC -fsanitize=undefined)
    # target_link_options(${EXE_NAME} PUBLIC -fsanitize=undefined)
    add_test(${EXE_NAME} ${EXE_NAME})
endforeach ()
