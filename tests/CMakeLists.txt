message("TESTS BUILDING STARTS\n")

file(GLOB_RECURSE TESTS test_*.cpp)

foreach (TEST IN LISTS TESTS)
    get_filename_component(FILE_NAME ${TEST} NAME_WE)
    set(EXE_NAME run_${FILE_NAME})
    add_executable(${EXE_NAME} ${TEST})
    message("TEST: ${EXE_NAME}")
    target_link_libraries(${EXE_NAME} ${PROJECT_NAME} gtest gtest_main)
endforeach ()