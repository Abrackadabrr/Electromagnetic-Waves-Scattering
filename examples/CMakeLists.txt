message("\nEXAMPLES BUILDING STARTS\n")
file(GLOB_RECURSE TASKS task_*.cpp)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")

foreach (TASK IN LISTS TASKS)
    get_filename_component(FILE_NAME ${TASK} NAME_WE)
    set(EXE_NAME run_${FILE_NAME})
    add_executable(${EXE_NAME} ${TASK}
            ../research/lattice/lattice_formation.cpp
            ../research/lattice/Equations.hpp
            ../research/lattice/FieldCalculation.hpp)
    target_link_libraries(${EXE_NAME} ${PROJECT_NAME} visual eigen)
    target_include_directories(${EXE_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../src")
    target_compile_options(${EXE_NAME} PRIVATE -Wconversion -fsanitize=leak -fsanitize=undefined)
    target_link_options(${EXE_NAME} PRIVATE -fsanitize=undefined)
endforeach ()
