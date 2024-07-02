project(electromagnetic_waves_scattering)
cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 23)
add_subdirectory(src)
add_subdirectory(third_party)
add_subdirectory(visualisation)

if (WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()


if (WITH_TASKS)
    add_subdirectory(examples)
endif ()


