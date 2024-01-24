project(electromagnetic_waves_scattering)
cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 17)
add_subdirectory(src)
add_subdirectory(third_party)

if (WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

if (WITH_EXAMPLES)
    add_subdirectory(examples)
endif ()


