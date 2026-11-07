# Find toml11 library

find_package(toml11 CONFIG)
if(NOT toml11_FOUND)
    include(FetchContent)
    set(toml11_INSTALL ${TOML11_INSTALL_FROM_SOURCE} CACHE BOOL "Build and install toml11 lib from sources")
    FetchContent_Declare(
            toml11
            GIT_REPOSITORY https://github.com/ToruNiina/toml11.git
            GIT_TAG        v3.8.1
            PREFIX ${PROJECT_SOURCE_DIR}/external/toml11
            SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/toml11
            BINARY_DIR ${PROJECT_BINARY_DIR}/external/toml11
            TMP_DIR ${PROJECT_BINARY_DIR}/external/toml11/tmp
            STAMP_DIR ${PROJECT_BINARY_DIR}/external/toml11/stamp
            INSTALL_DIR ${PROJECT_SOURCE_DIR}/external/toml11
            LOG_DIR ${PROJECT_BINARY_DIR}/external/toml11/log
            CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
            #EXCLUDE_FROM_ALL # since cmake 3.38
    )
    FetchContent_MakeAvailable(toml11)
    set(toml11_FOUND toml11_POPULATED)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(toml11 DEFAULT_MSG toml11_FOUND)
