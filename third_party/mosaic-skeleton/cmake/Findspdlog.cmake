# Find spdlog logging library

find_package(spdlog CONFIG)
if(NOT spdlog_FOUND)
    include(FetchContent)
    set(SPDLOG_INSTALL ${SPDLOG_INSTALL_FROM_SOURCE} CACHE BOOL "Build and install spdlog lib from sources")
    FetchContent_Declare(
            spdlog
            GIT_REPOSITORY https://github.com/gabime/spdlog.git
            GIT_TAG        v1.x
            PREFIX ${PROJECT_SOURCE_DIR}/external/spdlog
            SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/spdlog
            BINARY_DIR ${PROJECT_BINARY_DIR}/external/spdlog
            TMP_DIR ${PROJECT_BINARY_DIR}/external/spdlog/tmp
            STAMP_DIR ${PROJECT_BINARY_DIR}/external/spdlog/stamp
            INSTALL_DIR ${PROJECT_SOURCE_DIR}/external/spdlog
            LOG_DIR ${PROJECT_BINARY_DIR}/external/spdlog/log
            CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
            #EXCLUDE_FROM_ALL # since cmake 3.38
    )
    FetchContent_MakeAvailable(spdlog)
    set(spdlog_FOUND spdlog_POPULATED)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(spdlog DEFAULT_MSG spdlog_FOUND)
