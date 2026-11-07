include(GNUInstallDirs)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY "${PROJECT_BINARY_DIR}/include/" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
if(MSK_FORTRAN_MODULES)
    install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fmod/" DESTINATION "${CMAKE_INSTALL_LIBDIR}/msk/fmod")
endif()

install(TARGETS zaimsk
        EXPORT "zaimskTargets"
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

set("zaimsk_CMAKE_PACKAGE_INSTALL_SUBDIR" "share/zaimsk/cmake")

install(EXPORT "zaimskTargets"
        NAMESPACE zaimsk::
        DESTINATION "${zaimsk_CMAKE_PACKAGE_INSTALL_SUBDIR}"
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/zaimskConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMinorVersion
)

list(APPEND MSK_LIB_DEPS linalg_utilities spdlog Threads OpenMP)
list(APPEND MSK_LIB_DEPS_COMPONENTS "" "" "" "CXX")
if(MSK_USE_MPI)
        list(APPEND MSK_LIB_DEPS MPI)
        if(MSK_BUILD_FORTRAN_INTERFACE)
                list(APPEND MSK_LIB_DEPS_COMPONENTS "CXX C Fortran")
        elseif(MSK_BUILD_C_INTERFACE)
                list(APPEND MSK_LIB_DEPS_COMPONENTS "CXX C")
        else()
                list(APPEND MSK_LIB_DEPS_COMPONENTS "CXX")
        endif()
endif()
if(MSK_USE_FFT)
    list(APPEND MSK_LIB_DEPS FFTW)
    list(APPEND MSK_LIB_DEPS_COMPONENTS "Float Double")
endif()


configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/zaimskConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/zaimskConfig.cmake"
        INSTALL_DESTINATION "${zaimsk_CMAKE_PACKAGE_INSTALL_SUBDIR}"
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/zaimskConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/zaimskConfigVersion.cmake"
        DESTINATION "${zaimsk_CMAKE_PACKAGE_INSTALL_SUBDIR}"
)

install(FILES "${PROJECT_SOURCE_DIR}/cmake/FindFFTW.cmake"
              "${PROJECT_SOURCE_DIR}/cmake/Findspdlog.cmake"
              "${PROJECT_SOURCE_DIR}/cmake/Findtoml11.cmake"
        DESTINATION "${zaimsk_CMAKE_PACKAGE_INSTALL_SUBDIR}"
)
