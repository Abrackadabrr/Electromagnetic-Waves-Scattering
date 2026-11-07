include(GNUInstallDirs)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(MSK_USE_MPI)
    install(TARGETS linalg_utilities msk_cpp_blas msk_cpp_lapack msk_cpp_mpi
            EXPORT linalg_utilitiesTargets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
else()
    install(TARGETS linalg_utilities msk_cpp_blas msk_cpp_lapack
            EXPORT linalg_utilitiesTargets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()


set(linalg_utilities_CMAKE_PACKAGE_INSTALL_SUBDIR "share/linalg_utilities/cmake")

install(EXPORT linalg_utilitiesTargets
        NAMESPACE linalg_utilities::
        DESTINATION ${linalg_utilities_CMAKE_PACKAGE_INSTALL_SUBDIR}
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/linalg_utilitiesConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMinorVersion
)

configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/linalg_utilitiesConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/linalg_utilitiesConfig.cmake"
        INSTALL_DESTINATION ${linalg_utilities_CMAKE_PACKAGE_INSTALL_SUBDIR}
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/linalg_utilitiesConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/linalg_utilitiesConfigVersion.cmake"
        DESTINATION ${linalg_utilities_CMAKE_PACKAGE_INSTALL_SUBDIR})
