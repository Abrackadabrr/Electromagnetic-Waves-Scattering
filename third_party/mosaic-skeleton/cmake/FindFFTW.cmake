if(NOT FFTW_ROOT)
    if(DEFINED ENV{FFTWDIR})
        set(FFTW_ROOT $ENV{FFTWDIR})
    elseif(DEFINED ENV{FFTW_DIR})
        set(FFTW_ROOT $ENV{FFTW_DIR})
    endif()
endif()

find_package(FFTW CONFIG)

# Check if we can use PkgConfig
find_package(PkgConfig)

# Determine from PKG
if(PKG_CONFIG_FOUND AND NOT FFTW_ROOT)
    pkg_check_modules(PKG_FFTW QUIET "fftw3")
endif()

# Check whether to search static or dynamic libs
set(CMAKE_FIND_LIBRARY_SUFFIXES_SAVE ${CMAKE_FIND_LIBRARY_SUFFIXES})

if(${FFTW_USE_STATIC_LIBS})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAVE})
endif()

if(FFTW_ROOT)
    # Single precision
    find_library(
            FFTW_FLOAT_LIB
            NAMES "fftw3f" libfftw3f-3
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            FFTW_FLOAT_THREADS_LIB
            NAMES "fftw3f_threads"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            FFTW_FLOAT_OPENMP_LIB
            NAMES "fftw3f_omp"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            FFTW_FLOAT_MPI_LIB
            NAMES "fftw3f_mpi"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )

    # Double precision
    find_library(
            FFTW_DOUBLE_LIB
            NAMES "fftw3" libfftw3-3
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            FFTW_DOUBLE_THREADS_LIB
            NAMES "fftw3_threads"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            FFTW_DOUBLE_OPENMP_LIB
            NAMES "fftw3_omp"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            FFTW_DOUBLE_MPI_LIB
            NAMES "fftw3_mpi"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )

    # Quad precision
    find_library(
            FFTW_QUAD_LIB
            NAMES "fftw3l" libfftw3l-3
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            FFTW_QUAD_THREADS_LIB
            NAMES "fftw3l_threads"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            FFTW_QUAD_OPENMP_LIB
            NAMES "fftw3l_omp"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            FFTW_QUAD_MPI_LIB
            NAMES "fftw3l_mpi"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )

    # Includes paths
    find_path(FFTW_INCLUDE_DIRS
            NAMES "fftw3.h"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH
    )

else()

    # Single precision
    find_library(
            FFTW_FLOAT_LIB
            NAMES "fftw3f"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )
    find_library(
            FFTW_FLOAT_THREADS_LIB
            NAMES "fftw3f_threads"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )
    find_library(
            FFTW_FLOAT_OPENMP_LIB
            NAMES "fftw3f_omp"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )
    find_library(
            FFTW_FLOAT_MPI_LIB
            NAMES "fftw3f_mpi"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )

    # Double precision
    find_library(
            FFTW_DOUBLE_LIB
            NAMES "fftw3"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )
    find_library(
            FFTW_DOUBLE_THREADS_LIB
            NAMES "fftw3_threads"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )
    find_library(
            FFTW_DOUBLE_OPENMP_LIB
            NAMES "fftw3_omp"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )
    find_library(
            FFTW_DOUBLE_MPI_LIB
            NAMES "fftw3_mpi"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )

    # Quad precision
    find_library(
            FFTW_QUAD_LIB
            NAMES "fftw3l"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )
    find_library(
            FFTW_QUAD_THREADS_LIB
            NAMES "fftw3l_threads"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )
    find_library(FFTW_QUAD_OPENMP_LIB
            NAMES "fftw3l_omp"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )
    find_library(FFTW_QUAD_MPI_LIB
            NAMES "fftw3l_mpi"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )

    # Includes paths
    find_path(FFTW_INCLUDE_DIRS
            NAMES "fftw3.h"
            PATHS ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
    )

endif(FFTW_ROOT)

# Set CMake variables

# Single precision
if (FFTW_FLOAT_LIB)
    set(FFTW_Float_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_FLOAT_LIB})
    add_library(FFTW::Float INTERFACE IMPORTED)
    set_target_properties(FFTW::Float
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_FLOAT_LIB}"
    )
else()
    set(FFTW_Float_FOUND FALSE)
endif()
if (FFTW_FLOAT_THREADS_LIB)
    set(FFTW_Float_Threads_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_FLOAT_THREADS_LIB})
    add_library(FFTW::Float::Threads INTERFACE IMPORTED)
    set_target_properties(FFTW::Float::Threads
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_FLOAT_THREADS_LIB}"
    )
else()
    set(FFTW_Float_Threads_FOUND FALSE)
endif()
if (FFTW_FLOAT_OPENMP_LIB)
    set(FFTW_Float_OpenMP_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_FLOAT_OPENMP_LIB})
    add_library(FFTW::Float::OpenMP INTERFACE IMPORTED)
    set_target_properties(FFTW::Float::OpenMP
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_FLOAT_OPENMP_LIB}"
    )
else()
    set(FFTW_Float_OpenMP_FOUND FALSE)
endif()
if (FFTW_FLOAT_MPI_LIB)
    set(FFTW_Float_MPI_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_FLOAT_MPI_LIB})
    add_library(FFTW::Float::MPI INTERFACE IMPORTED)
    set_target_properties(FFTW::Float::MPI
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_FLOAT_MPI_LIB}"
    )
else()
    set(FFTW_Float_MPI_FOUND FALSE)
endif()

# Double precision
if (FFTW_DOUBLE_LIB)
    set(FFTW_Double_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_LIB})
    add_library(FFTW::Double INTERFACE IMPORTED)
    set_target_properties(FFTW::Double
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_LIB}"
    )
else()
    set(FFTW_Double_FOUND FALSE)
endif()
if (FFTW_DOUBLE_THREADS_LIB)
    set(FFTW_Double_Threads_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_THREADS_LIB})
    add_library(FFTW::Double::Threads INTERFACE IMPORTED)
    set_target_properties(FFTW::Double::Threads
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_THREADS_LIB}"
    )
else()
    set(FFTW_DDouble_Threads_FOUND FALSE)
endif()
if (FFTW_DOUBLE_OPENMP_LIB)
    set(FFTW_Double_OpenMP_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_OPENMP_LIB})
    add_library(FFTW::Double::OpenMP INTERFACE IMPORTED)
    set_target_properties(FFTW::Double::OpenMP
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_OPENMP_LIB}"
    )
else()
    set(FFTW_Double_OpenMP_FOUND FALSE)
endif()
if (FFTW_DOUBLE_MPI_LIB)
    set(FFTW_Double_MPI_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_MPI_LIB})
    add_library(FFTW::Double::MPI INTERFACE IMPORTED)
    set_target_properties(FFTW::Double::MPI
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_MPI_LIB}"
    )
else()
    set(FFTW_Double_MPI_FOUND FALSE)
endif()

# Quad precision
if (FFTW_QUAD_LIB)
    set(FFTW_Quad_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_QUAD_LIB})
    add_library(FFTW::Quad INTERFACE IMPORTED)
    set_target_properties(FFTW::Quad
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_QUAD_LIB}"
    )
else()
    set(FFTW_Quad_FOUND FALSE)
endif()
if (FFTW_QUAD_THREADS_LIB)
    set(FFTW_Quad_Threads_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_QUAD_THREADS_LIB})
    add_library(FFTW::Quad::Threads INTERFACE IMPORTED)
    set_target_properties(FFTW::Quad::Threads
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_QUAD_THREADS_LIB}"
    )
else()
    set(FFTW_Quad_Threads_FOUND FALSE)
endif()
if (FFTW_QUAD_OPENMP_LIB)
    set(FFTW_Quad_OpenMP_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_QUAD_OPENMP_LIB})
    add_library(FFTW::Quad::OpenMP INTERFACE IMPORTED)
    set_target_properties(FFTW::Quad::OpenMP
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_QUAD_OPENMP_LIB}"
    )
else()
    set(FFTW_Quad_OpenMP_FOUND FALSE)
endif()
if (FFTW_QUAD_MPI_LIB)
    set(FFTW_Quad_MPI_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_QUAD_MPI_LIB})
    add_library(FFTW::Quad::MPI INTERFACE IMPORTED)
    set_target_properties(FFTW::Quad::MPI
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW_QUAD_MPI_LIB}"
    )
else()
    set(FFTW_Quad_MPI_FOUND FALSE)
endif()


set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAVE})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(FFTW
        REQUIRED_VARS FFTW_LIBRARIES FFTW_INCLUDE_DIRS
        HANDLE_COMPONENTS
)

mark_as_advanced(
        FFTW_INCLUDE_DIRS
        FFTW_LIBRARIES
        FFTW_FLOAT_LIB
        FFTW_DOUBLE_LIB
        FFTW_QUAD_LIB
        FFTW_FLOAT_THREADS_LIB
        FFTW_DOUBLE_THREADS_LIB
        FFTW_QUAD_THREADS_LIB
        FFTW_FLOAT_OPENMP_LIB
        FFTW_DOUBLE_OPENMP_LIB
        FFTW_QUAD_OPENMP_LIB
        FFTW_FLOAT_MPI_LIB
        FFTW_DOUBLE_MPI_LIB
        FFTW_QUAD_MPI_LIB
)