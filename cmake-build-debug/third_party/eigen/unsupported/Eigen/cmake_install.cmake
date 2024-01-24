# Install script for directory: /media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/AdolcForward"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/AlignedVector3"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/ArpackSupport"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/AutoDiff"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/BVH"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/EulerAngles"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/FFT"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/IterativeSolvers"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/KroneckerProduct"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/LevenbergMarquardt"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/MatrixFunctions"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/MoreVectorization"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/MPRealSupport"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/NNLS"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/NonLinearOptimization"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/NumericalDiff"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/OpenGLSupport"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/Polynomials"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/Skyline"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/SparseExtra"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/SpecialFunctions"
    "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/third_party/eigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/media/evgen/SecondLinuxDisk1/4_level/Electromagnetic-Waves-Scattering/cmake-build-debug/third_party/eigen/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

