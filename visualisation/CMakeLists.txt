add_library(${EMW_VISUAL_LIB} VTKFunctions.hpp VTKFunctions.cpp)

set(VTK_DIR "/opt/lib/vtk/lib/cmake/vtk-9.3")

find_package(VTK REQUIRED)

if (${VTK_FOUND})
    message("VTK lib found")
else()
    message(FATAL_ERROR "visualisation/CMakeLists.txt: You use visualisation in EMW, but CMake cannot find VTK library")
endif()

target_link_libraries(${EMW_VISUAL_LIB} PUBLIC ${EMW_LIBRARY_NAME} ${VTK_LIBRARIES})
target_include_directories(${EMW_VISUAL_LIB} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
