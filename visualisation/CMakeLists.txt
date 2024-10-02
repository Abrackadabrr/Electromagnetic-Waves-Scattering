add_library(visual VTKFunctions.hpp VTKFunctions.cpp)

set(VTK_DIR "/opt/lib/vtk/lib/cmake/vtk-9.3")

find_package(VTK REQUIRED)

target_link_libraries(visual ${PROJECT_NAME} ${VTK_LIBRARIES})
target_include_directories(visual PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
