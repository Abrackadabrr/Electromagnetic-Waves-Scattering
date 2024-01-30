add_library(visual VTKFunctions.hpp VTKFunctions.cpp)

set(VTK_DIR "/media/evgen/SecondLinuxDisk/Libraries/vtkbase/lib/cmake/vtk-9.1")
find_package(VTK REQUIRED)

target_link_libraries(visual ${PROJECT_NAME} ${VTK_LIBRARIES})
target_include_directories(visual PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})