# Cmake for a calculation library compilation

message("Building a library: " ${PROJECT_NAME})

file(GLOB_RECURSE SOURCES *.hpp *.cpp)

add_library(${PROJECT_NAME} STATIC ${SOURCES})

set(TARGET_LIBS eigen)

find_package(OpenMP)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

target_link_libraries(${PROJECT_NAME} PUBLIC ${TARGET_LIBS})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/..")
