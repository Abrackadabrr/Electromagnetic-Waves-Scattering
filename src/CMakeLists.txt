# Cmake for a calculation library compilation

message("Building a library for" ${PROJECT_NAME})

file(GLOB_RECURSE SOURCES *.hpp *.cpp)

add_library(${PROJECT_NAME} INTERFACE ${SOURCES})

set(TARGET_LIBS Eigen)

target_link_libraries(${PROJECT_NAME} INTERFACE ${TARGET_LIBS})
target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/..")
