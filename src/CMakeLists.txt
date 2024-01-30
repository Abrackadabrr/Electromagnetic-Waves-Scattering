# Cmake for a calculation library compilation

message("Building a library: " ${PROJECT_NAME})

file(GLOB_RECURSE SOURCES *.hpp *.cpp)

add_library(${PROJECT_NAME} STATIC ${SOURCES})

set(TARGET_LIBS eigen)

target_link_libraries(${PROJECT_NAME} PUBLIC ${TARGET_LIBS})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/..")
