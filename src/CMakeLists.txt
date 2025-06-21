# Cmake for a calculation library compilation

message("Building a library: " ${PROJECT_NAME})

file(GLOB_RECURSE SOURCES *.cpp)

add_library(${PROJECT_NAME} STATIC ${SOURCES})

set(TARGET_LIBS eigen)

find_package(OpenMP)
if (OPENMP_FOUND)
    message("OpenMP Found")
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(${PROJECT_NAME} PUBLIC ${TARGET_LIBS})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/..")

target_compile_options(${PROJECT_NAME} PUBLIC)

# ENABLING PGO
# target_compile_options(${PROJECT_NAME} PUBLIC -fprofile-generate)
# target_link_options(${PROJECT_NAME} PUBLIC -fprofile-generate)
