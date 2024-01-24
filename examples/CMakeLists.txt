file(GLOB_RECURSE TASKS task_*.cpp)

foreach (TASK IN LISTS TASKS)
    get_filename_component(FILE_NAME ${TASK} NAME_WE)
    set(EXE_NAME run_${FILE_NAME})
    add_executable(EXE_NAME ${TASK})
    target_link_libraries(EXE_NAME ${PROJECT_NAME})
endforeach ()