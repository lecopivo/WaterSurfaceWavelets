add_library(WaveGrid WaveGrid.cpp Enviroment.cpp)

add_executable(main main.cpp)

target_link_libraries(main WaveGrid)

add_subdirectory(tests)
add_subdirectory(visualization)

