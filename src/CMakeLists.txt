find_package(Eigen3 CONFIG REQUIRED)

add_library(WaveGrid WaveGrid.cpp Enviroment.cpp Grid.cpp ProfileBuffer.cpp Spectrum.cpp ${Environment_RCS})
target_link_libraries(WaveGrid Eigen3::Eigen)

add_subdirectory(tests)
add_subdirectory(visualization)

