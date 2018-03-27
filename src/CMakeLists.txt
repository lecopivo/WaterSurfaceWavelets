add_library(WaveGrid WaveGrid.cpp Enviroment.cpp Grid.cpp ProfileBuffer.cpp Spectrum.cpp)

add_subdirectory(tests)
add_subdirectory(visualization)

