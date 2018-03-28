add_library(WaveGrid WaveGrid.cpp Enviroment.cpp Grid.cpp ProfileBuffer.cpp Spectrum.cpp ${Environment_RCS})

add_subdirectory(tests)
add_subdirectory(visualization)

