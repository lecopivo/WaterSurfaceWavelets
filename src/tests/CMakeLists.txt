find_package(Magnum REQUIRED Sdl2Application Shaders)

add_executable(sandbox sandbox.cpp)

target_link_libraries(sandbox
  WaveGrid
  Magnum::Application
  Magnum::Shaders
  Magnum::Magnum)
