find_package(Magnum REQUIRED
  MeshTools
  Primitives
  Shaders
  SceneGraph
  Sdl2Application)

add_executable(sandbox sandbox.cpp DrawingPrimitives.cpp)

target_link_libraries(sandbox
  WaveGrid
  Magnum::MeshTools
  Magnum::Primitives
  Magnum::Application
  Magnum::Magnum
  Magnum::Shaders
  Magnum::SceneGraph
  MagnumImGui)
