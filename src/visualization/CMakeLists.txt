add_subdirectory(drawables)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

find_package(Magnum REQUIRED
  MeshTools
  Primitives
  Shaders
  SceneGraph
  Sdl2Application)

corrade_add_resource(WaterSurfaceShaders_RCS resources.conf)

# this might be important
# set_target_properties(MagnumShaders_RCS-dependencies PROPERTIES FOLDER ".")


add_executable(MyApplication MyApplication.cpp WaterSurfaceShader.cpp WaterSurfaceMesh.cpp ${WaterSurfaceShaders_RCS})

target_link_libraries(MyApplication
  WaveGrid
  drawables
  Magnum::MeshTools
  Magnum::Primitives
  Magnum::Application
  Magnum::Magnum
  Magnum::Shaders
  Magnum::SceneGraph
  MagnumImGui)
