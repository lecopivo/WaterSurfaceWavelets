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

find_library(MAGNUM_IMGUI_LIB MagnumImGui)
find_library(IMGUI_LIB ImGui)

set(VIS_LIBS WaveGrid
  ${MAGNUM_IMGUI_LIB}
  ${IMGUI_LIB}
  drawables
  Magnum::MeshTools
  Magnum::Primitives
  Magnum::Application
  Magnum::Magnum
  Magnum::Shaders
  Magnum::SceneGraph)

add_executable(MyApplication MyApplication.cpp WaterSurfaceShader.cpp WaterSurfaceMesh.cpp ${WaterSurfaceShaders_RCS})
target_link_libraries(MyApplication ${VIS_LIBS})

add_executable(ProfileVisualization ProfileVisualization.cpp )
target_link_libraries(ProfileVisualization ${VIS_LIBS})






