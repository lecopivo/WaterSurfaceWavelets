#  _     _ _         
# | |   (_) |__  ___ 
# | |   | | '_ \/ __|
# | |___| | |_) \__ \
# |_____|_|_.__/|___/
#


#Find Magnum 
find_package(Magnum REQUIRED
  GL
  MeshTools
  Primitives
  Shaders
  SceneGraph
  Sdl2Application)

find_package(MagnumIntegration CONFIG REQUIRED ImGui)

message("ImGui Found: ${ImGUi_FOUND}")
message("Integration Found: ${MagnumIntegration_ImGui_FOUND}")

# include_directories(${SDL2_INCLUDE_DIR}) # Why do I have to have this line?

add_subdirectory(base)
add_subdirectory(waterSurface)
add_subdirectory(drawables)

# Set up library to link to 
set(VIS_LIBS
  Magnum::GL
  Magnum::MeshTools
  Magnum::Primitives
  Magnum::Application
  Magnum::Magnum
  Magnum::Shaders
  Magnum::SceneGraph
  MagnumIntegration::ImGui
  WaveGrid
  Drawables
  WaterSurface
  SceneBase3D)

#  _____  _    ____   ____ _____ _____ ____  
# |_   _|/ \  |  _ \ / ___| ____|_   _/ ___| 
#   | | / _ \ | |_) | |  _|  _|   | | \___ \ 
#   | |/ ___ \|  _ <| |_| | |___  | |  ___) |
#   |_/_/   \_\_| \_\\____|_____| |_| |____/ 
#                                           

corrade_add_resource(WaterSurfaceShaders_RCS waterSurface/resources.conf)

# Build Demo
add_executable(Demo Demo.cpp ${WaterSurfaceShaders_RCS})
target_link_libraries(Demo ${VIS_LIBS})

# Build Profile Visualizer
# add_executable(ProfileVisualization ProfileVisualization.cpp )
# target_link_libraries(ProfileVisualization ${VIS_LIBS})




