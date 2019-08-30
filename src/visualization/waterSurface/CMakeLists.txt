add_library(WaterSurface WaterSurfaceShader.cpp WaterSurfaceMesh.cpp ) # ${WaterSurfaceShaders_RCS}
target_include_directories(WaterSurface PUBLIC ${MAGNUM_INCLUDE_DIR})
