# Right now I compile all resource directly to the executable
# I should use: http://doc.magnum.graphics/corrade/Resource_8h.html#afe3aa100fc67ff4c37d00436c194f858
# corrade_add_resource(WaterSurfaceShaders_RCS waterSurface/resources.conf)

add_library(WaterSurface WaterSurfaceShader.cpp WaterSurfaceMesh.cpp ) # ${WaterSurfaceShaders_RCS}
