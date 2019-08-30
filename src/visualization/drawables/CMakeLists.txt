add_library(Drawables Primitives3D.cpp)
target_include_directories(Drawables PUBLIC ${MAGNUM_INCLUDE_DIR})
target_link_libraries(Drawables PUBLIC
  Magnum::GL
  Magnum::Primitives
  Magnum::Shaders)



