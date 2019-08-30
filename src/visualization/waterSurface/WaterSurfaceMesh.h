#pragma once

#include <Magnum/GL/Buffer.h>
#include <Magnum/DefaultFramebuffer.h>
#include <Magnum/Image.h>
#include <Magnum/ImageView.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Mesh.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/Icosphere.h>
#include <Magnum/Primitives/Plane.h>
#include <Magnum/Primitives/UVSphere.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/GL/Shader.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/GL/Texture.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/Trade/MeshData3D.h>

#include <iostream>

#include "../base/SceneBase3D.h"
#include "WaterSurfaceShader.h"
#include "../../ProfileBuffer.h"

namespace Magnum {

class WaterSurfaceMesh : public SceneBase3D::Object3D, public SceneGraph::Drawable3D {
public:
  struct VertexData {
    Vector3                      position;
    Math::Vector<DIR_NUM, Float> amplitude;
  };

public:
  explicit WaterSurfaceMesh(SceneBase3D::Object3D *                   parent,
                            SceneGraph::DrawableGroup3D *group, int n);

public:
  template <class Fun> void setVertices(Fun fun) {
  // std::vector<VertexData> newData = _data;

#pragma omp parallel for
    for (size_t i = 0; i < _data.size(); i++) {
      fun(i, _data[i]);
    }
    bindBuffers(_data);
  }

  void loadProfile(WaterWavelets::ProfileBuffer const& profileBuffer);

  void showTriangulationToggle();

protected:
  void bindBuffers(std::vector<VertexData> const &data);
  void bindTexture();

private:
  void draw(const Matrix4 &       transformationMatrix,
            SceneGraph::Camera3D &camera) override;

public:
  GL::Mesh   _mesh;
  GL::Buffer _vertexBuffer, _indexBuffer;
  bool _showTriangulation = false;

  Shaders::WaterSurfaceShader _shader;

  std::vector<VertexData>  _data;
  std::vector<UnsignedInt> _indices;

  GL::Texture1D _profileTexture;
};

} // namespace Magnum
