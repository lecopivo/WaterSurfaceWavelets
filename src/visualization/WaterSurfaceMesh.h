#pragma once

#include <Magnum/Buffer.h>
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
#include <Magnum/Renderer.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/Shader.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Texture.h>
#include <Magnum/TextureFormat.h>
#include <Magnum/Trade/MeshData3D.h>

#include <iostream>

#include "WaterSurfaceShader.h"

namespace Magnum {

using Object3D = SceneGraph::Object<SceneGraph::MatrixTransformation3D>;
using Scene3D  = SceneGraph::Scene<SceneGraph::MatrixTransformation3D>;

class WaterSurfaceMesh : public Object3D, public SceneGraph::Drawable3D {
public:
  struct VertexData {
    Vector3                      position;
    Math::Vector<DIR_NUM, Float> amplitude;
  };

public:
  explicit WaterSurfaceMesh(Object3D *                   parent,
                            SceneGraph::DrawableGroup3D *group, int n);

public:
  template <class Fun> void setVertices(Fun fun) {
    std::vector<VertexData> newData = _data;

    for (size_t i = 0; i < newData.size(); i++) {
      fun(i, newData[i]);
    }
    bindBuffers(newData);
  }

  template <class Fun> void setHeightField(Fun fun) {

    float dx = 1.0f / _heightTextureSize;

    for (int i = 0; i < _heightTextureSize; i++) {
      float x = (i + 0.5f) * dx;
      fun(x, _heightData[i]);
    }

    bindTexture();
  }

protected:
  void bindBuffers(std::vector<VertexData> const &data);
  void bindTexture();

private:
  void draw(const Matrix4 &       transformationMatrix,
            SceneGraph::Camera3D &camera) override;

public:
  Mesh   _mesh;
  Buffer _vertexBuffer, _indexBuffer;

  Shaders::WaterSurfaceShader _shader;

  std::vector<VertexData>  _data;
  std::vector<UnsignedInt> _indices;

  const int            _heightTextureSize = 4096;
  std::vector<Vector4> _heightData{size_t{_heightTextureSize}};
  Texture1D            _heightTexture;
};

} // namespace Magnum
