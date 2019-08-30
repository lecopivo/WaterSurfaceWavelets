#pragma once

#include <Magnum/GL/Buffer.h>
#include <Magnum/DefaultFramebuffer.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Mesh.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/MeshTools/Interleave.h>
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
#include <Magnum/Trade/MeshData3D.h>

#include <variant>

#include "../base/SceneBase3D.h"

namespace Magnum {
namespace Drawables {

class DrawableMesh : public SceneBase3D::Object3D, public SceneGraph::Drawable3D {
public:
  // using Shader = std::variant<Shaders::Phong, Shaders::VertexColor3D,
  //                             Shaders::Flat3D, Shaders::MeshVisualizer>;
  using Shader = Shaders::Phong;

  struct VertexData {
    Vector3 position;
    Vector3 normal;
    Color4  color;
  };

protected:
  explicit DrawableMesh(SceneBase3D::Object3D *parent, SceneGraph::DrawableGroup3D *group);

public:
  template <class Fun> void setVertices(Fun fun) {
    std::vector<VertexData> newData = _data;

    for (size_t i = 0; i < newData.size(); i++) {
      fun(i, newData[i]);
    }
    recomputeNormals(newData);
    bindBuffers(newData);
  }

protected:
  void bindBuffers(std::vector<VertexData> const &data);

private:
  void draw(const Matrix4 &       transformationMatrix,
            SceneGraph::Camera3D &camera) override;

  void recomputeNormals(std::vector<VertexData> &data);

public:
  GL::Mesh   _mesh;
  GL::Buffer _vertexBuffer, _indexBuffer{Buffer::TargetHint::ElementArray};
  Shader _shader;

  std::vector<VertexData>  _data;
  std::vector<UnsignedInt> _indices;
};

class Plane : public DrawableMesh {
public:
  explicit Plane(SceneBase3D::Object3D *parent, SceneGraph::DrawableGroup3D *group,
                         int nx, int ny);
};

class Sphere : public DrawableMesh {
public:
  explicit Sphere(SceneBase3D::Object3D *parent, SceneGraph::DrawableGroup3D *group,
                          int rings, int segments);
};

class DrawableLine : public SceneBase3D::Object3D, public SceneGraph::Drawable3D {
public:
  using Shader = Shaders::VertexColor3D;

  struct VertexData {
    Vector3 position;
    Color4  color;
  };

public:
  explicit DrawableLine(SceneBase3D::Object3D *parent, SceneGraph::DrawableGroup3D *group,
                        int n);

  void resize(int n);

  template <class Fun> void setVertices(Fun fun) {
    std::vector<VertexData> newData = _data;

    for (size_t i = 0; i < newData.size(); i++) {
      fun(i, newData[i]);
    }
    bindBuffers(newData);
  }

protected:
  void bindBuffers(std::vector<VertexData> const &data);

private:
  void draw(const Matrix4 &       transformationMatrix,
            SceneGraph::Camera3D &camera) override;

public:
  GL::Mesh   _mesh;
  GL::Buffer _vertexBuffer;
  Shader _shader;

  std::vector<VertexData> _data;
};
} // namespace Drawables
} // namespace Magnum
