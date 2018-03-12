#include "VisualizationPrimitives.h"

#include <iostream>

namespace Magnum {

using namespace Magnum::Math::Literals;

DrawableMesh::DrawableMesh(Object3D *parent, SceneGraph::DrawableGroup3D *group)
    : Object3D{parent}, SceneGraph::Drawable3D{*this, group},
      _shader(Shaders::Phong{}) {

  auto &s = std::get<Shaders::Phong>(_shader);
  s.setDiffuseColor(Color4{0.4f, 0.4f, 0.8f, 1.f})
      .setAmbientColor(Color3{0.25f, 0.2f, 0.23f});

  using GShader = Shaders::Generic3D;

  _mesh.setPrimitive(MeshPrimitive::Triangles)
      .setCount(_indices.size())
      .addVertexBuffer(_vertexBuffer, 0, GShader::Position{}, GShader::Normal{},
                       GShader::Color{GShader::Color::Components::Four})
      .setIndexBuffer(_indexBuffer, 0, Mesh::IndexType::UnsignedInt);
}

void DrawableMesh::bindBuffers(std::vector<VertexData> const &data) {

  _vertexBuffer.setData(data, BufferUsage::DynamicDraw);
  _indexBuffer.setData(_indices, BufferUsage::StaticDraw);
  _mesh.setCount(_indices.size());
}

void DrawableMesh::draw(const Matrix4 &       transformationMatrix,
                        SceneGraph::Camera3D &camera) {

  if (std::holds_alternative<Shaders::Phong>(_shader)) {
    auto &s = std::get<Shaders::Phong>(_shader);

    s.setLightPosition(camera.cameraMatrix().transformPoint({5.0f, 5.0f, 7.0f}))
        .setTransformationMatrix(transformationMatrix)
        .setNormalMatrix(transformationMatrix.rotation())
        .setProjectionMatrix(camera.projectionMatrix());
    _mesh.draw(s);
  }

  if (std::holds_alternative<Shaders::VertexColor3D>(_shader)) {
    auto &s = std::get<Shaders::VertexColor3D>(_shader);

    s.setTransformationProjectionMatrix(camera.projectionMatrix() *
                                        transformationMatrix);
    _mesh.draw(s);
  }

  if (std::holds_alternative<Shaders::Flat3D>(_shader)) {
    auto &s = std::get<Shaders::Flat3D>(_shader);

    s.setTransformationProjectionMatrix(camera.projectionMatrix() *
                                        transformationMatrix);
    _mesh.draw(s);
  }

  if (std::holds_alternative<Shaders::MeshVisualizer>(_shader)) {
    auto &s = std::get<Shaders::MeshVisualizer>(_shader);

    s.setTransformationProjectionMatrix(camera.projectionMatrix() *
                                        transformationMatrix)
        .setViewportSize(Vector2{defaultFramebuffer.viewport().size()});
    _mesh.draw(s);
  }
}

void DrawableMesh::recomputeNormals(std::vector<VertexData> &data) {
  for (auto &d : data) {
    d.normal = Vector3{0.f, 0.f, 0.f};
  }

  for (size_t i = 0; i < _indices.size();) {

    auto id1 = _indices[i++];
    auto id2 = _indices[i++];
    auto id3 = _indices[i++];

    auto v1 = data[id1].position;
    auto v2 = data[id2].position;
    auto v3 = data[id3].position;

    // This does weighted area based on triangle area
    auto n = cross((v2 - v1), (v3 - v1));

    data[id1].normal += n;
    data[id2].normal += n;
    data[id3].normal += n;
  }

  for (auto &d : data)
    d.normal = d.normal.normalized();
}

DrawablePlane::DrawablePlane(Object3D *                   parent,
                             SceneGraph::DrawableGroup3D *group, int nx, int ny)
    : DrawableMesh(parent, group) {

  float dx = 2.0f / nx;
  float dy = 2.0f / ny;

  _data.clear();
  for (int i = 0; i <= nx; i++) {
    for (int j = 0; j <= ny; j++) {
      VertexData vertex;
      vertex.position = Vector3{-1.0f + i * dx, -1.0f + j * dy, 0.0f};
      vertex.normal   = Vector3{0.f, 0.f, 1.0f};
      vertex.color    = Color4{1.f, 1.f, 1.f, 1.f};

      _data.push_back(vertex);
    }
  }

  _indices.clear();
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      const int idx = j + i * (ny + 1);
      const int J   = 1;
      const int I   = ny + 1;

      _indices.push_back(idx);
      _indices.push_back(idx + I);
      _indices.push_back(idx + J);

      _indices.push_back(idx + I);
      _indices.push_back(idx + I + J);
      _indices.push_back(idx + J);
    }
  }

  bindBuffers(_data);
}

DrawableSphere::DrawableSphere(Object3D *                   parent,
                               SceneGraph::DrawableGroup3D *group, int rings,
                               int segments)
    : DrawableMesh(parent, group) {

  auto meshData = Primitives::UVSphere::solid(rings, segments);

  _data.resize(meshData.positions(0).size());
  for (size_t i = 0; i < meshData.positions(0).size(); i++) {
    _data[i].position = Matrix4::rotationX(Rad{M_PI / 2})
                            .transformPoint(meshData.positions(0)[i]);
    _data[i].normal = meshData.normals(0)[i];
  }

  _indices = meshData.indices();

  bindBuffers(_data);
}

// DrawableLine::DrawableLine(Object3D *parent, SceneGraph::DrawableGroup3D
// *group,
//                            int n)
//     : Object3D{parent}, SceneGraph::Drawable3D{*this, group} {
//   resize(n);
//   using Shader = Shaders::Generic3D;

//   _mesh.setPrimitive(MeshPrimitive::LineStrip)
//       .addVertexBuffer(_vertexBuffer, 0, Shader::Position{},
//                        Shader::Color{Shader::Color::Components::Four});
// }

// void DrawableLine::resize(int n) {
//   if (n != _data.size()) {
//     _data.resize(n);

//     float dx = 1.0 / (n - 1);
//     for (int i = 0; i < n; i++) {
//       _data[i].position = Vector3{i * dx, 0.f, 0.f};
//       _data[i].color    = Color4{1.f, 1.f, 1.f, 1.f};
//     }

//     bindBuffers(_data);
//   }
// }

// void DrawableLine::bindBuffers(std::vector<VertexData> const &data) {

//   _vertexBuffer.setData(data, BufferUsage::DynamicDraw);
//   _mesh.setCount(data.size());
// }

// void DrawableLine::draw(const Matrix4 &       transformationMatrix,
//                         SceneGraph::Camera3D &camera) {

//   _shader.setTransformationProjectionMatrix(camera.projectionMatrix() *
//                                             transformationMatrix);

//   _mesh.draw(_shader);
// }

} // namespace Magnum
