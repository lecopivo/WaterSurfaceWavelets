#include "WaterSurfaceMesh.h"

#include <iostream>

namespace Magnum {

using namespace Magnum::Math::Literals;

template <typename Mesh, typename Buffer, size_t... I>
void addVertexBuffer(Mesh &mesh, Buffer &buffer, std::index_sequence<I...>) {

  using Shader = Shaders::WaterSurfaceShader;

  mesh.addVertexBuffer(buffer, 0, Shader::Position{},
                       Shader::Amplitude<I>{}...);
}

WaterSurfaceMesh::WaterSurfaceMesh(SceneBase3D::Object3D *                   parent,
                                   SceneGraph::DrawableGroup3D *group, int n)
    : SceneBase3D::Object3D{parent}, SceneGraph::Drawable3D{*this, group} {

  _shader.setColor(Color4{0.4f, 0.4f, 0.8f, 1.f})
      .setAmbientColor(Color3{0.25f, 0.2f, 0.23f})
      .setTime(0)
      .setGerstnerParameter(1)
      .setWaveDirectionToShow(-1);

  _mesh.setPrimitive(MeshPrimitive::Triangles)
      .setCount(_indices.size())
      .setIndexBuffer(_indexBuffer, 0, Mesh::IndexType::UnsignedInt);

  addVertexBuffer(_mesh, _vertexBuffer,
                  std::make_index_sequence<DIR_NUM / 4>{});

  int   nx = n;
  int   ny = n;
  float dx = 2.0f / nx;
  float dy = 2.0f / ny;

  _data.clear();
  for (int i = 0; i <= nx; i++) {
    for (int j = 0; j <= ny; j++) {
      VertexData vertex;
      vertex.position = Vector3{-1.0f + i * dx, -1.0f + j * dy, 0.0f};
      for (int k = 0; k < DIR_NUM; k++)
        vertex.amplitude[k] = 0;
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

void WaterSurfaceMesh::loadProfile(
    WaterWavelets::ProfileBuffer const &profileBuffer) {

  Containers::ArrayView<const void> data(profileBuffer.m_data.data(),
                                         profileBuffer.m_data.size() *
                                             sizeof(std::array<float, 4>));

  ImageView1D image(GL::PixelFormat::RGBA, PixelType::Float, profileBuffer.m_data.size(),
                    data);

  _profileTexture.setWrapping(Sampler::Wrapping::Repeat)
      .setMagnificationFilter(Sampler::Filter::Linear)
      .setMinificationFilter(Sampler::Filter::Linear)
      .setStorage(1, TextureFormat::RGBA32F, profileBuffer.m_data.size())
      .setSubImage(0, {}, image);

  _shader.bindTexture(_profileTexture).setProfilePeriod(profileBuffer.m_period);
}

void WaterSurfaceMesh::showTriangulationToggle() {
  _showTriangulation = !_showTriangulation;
}

void WaterSurfaceMesh::bindBuffers(std::vector<VertexData> const &data) {

  _vertexBuffer.setData(data, BufferUsage::DynamicDraw);
  _indexBuffer.setData(_indices, BufferUsage::StaticDraw);
  _mesh.setCount(_indices.size());
}

void WaterSurfaceMesh::draw(const Matrix4 &       transformationMatrix,
                            SceneGraph::Camera3D &camera) {

  _shader
      .setLightPosition(
          camera.cameraMatrix().transformPoint({15.0f, 15.0f, 30.0f}))
      //    .setLightPosition(Vector3{0.0,5.0,5.0})
      .setTransformationMatrix(transformationMatrix)
      .setNormalMatrix(transformationMatrix.rotation())
      .setProjectionMatrix(camera.projectionMatrix());

  if (_showTriangulation)
    Renderer::setPolygonMode(Renderer::PolygonMode::Line);
  _mesh.draw(_shader);
  if (_showTriangulation)
    Renderer::setPolygonMode(Renderer::PolygonMode::Fill);
}

} // namespace Magnum
