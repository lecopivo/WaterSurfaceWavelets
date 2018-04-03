#include "WaterSurfaceMesh.h"

#include <iostream>

namespace Magnum {

using namespace Magnum::Math::Literals;

WaterSurfaceMesh::WaterSurfaceMesh(SceneBase3D::Object3D *      parent,
                                   SceneGraph::DrawableGroup3D *group, int n)
    : SceneBase3D::Object3D{parent}, SceneGraph::Drawable3D{*this, group} {

  using Shader = Shaders::WaterSurfaceShader;

  _shader.setColor(Color4{0.4f, 0.4f, 0.8f, 1.f})
      .setAmbientColor(Color3{0.25f, 0.2f, 0.23f});

  int   nx = n;
  int   ny = n;
  float dx = 2.0f / nx;
  float dy = 2.0f / ny;

  _data.clear();
  for (int i = 0; i <= nx; i++) {
    for (int j = 0; j <= ny; j++) {
      VertexData vertex;
      vertex.position = Vector3{-1.0f + i * dx, -1.0f + j * dy, 0.0f};
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

  _vertexBuffer.setData(_data, BufferUsage::DynamicDraw);
  _indexBuffer.setData(_indices, BufferUsage::StaticDraw);

  _mesh.setPrimitive(MeshPrimitive::Triangles)
      .setCount(_indices.size())
      .addVertexBuffer(_vertexBuffer, 0, Shader::Position{})
      .setIndexBuffer(_indexBuffer, 0, Mesh::IndexType::UnsignedInt);
}

void WaterSurfaceMesh::loadProfile(
    WaterWavelets::ProfileBuffer const &profileBuffer) {

  Containers::ArrayView<const void> data(profileBuffer.m_data.data(),
                                         profileBuffer.m_data.size() *
                                             sizeof(std::array<float, 4>));

  ImageView1D image(PixelFormat::RGBA, PixelType::Float,
                    profileBuffer.m_data.size(), data);

  _profileTexture.setWrapping(Sampler::Wrapping::Repeat)
      .setMagnificationFilter(Sampler::Filter::Linear)
      .setMinificationFilter(Sampler::Filter::Linear)
      .setStorage(1, TextureFormat::RGBA32F, profileBuffer.m_data.size())
      .setSubImage(0, {}, image);

  _shader.bindProfileTexture(_profileTexture)
      .setProfilePeriod(profileBuffer.m_period);
}

void WaterSurfaceMesh::loadAmplitude(WaterWavelets::WaveGrid const &grid) {

  Vector3i size = {grid.gridDim(2), grid.gridDim(1), grid.gridDim(0)};

  assert(grid.gridDim(3) ==
         1); // Right now this works only for discretization in one wave number

  Containers::ArrayView<const void> data(grid.m_amplitude.m_data.data(),
                                         size[0] * size[1] * size[2] *
                                             sizeof(float));

  ImageView3D image(PixelFormat::Red, PixelType::Float, size, data);

  _amplitudeTexture
      .setWrapping({Sampler::Wrapping::Repeat, Sampler::Wrapping::ClampToBorder,
                    Sampler::Wrapping::ClampToBorder})
      .setMagnificationFilter(Sampler::Filter::Linear)
      .setMinificationFilter(Sampler::Filter::Linear)
      .setStorage(1, TextureFormat::R32F, size)
      .setSubImage(0, {}, image);

  _shader.bindAmplitudeTexture(_amplitudeTexture)
      .setDomainSize(grid.m_xmax[0])
      .setDirectionNumber(grid.gridDim(2));
}

void WaterSurfaceMesh::showTriangulationToggle() {
  _showTriangulation = !_showTriangulation;
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
