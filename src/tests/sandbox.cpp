#include <Magnum/Buffer.h>
#include <Magnum/DefaultFramebuffer.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Mesh.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/Shaders/VertexColor.h>

#include "../WaveGrid.hpp"
#include <iostream>

namespace Magnum {
namespace Examples {

struct ColorVertex {
  Vector2 position;
  Color3 color;
};

template <class LineVertex, class Shader> class DrawingLine {
public:
  DrawingLine(int n) {
    vertexData.resize(n + 1);

    vertexBuffer.setData(vertexData, BufferUsage::DynamicDraw);
    mesh.setPrimitive(MeshPrimitive::LineStrip)
        .setCount(n + 1)
        .addVertexBuffer(vertexBuffer, 0, typename Shader::Position{},
                         typename Shader::Color{
                             Shaders::VertexColor2D::Color::Components::Three});
  }

  template <class Fun> void setVertexData(Fun fun) {
    for (int i = 0; i < vertexData.size(); i++) {
      fun(i, vertexData[i]);
    }
    vertexBuffer.setData(vertexData, BufferUsage::DynamicDraw);
  }

  void draw() { mesh.draw(shader); }

public:
  std::vector<LineVertex> vertexData;
  Buffer vertexBuffer;
  Mesh mesh;
  Shader shader;
};

template <class Vertex, class Shader> class DrawingPlane {
public:
  DrawingPlane(int _nx, int _ny) {
    nx = _nx;
    ny = _ny;

    vertexData.resize((nx + 1) * (ny + 1));
    indexData.clear();
    indexData.reserve(6 * nx * ny);

    for (int i = 0; i < nx + 1; i++) {
      for (int j = 0; j < ny + 1; j++) {

        if (i < nx && j < ny) {
          indexData.push_back(j + i * (ny + 1));
          indexData.push_back((j + 1) + i * (ny + 1));
          indexData.push_back(j + (i + 1) * (ny + 1));

          indexData.push_back(j + (i + 1) * (ny + 1));
          indexData.push_back((j + 1) + (i + 1) * (ny + 1));
          indexData.push_back((j + 1) + i * (ny + 1));
        }
      }
    }

    vertexBuffer.setData(vertexData, BufferUsage::DynamicDraw);
    indexBuffer.setData(indexData, BufferUsage::StaticDraw);
    mesh.setPrimitive(MeshPrimitive::Triangles)
        .setCount(6 * nx * ny)
        .addVertexBuffer(vertexBuffer, 0, typename Shader::Position{},
                         typename Shader::Color{
                             Shaders::VertexColor2D::Color::Components::Three})
        .setIndexBuffer(indexBuffer, 0, Mesh::IndexType::UnsignedInt);
  }

  void draw() { mesh.draw(shader); }

  //  Fun: void(float x, float y, PlaneVertex & vertex)
  template <class Fun> void setVertexData(Fun fun) {

    const int NX = nx;
    const int NY = ny;

#pragma omp parallel for collapse(2)
    for (int i = 0; i <= NX; i++) {
      for (int j = 0; j <= NY; j++) {
        fun(i, j, vertexData[j + i * (ny + 1)]);
      }
    }
    vertexBuffer.setData(vertexData, BufferUsage::DynamicDraw);
  }

  std::vector<Vertex> vertexData;
  std::vector<UnsignedInt> indexData;

  Buffer vertexBuffer, indexBuffer;
  Mesh mesh;
  Shader shader;

  float nx, ny;
};

class Scene {
public:
  Scene() : grid(initSettings()) {

    for (int b_theta = 0; b_theta < grid.gridDim(2); b_theta++) {
      grid.m_amplitude(20, 50, b_theta, 0) = 100.0;
    }

    grid.islandCenter << 20.0, 0.0;
    grid.islandRadius = 20.0;
  }

  WaveGrid::Settings initSettings() {
    WaveGrid::Settings s;
    s.xmin = -50;
    s.xmax = 50;
    s.ymin = -50;
    s.ymax = 50;

    s.n_x = 200;
    s.n_y = 200;
    s.n_theta = 32;
    s.n_k = 1;

    s.spectrumType = WaveGrid::Settings::HAHA_SPECTRUM;

    return s;
  }

  void timeStep(float lambda) {
    double dt = lambda * grid.cflTimeStep();
    t += dt;
    grid.timeStep(dt, t);
  }

  Color3 amplitudeColor(float x, float y) {

    Color3 out = {0, 0, 0};

    Vec4 pos4;
    double k = grid.waveNumber(0);

    pos4 << x, y, 0, k;

    for (int b_theta = 0; b_theta < grid.gridDim(2); b_theta++) {
      float theta = grid.waveVectorAngle(b_theta);

      pos4(2) = theta;
      out = out +
            Color3::fromHsv(Deg{(float)(360.0 * (theta + 3.0) / (2.0 * M_PI))},
                            1.0, grid.amplitude(pos4));
    }
    return out;
  }

  void drawAmplitude() {
    DrawingPlane<ColorVertex, Shaders::VertexColor2D> drawingPlane(
        grid.gridDim(0), grid.gridDim(1));

    drawingPlane.setVertexData([this](int i, int j, ColorVertex &vertex) {
      Vec2 pos = grid.nodePosition(i, j);
      vertex.position = Vector2{(float)pos(0), (float)pos(1)};
      vertex.color = amplitudeColor(pos(0), pos(1));
    });
    drawingPlane.shader.setTransformationProjectionMatrix(
        Matrix3::scaling({0.02, 0.02}));
    drawingPlane.draw();
  }

  void drawTrajectory(Vec2 pos, float theta) {
    int n = 100;
    DrawingLine<ColorVertex, Shaders::VertexColor2D> drawingLine(n);

    Vec4 pos4;
    pos4 << pos(0), pos(1), theta, grid.waveNumber(0);

    float dt = 0.8 * grid.cflTimeStep();

    drawingLine.setVertexData([&pos4, dt, this](int i, ColorVertex &vertex) {
      vertex.position = {(float)pos4(0), (float)pos4(1)};
      vertex.color = Color3{1.0f, 1.0f, 1.0f};

      Vec2 cg = grid.groupVelocity(pos4(2), pos4(3));

      pos4.segment<2>(0) += dt * cg;
      pos4 = grid.boundaryReflection(pos4);
    });
    drawingLine.shader.setTransformationProjectionMatrix(
        Matrix3::scaling({0.02, 0.02}));
    drawingLine.draw();
  }

public:
  WaveGrid grid;
  double t = 0.0;
};

class TriangleExample : public Platform::Application {
public:
  explicit TriangleExample(const Arguments &arguments);

private:
  void drawEvent() override;

  Scene scene;
};

TriangleExample::TriangleExample(const Arguments &arguments)
    : Platform::Application{
          arguments, Configuration{}
                         .setTitle("Magnum Triangle Example")
                         .setWindowFlags(Configuration::WindowFlag::Resizable)
                         .setSize({1100, 1100})} {}

float time = 0.0;

void TriangleExample::drawEvent() {
  defaultFramebuffer.clear(FramebufferClear::Color);

  time += 0.1;


  scene.drawAmplitude();

  // for(float x=-20.0;x<=20.0;x+=1.0){
  //   Vec2 pos; pos << -50.0,x+0.1;
  //   scene.drawTrajectory(pos,0.0);
  // }
  
  scene.timeStep(1.0);

  swapBuffers();

  redraw();
}

} // namespace Examples
} // namespace Magnum

MAGNUM_APPLICATION_MAIN(Magnum::Examples::TriangleExample)
