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

class Scene {
public:
  Scene() : grid(initSettings()) {}

  WaveGrid::Settings initSettings() {
    WaveGrid::Settings s;
    s.xmin = -50;
    s.xmax = 50;
    s.ymin = -50;
    s.ymax = 50;

    s.n_x = 100;
    s.n_y = 100;
    s.n_theta = 8;
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

    Color3 out{0,0,0};

    Vec4 pos4;
    double k = grid.waveNumber(0);

    pos4 << x, y, 0, k;

    for (int c_theta = 0; c_theta < grid.gridDim(2); c_theta++) {
      double theta = grid.waveVectorAngle(c_theta);

      pos4(2) = theta;
#error I stoped here!
      out += grid.amplitude(pos4);
    }
  }

  WaveGrid grid;
  double t = 0.0;
};

template <class PlaneVertex, class Shader> class DrawingPlane {
public:
  DrawingPlane() {}

  void init(float _lx, float _ly, int _nx, int _ny) {
    lx = _lx;
    ly = _ly;
    nx = _nx;
    ny = _ny;
    dx = _lx / _nx;
    dy = _ly / _ny;

    vertexData.reserve((nx + 1) * (ny + 1));
    indexData.reserve(6 * nx * ny);

    for (int i = 0; i < nx + 1; i++) {
      for (int j = 0; j < ny + 1; j++) {
        float x = -0.5 * lx + i * dx;
        float y = -0.5 * ly + j * dy;

        vertexData.push_back(PlaneVertex(x, y));

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
  template <class Fun> void changeVertexData(Fun fun) {
    for (int i = 0; i < nx + 1; i++) {
      for (int j = 0; j < ny + 1; j++) {
        float x = -0.5 * lx + i * dx;
        float y = -0.5 * ly + j * dy;

        fun(x, y, vertexData[j + i * (ny + 1)]);
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

  std::vector<PlaneVertex> vertexData;
  std::vector<UnsignedInt> indexData;

  Buffer vertexBuffer, indexBuffer;
  Mesh mesh;
  Shader shader;

  float lx, ly, nx, ny, dx, dy;
};

class TriangleExample : public Platform::Application {
public:
  explicit TriangleExample(const Arguments &arguments);

private:
  void drawEvent() override;

  struct PlaneVertex {
    PlaneVertex(float x, float y)
        : position(x, y), color(sin(10 * x), cos(10 * y), 1.0) {}

    Vector2 position;
    Color3 color;
  };

  DrawingPlane<PlaneVertex, Shaders::VertexColor2D> drawingPlane;

  Scene scene;
};

TriangleExample::TriangleExample(const Arguments &arguments)
    : Platform::Application{
          arguments, Configuration{}
                         .setTitle("Magnum Triangle Example")
                         .setWindowFlags(Configuration::WindowFlag::Resizable)
                         .setSize({800, 800})} {
  drawingPlane.init(2.0, 2.0, 800, 800);
}

float time = 0.0;

void TriangleExample::drawEvent() {
  defaultFramebuffer.clear(FramebufferClear::Color);

  // drawingPlane.shader.setTransformationProjectionMatrix(Matrix3::scaling(1.0,));
  drawingPlane.draw();

  time += 0.1;

  drawingPlane.changeVertexData([](float x, float y, PlaneVertex &vertex) {
    vertex.color = Color3{sin(time + 10 * x), cos(time + 10 * y), 1.0};
  });

  swapBuffers();

  redraw();
}

} // namespace Examples
} // namespace Magnum

MAGNUM_APPLICATION_MAIN(Magnum::Examples::TriangleExample)
