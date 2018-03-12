#include <Magnum/Buffer.h>
#include <Magnum/DefaultFramebuffer.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Mesh.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/Platform/Sdl2Application.h>
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
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Trade/MeshData3D.h>

#include <MagnumImGui/MagnumImGui.h>
#include <MagnumImGui/imgui.h>

#include <iostream>

#include "../WaveGrid.hpp"
#include "DrawingPrimitives.hpp"

namespace Magnum {

using namespace Magnum::Math::Literals;

typedef SceneGraph::Object<SceneGraph::MatrixTransformation3D> Object3D;
typedef SceneGraph::Scene<SceneGraph::MatrixTransformation3D>  Scene3D;

auto settings = []() {
  WaveGrid::Settings s;
  s.xmin = -50.0;
  s.xmax = +50.0;
  s.ymin = -50.0;
  s.ymax = +50.0;

  s.n_x     = 100;
  s.n_y     = 100;
  s.n_theta = 32;
  s.n_k     = 1;

  return s;
}();

Color4 gridAmplitudeColor(float x, float y, WaveGrid &grid) {
  Magnum::Color4 color = {0.f, 0.f, 0.f, 1.f};

  for (int b_theta = 0; b_theta < grid.gridDim(2); b_theta++) {
    Vec4 pos;
    Real theta = grid.idxToPos(b_theta, 2);
    Real k     = grid.idxToPos(0, 3);
    pos << x, y, theta, k;
    Real a = grid.amplitude(pos);

    color += Color4::fromHsv(Rad{(float)theta}, 1.f, a);
  }
  return color;
}

Vector3 simToWorld(Vec2 v) {
  float x = (v[0] - settings.xmin) / (settings.xmax - settings.xmin);
  float y = (v[1] - settings.ymin) / (settings.ymax - settings.ymin);
  x       = 2.0 * x - 1.0;
  y       = 2.0 * y - 1.0;
  return Vector3{x, y, 0};
}

Vec2 worldToSim(Vector3 v) {
  float x = 0.5 * (v[0] + 1.f);
  float y = 0.5 * (v[1] + 1.f);
  x       = x * (settings.xmax - settings.xmin) + settings.xmin;
  y       = y * (settings.ymax - settings.ymin) + settings.ymin;
  Vec2 out;
  out << x, y;
  return out;
}

class PrimitivesExample : public Platform::Application {
public:
  explicit PrimitivesExample(const Arguments &arguments);

private:
  void drawEvent() override;
  void drawGui();
  void update();

  void viewportEvent(const Vector2i &size) override;

  void keyPressEvent(KeyEvent &event) override;
  void keyReleaseEvent(KeyEvent &event) override;
  void mousePressEvent(MouseEvent &event) override;
  void mouseReleaseEvent(MouseEvent &event) override;
  void mouseMoveEvent(MouseMoveEvent &event) override;
  void mouseScrollEvent(MouseScrollEvent &event) override;
  void textInputEvent(TextInputEvent &event) override;

  void mouseRotation(MouseMoveEvent const &event, Vector2 delta);
  void mouseZoom(MouseMoveEvent const &event, Vector2 delta);
  void mousePan(MouseMoveEvent const &event, Vector2 delta);

  Vector3 mousePlanePosition(Vector2i mouseScreenPos);

  Scene3D                     _scene;
  Object3D *                  _cameraObject;
  SceneGraph::Camera3D *      _camera;
  SceneGraph::DrawableGroup3D _drawables;

  Vector2i _previousMousePosition;
  Vector3  _mousePlanePosition[2];

  MagnumImGui _gui;

  DrawablePlane * plane;
  DrawableSphere *sphere;
  DrawableLine *  line;

  WaveGrid _waveGrid;
  float    time = 0.f;

  Color3 _color;
};

PrimitivesExample::PrimitivesExample(const Arguments &arguments)
    : Platform::Application{arguments,
                            Configuration{}
                                .setTitle("Magnum object picking example")
                                .setWindowFlags(Sdl2Application::Configuration::
                                                    WindowFlag::Resizable)},
      _waveGrid(settings) {

  Renderer::enable(Renderer::Feature::DepthTest);
  Renderer::enable(Renderer::Feature::FaceCulling);

  //   /* Configure camera */
  _cameraObject = new Object3D{&_scene};
  _cameraObject->translate(Vector3::zAxis(4.0f)).rotateX(Rad{M_PI / 4});
  _camera = new SceneGraph::Camera3D{*_cameraObject};
  viewportEvent(defaultFramebuffer.viewport().size()); // set up camera

  /* TODO: Prepare your objects here and add them to the scene */
  sphere = (new DrawableSphere(&_scene, &_drawables, 10, 10));

  sphere->setVertices([settings](int i, DrawableSphere::VertexData &v) {
    v.color = Color4{1.f, 0.f, 0.f, 1.f};
    v.position *= 1.0 / 50.0f;
  });

  plane = (new DrawablePlane(&_scene, &_drawables, settings.n_x + 1,
                             settings.n_y + 1));

  line = (new DrawableLine(&_scene, &_drawables, 100));
}

void PrimitivesExample::update() {

  // times step simulation
  float dt = _waveGrid.cflTimeStep();
  _waveGrid.timeStep(dt, time);
  time += dt;

  // Set up wave grid visualization
  plane->setVertices([this](int i, DrawableMesh::VertexData &v) {
    Vec2 p  = worldToSim(v.position);
    v.color = gridAmplitudeColor(p.x(), p.y(), _waveGrid);
  });

  // Set up trajectory
  Vec4    pos4;
  Vector3 dir   = _mousePlanePosition[0] - _mousePlanePosition[1];
  Vec2    x0 = worldToSim(_mousePlanePosition[0]);
  double  theta = atan2(dir.y(), dir.x());
  double  k     = _waveGrid.idxToPos(0, 3);
  pos4 << x0.x(), x0.y(), theta, k;

  auto trajectory = _waveGrid.trajectory(pos4, settings.xmax-settings.xmin);
  line->resize(trajectory.size());
  line->setVertices([&](int i, DrawableLine::VertexData &v) {
    auto pos4     = trajectory[i];
    v.position    = simToWorld(trajectory[i].segment<2>(0));
    v.position[2] = 0.01f;
  });

  // mouse location visualization
  sphere->translate(
      _mousePlanePosition[0] -
      sphere->transformation().transformPoint(Vector3{0.f, 0.f, 0.f}));
}

void PrimitivesExample::drawEvent() {
  defaultFramebuffer.clear(FramebufferClear::Color | FramebufferClear::Depth);

  update();

  _camera->draw(_drawables);

  drawGui();

  swapBuffers();
}

void PrimitivesExample::drawGui() {
  _gui.newFrame(windowSize(), defaultFramebuffer.viewport().size());

  ImGui::ColorEdit3("Box color", &(_color[0]));

  _gui.drawFrame();

  redraw();
}

void PrimitivesExample::viewportEvent(const Vector2i &size) {
  defaultFramebuffer.setViewport({{}, size});

  _camera->setProjectionMatrix(Matrix4::perspectiveProjection(
      60.0_degf, Vector2{size}.aspectRatio(), 0.001f, 10000.0f));
}

void PrimitivesExample::keyPressEvent(KeyEvent &event) {
  if (_gui.keyPressEvent(event)) {
    redraw();
    return;
  }

  if (event.key() == KeyEvent::Key::Esc) {
    exit();
  }
}

void PrimitivesExample::keyReleaseEvent(KeyEvent &event) {
  if (_gui.keyReleaseEvent(event)) {
    redraw();
    return;
  }
}

void PrimitivesExample::mousePressEvent(MouseEvent &event) {
  if (_gui.mousePressEvent(event)) {
    redraw();
    return;
  }

  if (event.button() != MouseEvent::Button::Left)
    return;

  _previousMousePosition = event.position();
  event.setAccepted();
}

void PrimitivesExample::mouseReleaseEvent(MouseEvent &event) {
  if (_gui.mouseReleaseEvent(event)) {
    redraw();
    return;
  }

  event.setAccepted();
  redraw();
}

void PrimitivesExample::mouseMoveEvent(MouseMoveEvent &event) {
  if (_gui.mouseMoveEvent(event)) {
    redraw();
    return;
  }

  float lambda = 0.02;
  _mousePlanePosition[1] =
      lambda * _mousePlanePosition[0] + (1 - lambda) * _mousePlanePosition[1];
  _mousePlanePosition[0] = mousePlanePosition(event.position());

  const Vector2 delta = Vector2{event.position() - _previousMousePosition} /
                        Vector2{defaultFramebuffer.viewport().size()};

  if (event.buttons() & MouseMoveEvent::Button::Left)
    mouseRotation(event, delta);

  if (event.buttons() & MouseMoveEvent::Button::Right)
    mouseZoom(event, delta);

  if (event.buttons() & MouseMoveEvent::Button::Middle)
    mousePan(event, delta);

  _previousMousePosition = event.position();
  event.setAccepted();
  redraw();
}

void PrimitivesExample::mouseScrollEvent(MouseScrollEvent &event) {
  if (_gui.mouseScrollEvent(event)) {
    redraw();
    return;
  }
}

void PrimitivesExample::textInputEvent(TextInputEvent &event) {
  if (_gui.textInputEvent(event)) {
    redraw();
    return;
  }
}

void PrimitivesExample::mouseRotation(MouseMoveEvent const &event,
                                      Vector2               delta) {

  auto camPos =
      _cameraObject->transformation().transformPoint(Vector3{0.0, 0.0, 0.0});

  auto axis = cross(Vector3{0.f, 0.f, 1.f}, camPos.normalized()).normalized();

  _cameraObject->rotate(Rad{-3.0f * delta.y()}, axis);
  _cameraObject->rotateZ(Rad{-3.0f * delta.x()});
}

void PrimitivesExample::mouseZoom(MouseMoveEvent const &event, Vector2 delta) {
  auto dir =
      _cameraObject->transformation().transformVector(Vector3{0.0, 0.0, 1.0});

  _cameraObject->translate(30.0f * delta.y() * dir);
}

void PrimitivesExample::mousePan(MouseMoveEvent const &event, Vector2 delta) {}

Vector3 PrimitivesExample::mousePlanePosition(Vector2i mouseScreenPos) {
  Vector2 mpos = 2.0f * (Vector2{mouseScreenPos} /
                             Vector2{defaultFramebuffer.viewport().size()} -
                         Vector2{.5f, 0.5f});
  mpos[1] *= -1.f;

  Vector3 mpos3 = {mpos[0], mpos[1], -1.f};
  auto    trans =
      _cameraObject->transformation() * _camera->projectionMatrix().inverted();
  mpos3 = trans.transformPoint(mpos3);
  Vector3 camOrg =
      _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  Vector3 dir    = mpos3 - camOrg;
  float   lambda = -camOrg.z() / dir.z();
  mpos3          = camOrg + lambda * dir;

  return mpos3;
}

} // namespace Magnum

MAGNUM_APPLICATION_MAIN(Magnum::PrimitivesExample)
