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

#include <MagnumImGui.h>
#include <imgui.h>

#include <algorithm>
#include <iostream>
#include <tuple>

#include "WaterSurfaceMesh.h" 1
#include "WaterSurfaceShader.h"
#include "drawables/VisualizationPrimitives.h"

#include "../WaveGrid.h"

using namespace Magnum;
using namespace Magnum::Math::Literals;
using namespace WaterWavelets;

using Object3D = SceneGraph::Object<SceneGraph::MatrixTransformation3D>;
using Scene3D  = SceneGraph::Scene<SceneGraph::MatrixTransformation3D>;

auto settings = []() {
  WaveGrid::Settings s;

  s.size      = 50;
  s.windSpeed = 20;

  s.n_x     = 100;
  s.n_theta = DIR_NUM;
  s.n_zeta  = 1;

  s.spectrumType = WaveGrid::Settings::PiersonMoskowitz;

  s.initial_time = 100;

  return s;
}();

constexpr float pi = 3.14159265359f;

struct CameraParameters {
  Vector3 target         = {0.f, 0.f, 0.f};
  float   longitude      = pi / 4;
  float   latitude       = pi / 4;
  float   targetDistance = 100.0f;

  Matrix4 getCameraTransformation() {
    Matrix4 trans = Matrix4::translation(target) *
                    Matrix4::rotationZ(Rad{longitude}) *
                    Matrix4::rotationX(-Rad{latitude}) *
                    Matrix4::translation(Vector3{0.f, -targetDistance, 0.f}) *
                    Matrix4::rotationX(Rad{pi / 2});
    return trans;
  }
};

template <int N, typename T> void printv(Magnum::Math::Vector<N, T> v) {
  for (int i = 0; i < N; i++) {
    std::cout << v[i] << " ";
  }
  std::cout << std::endl;
}

class MyApplication : public Platform::Application {
public:
  explicit MyApplication(const Arguments &arguments);

private:
  void drawEvent() override;
  void drawGui();

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

  /*
   * Sends a ray from a camera pixel
   *
   * \param pixel pixel to send a ray from
   * \return direction vector and camera position(in that order)
   */
  std::tuple<Vector3, Vector3> cameraRayCast(Vector2i pixel) const;
  std::tuple<Vector3, Vector3> cameraRayCast(Vector2 mouseScreenPos) const;

  Scene3D                     _scene;
  Object3D *                  _cameraObject;
  SceneGraph::Camera3D *      _camera;
  SceneGraph::DrawableGroup3D _drawables;
  CameraParameters            _cameraParams;

  Vector2i _previousMousePosition;

  MagnumImGui _gui;

  // Example objects to draw
  DrawablePlane * plane;
  DrawableSphere *sphere;
  // DrawableLine *  line;

  WaterSurfaceMesh *water_surface;

  // Stokes wave
  float logdt              = 0.1;
  float amplitude          = 3.0;
  float gerstnerParameter  = 1.0;
  int   directionToShow    = -1;
  int   gridResolution     = 200;
  bool  update_screen_grid = true;

  WaveGrid _waveGrid;
};

MyApplication::MyApplication(const Arguments &arguments)
    : Platform::Application{arguments,
                            Configuration{}
                                .setTitle("Magnum object picking example")
                                .setWindowFlags(Sdl2Application::Configuration::
                                                    WindowFlag::Resizable)},
      _waveGrid(settings) {
  /* Configure OpenGL state */
  Renderer::enable(Renderer::Feature::DepthTest);
  Renderer::enable(Renderer::Feature::FaceCulling);

  /* Configure camera */
  _cameraObject = new Object3D{&_scene};
  _cameraObject->setTransformation(_cameraParams.getCameraTransformation());
  _camera = new SceneGraph::Camera3D{*_cameraObject};
  viewportEvent(defaultFramebuffer.viewport().size()); // set up camera

  /* Set up object to draw */
  sphere        = new DrawableSphere(&_scene, &_drawables, 10, 10);
  plane         = new DrawablePlane(&_scene, &_drawables, 200, 200);
  water_surface = new WaterSurfaceMesh(&_scene, &_drawables, gridResolution);

  plane->setVertices([&](int, DrawablePlane::VertexData &v) {
    v.position *= 50;
    Vec2 p{v.position[0], v.position[1]};
    v.position[2] = -_waveGrid.m_enviroment.levelset(p);
  });
  std::get<Shaders::Phong>(plane->_shader)
      .setDiffuseColor(Color4{0.4, 0.4, 0.4, 1.0})
      .setAmbientColor(Color3{0.25f, 0.2f, 0.23f})
      .setShininess(10)
      .setSpecularColor(Color4{0.2, 0.2, 0.2, 1.0});

  std::get<Shaders::Phong>(sphere->_shader)
      .setDiffuseColor(Color4{0.8, 0.2, 0.2, 1.0})
      .setAmbientColor(Color3{0.8f, 0.2f, 0.23f})
      .setSpecularColor(Color4{0.2, 0.2, 0.2, 1.0});
}

void MyApplication::drawEvent() {
  defaultFramebuffer.clear(FramebufferClear::Color | FramebufferClear::Depth);

  if (update_screen_grid) {
    water_surface->setVertices([&](int i, WaterSurfaceMesh::VertexData &v) {

      int     ix        = i / (gridResolution + 1);
      int     iy        = i % (gridResolution + 1);
      Vector2 screenPos = Vector2{(2.0f * ix) / gridResolution - 1.0f,
                                  (2.0f * iy) / gridResolution - 1.0f};

      auto[dir, camPos] = cameraRayCast(screenPos);
      dir               = dir.normalized();
      float t           = -camPos.z() / dir.z();
      t                 = t < 0 ? 1000 : t;
      v.position        = camPos + t * dir;
      v.position.z()    = 0;

      for (int itheta = 0; itheta < DIR_NUM; itheta++) {
        float theta = _waveGrid.idxToPos(itheta, WaveGrid::Theta);
        // float theta = 2*pi*itheta/DIR_NUM;
        Vec4 pos4{v.position.x(), v.position.y(), theta,
                  _waveGrid.idxToPos(0, WaveGrid::Zeta)};
        v.amplitude[itheta] = amplitude * _waveGrid.amplitude(pos4);
      }
    });
  }

  water_surface->loadProfile(_waveGrid.m_profileBuffers[0]);

  _waveGrid.timeStep(_waveGrid.cflTimeStep() * pow(10, logdt));

  _camera->draw(_drawables);

  drawGui();

  swapBuffers();
}

void MyApplication::drawGui() {
  _gui.newFrame(windowSize(), defaultFramebuffer.viewport().size());

  ImGui::SliderFloat("amplitude", &amplitude, 0, 10);
  ImGui::SliderFloat("log10(dt)", &logdt, -3, 3);
  if (ImGui::SliderInt("direction", &directionToShow, -1, DIR_NUM)) {
    water_surface->_shader.setWaveDirectionToShow(directionToShow);
  }
  if (ImGui::SliderFloat("Gerstner wave", &gerstnerParameter, 0, 1)) {
    water_surface->_shader.setGerstnerParameter(gerstnerParameter);
  }
  ImGui::Checkbox("update screen grid", &update_screen_grid);
  if (ImGui::Button("Show triangulation"))
    water_surface->showTriangulationToggle();

  _gui.drawFrame();

  redraw();
}

void MyApplication::viewportEvent(const Vector2i &size) {
  defaultFramebuffer.setViewport({{}, size});

  _camera->setProjectionMatrix(Matrix4::perspectiveProjection(
      60.0_degf, Vector2{size}.aspectRatio(), 0.1f, 1000.0f));
}

void MyApplication::keyPressEvent(KeyEvent &event) {
  if (_gui.keyPressEvent(event)) {
    redraw();
    return;
  }

  if (event.key() == KeyEvent::Key::Esc) {
    exit();
  }

  if (event.key() == KeyEvent::Key::F) {
    _cameraParams.target = Vector3{0.f, 0.f, 0.f};
    _cameraObject->setTransformation(_cameraParams.getCameraTransformation());
  }

  redraw();
}

void MyApplication::keyReleaseEvent(KeyEvent &event) {
  if (_gui.keyReleaseEvent(event)) {
    redraw();
    return;
  }
}

void MyApplication::mousePressEvent(MouseEvent &event) {
  if (_gui.mousePressEvent(event)) {
    redraw();
    return;
  }

  if (event.button() == MouseEvent::Button::Left) {
    _previousMousePosition = event.position();
    event.setAccepted();
  }
}

void MyApplication::mouseReleaseEvent(MouseEvent &event) {
  if (_gui.mouseReleaseEvent(event)) {
    redraw();
    return;
  }

  event.setAccepted();
  redraw();
}

void MyApplication::mouseMoveEvent(MouseMoveEvent &event) {
  if (_gui.mouseMoveEvent(event)) {
    redraw();
    return;
  }

  auto[dir, camPos] = cameraRayCast(event.position());
  double  t         = -camPos.z() / dir.z();
  Vector3 spherePos = sphere->transformation().transformPoint(Vector3{0, 0, 0});
  Vector3 sphereNewPos = (camPos + t * dir);
  sphere->translate(sphereNewPos - spherePos);
  sphere->scale({0.01, 0.01, 0.01});
  Vec2 pos{sphereNewPos.x(), sphereNewPos.y()};
  if ((event.buttons() & MouseMoveEvent::Button::Left) &&
      !(event.modifiers() & MouseMoveEvent::Modifier::Alt)) {
    _waveGrid.addPointDisturbance(pos, 0.1);
  }

  const Vector2 delta = Vector2{event.position() - _previousMousePosition} /
                        Vector2{defaultFramebuffer.viewport().size()};

  if ((event.modifiers() & MouseMoveEvent::Modifier::Alt) &&
      (event.buttons() & MouseMoveEvent::Button::Left))
    mouseRotation(event, delta);

  if (event.modifiers() & MouseMoveEvent::Modifier::Alt &&
      event.buttons() & MouseMoveEvent::Button::Right)
    mouseZoom(event, delta);

  if (event.modifiers() & MouseMoveEvent::Modifier::Alt &&
      event.buttons() & MouseMoveEvent::Button::Middle)
    mousePan(event, delta);

  _previousMousePosition = event.position();
  event.setAccepted();
  redraw();
}

void MyApplication::mouseScrollEvent(MouseScrollEvent &event) {
  if (_gui.mouseScrollEvent(event)) {
    redraw();
    return;
  }
}

void MyApplication::textInputEvent(TextInputEvent &event) {
  if (_gui.textInputEvent(event)) {
    redraw();
    return;
  }
}

void MyApplication::mouseRotation(MouseMoveEvent const &event, Vector2 delta) {

  _cameraParams.longitude -= 3.0f * delta.x();
  _cameraParams.latitude += 3.0f * delta.y();

  _cameraParams.latitude = std::clamp(_cameraParams.latitude, -pi / 2, pi / 2);
  _cameraObject->setTransformation(_cameraParams.getCameraTransformation());
}

void MyApplication::mouseZoom(MouseMoveEvent const &event, Vector2 delta) {
  _cameraParams.targetDistance *= 1.0 - 2.0 * delta.y();
  _cameraObject->setTransformation(_cameraParams.getCameraTransformation());
}

void MyApplication::mousePan(MouseMoveEvent const &event, Vector2 delta) {

  Vector2i pmp = _previousMousePosition;
  Vector2i mp  = event.position();

  auto point_from_camera = [&](Vector2i screen_pos, float dist) -> Vector3 {
    auto[dir, cam_pos] = cameraRayCast(screen_pos);
    return cam_pos + dist * dir;
  };

  float dist = _cameraParams.targetDistance;
  _cameraParams.target +=
      point_from_camera(pmp, dist) - point_from_camera(mp, dist);
  _cameraObject->setTransformation(_cameraParams.getCameraTransformation());
}

std::tuple<Vector3, Vector3>
MyApplication::cameraRayCast(Vector2i pixel) const {

  Vector2 screenPos =
      2.0f * (Vector2{pixel} / Vector2{defaultFramebuffer.viewport().size()} -
              Vector2{.5f, 0.5f});
  screenPos[1] *= -1.f;

  return cameraRayCast(screenPos);
}

std::tuple<Vector3, Vector3>
MyApplication::cameraRayCast(Vector2 screenPos) const {
  Matrix4 camTrans = _cameraObject->transformation();
  Matrix4 camProj  = _camera->projectionMatrix();
  Matrix4 trans    = camTrans * camProj.inverted();

  Vector3 point = Vector3{screenPos[0], screenPos[1], 0} +
                  camProj.transformPoint(Vector3{0, 0, -1});

  point = trans.transformPoint(point);

  Vector3 camPos =
      _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  Vector3 dir = (point - camPos).normalized();
  return {dir, camPos};
}

MAGNUM_APPLICATION_MAIN(MyApplication)
