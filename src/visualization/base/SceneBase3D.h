#pragma once

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

#include "Camera.h"
#include "SimulationLoop.h"

#include <iostream>

namespace Magnum {

template <typename Simulation>
class SceneBase3D : public Platform::Application {
public:
  using Object3D = SceneGraph::Object<SceneGraph::MatrixTransformation3D>;
  using Scene3D  = SceneGraph::Scene<SceneGraph::MatrixTransformation3D>;

public:
  explicit SceneBase3D(const Arguments &arguments);

  virtual ~SceneBase3D(){};

protected:
  void drawEvent() final;

  virtual void update();
  virtual void drawGui();

  void showNavigationHelp() const;

  void viewportEvent(const Vector2i &size) final;

  void keyPressEvent(KeyEvent &event) final;
  void keyReleaseEvent(KeyEvent &event) final;
  void mousePressEvent(MouseEvent &event) final;
  void mouseReleaseEvent(MouseEvent &event) final;
  void mouseMoveEvent(MouseMoveEvent &event) final;
  void mouseScrollEvent(MouseScrollEvent &event) final;
  void textInputEvent(TextInputEvent &event) final;

  Scene3D                     _scene;
  SceneGraph::DrawableGroup3D _drawables;

  Vector2i                  _previousMousePosition;
  Camera<Object3D, Scene3D> _camera;

  MagnumImGui _gui;

  Simulation                 _sim;
  SimulationLoop<Simulation> _loop;
  std::mutex                 _updateMutex;
};

template <typename Simulation>
SceneBase3D<Simulation>::SceneBase3D(const Arguments &arguments)
    : Platform::Application{arguments,
                            Configuration{}
                                .setTitle("SceneBase3D")
                                .setWindowFlags(Sdl2Application::Configuration::
                                                    WindowFlag::Resizable)
                                .setSize({1200, 800})},
      _camera(_scene) {

  /* Configure OpenGL state */
  Renderer::enable(Renderer::Feature::DepthTest);
  Renderer::enable(Renderer::Feature::FaceCulling);

  /* Configure camera */
  viewportEvent(defaultFramebuffer.viewport().size()); // set up camera

  /* Start simualation */
  _loop.start(_sim);

  /* Set up loading of drawing state */
  _loop.addTaskPerLoop([this](Simulation &sim) {
    std::lock_guard<std::mutex> guard(_updateMutex);
    update();
  });
}

template <typename Simulation> void SceneBase3D<Simulation>::drawEvent() {
  defaultFramebuffer.clear(FramebufferClear::Color | FramebufferClear::Depth);

  // Not allow updating from simulation
  std::lock_guard<std::mutex> guard(_updateMutex);

  _camera.draw(_drawables);

  drawGui();

  swapBuffers();
  redraw();
}

template <typename Simulation> void SceneBase3D<Simulation>::update() {
  std::cout << "You should override update function!" << std::endl;
}

template <typename Simulation> void SceneBase3D<Simulation>::drawGui() {
  _gui.newFrame(windowSize(), defaultFramebuffer.viewport().size());

  ImGui::Text("You should override drawGui function!");

  _gui.drawFrame();
}

template <typename Simulation>
void SceneBase3D<Simulation>::showNavigationHelp() const {

  if (ImGui::CollapsingHeader("Navigation Help")) {
    ImGui::Text("Rotate:   Alt + LMB");
    ImGui::Text("Zoom:     Alt + RMB");
    ImGui::Text("Pan:      Alt + MMB");
    ImGui::Text("Recenter: f");
  }
}

template <typename Simulation>
void SceneBase3D<Simulation>::viewportEvent(const Vector2i &size) {
  defaultFramebuffer.setViewport({{}, size});

  _camera.setProjectionMatrix(Matrix4::perspectiveProjection(
      Deg{60.0}, Vector2{size}.aspectRatio(), 0.1f, 1000.0f));
}

template <typename Simulation>
void SceneBase3D<Simulation>::keyPressEvent(KeyEvent &event) {
  if (_gui.keyPressEvent(event)) {
    redraw();
    return;
  }

  if (event.key() == KeyEvent::Key::Esc) {
    exit();
  }

  if (event.key() == KeyEvent::Key::F) {
    _camera.centerToOrigin();
  }

  redraw();
}

template <typename Simulation>
void SceneBase3D<Simulation>::keyReleaseEvent(KeyEvent &event) {
  if (_gui.keyReleaseEvent(event)) {
    redraw();
    return;
  }
}

template <typename Simulation>
void SceneBase3D<Simulation>::mousePressEvent(MouseEvent &event) {
  if (_gui.mousePressEvent(event)) {
    std::cout.flush();
    redraw();
    return;
  }

  if (event.button() == MouseEvent::Button::Left) {
    _previousMousePosition = event.position();
    event.setAccepted();
  }
}

template <typename Simulation>
void SceneBase3D<Simulation>::mouseReleaseEvent(MouseEvent &event) {
  if (_gui.mouseReleaseEvent(event)) {
    redraw();
    return;
  }

  event.setAccepted();
  redraw();
}

template <typename Simulation>
void SceneBase3D<Simulation>::mouseMoveEvent(MouseMoveEvent &event) {
  if (_gui.mouseMoveEvent(event)) {
    redraw();
    return;
  }

  Vector2i mousePosNew = event.position();
  Vector2i mousePosOld = _previousMousePosition;
  Vector2i screenSize  = defaultFramebuffer.viewport().size();

  if ((event.modifiers() & MouseMoveEvent::Modifier::Alt) &&
      (event.buttons() & MouseMoveEvent::Button::Left))
    _camera.rotate(mousePosNew, mousePosOld, screenSize);

  if (event.modifiers() & MouseMoveEvent::Modifier::Alt &&
      event.buttons() & MouseMoveEvent::Button::Right)
    _camera.zoom(mousePosNew, mousePosOld, screenSize);

  if (event.modifiers() & MouseMoveEvent::Modifier::Alt &&
      event.buttons() & MouseMoveEvent::Button::Middle)
    _camera.pan(mousePosNew, mousePosOld, screenSize);

  _previousMousePosition = event.position();
  event.setAccepted();
  redraw();
}

template <typename Simulation>
void SceneBase3D<Simulation>::mouseScrollEvent(MouseScrollEvent &event) {
  if (_gui.mouseScrollEvent(event)) {
    redraw();
    return;
  }
}

template <typename Simulation>
void SceneBase3D<Simulation>::textInputEvent(TextInputEvent &event) {
  if (_gui.textInputEvent(event)) {
    redraw();
    return;
  }
}

} // namespace Magnum
