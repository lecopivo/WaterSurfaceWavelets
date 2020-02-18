
#include "SceneBase3D.h"

namespace Magnum {

SceneBase3D::SceneBase3D(const Arguments &arguments)
    : Platform::Application{arguments,
                            Configuration{}
                                .setTitle("SceneBase3D")
                                .setWindowFlags(Sdl2Application::Configuration::
                                                    WindowFlag::Resizable)
                                .setSize({1200, 800}),
                            GLConfiguration{}.setSampleCount(16)},
      _camera(_scene) {

  /* Configure OpenGL state */
  Renderer::enable(Renderer::Feature::DepthTest);
  Renderer::enable(Renderer::Feature::FaceCulling);
  Renderer::enable(Renderer::Feature::Multisampling);

  Renderer::setBlendEquation(Renderer::BlendEquation::Add, Renderer::BlendEquation::Add);
  Renderer::setBlendFunction(Renderer::BlendFunction::SourceAlpha, Renderer::BlendFunction::OneMinusSourceAlpha);

  /* Configure camera */
  viewportEvent(GL::defaultFramebuffer.viewport().size()); // set up camera
  
  _gui = ImGuiIntegration::Context(Vector2{windowSize()}/dpiScaling(), windowSize(), framebufferSize());
}

void SceneBase3D::drawEvent() {
  GL::defaultFramebuffer.clear(FramebufferClear::Color | FramebufferClear::Depth);

  update();

  _camera.draw(_drawables);

  Renderer::enable(Renderer::Feature::Blending);
  Renderer::enable(Renderer::Feature::ScissorTest);
  Renderer::disable(Renderer::Feature::FaceCulling);
  Renderer::disable(Renderer::Feature::DepthTest);

  drawGui();

  Renderer::enable(Renderer::Feature::DepthTest);
  Renderer::enable(Renderer::Feature::FaceCulling);
  Renderer::disable(Renderer::Feature::ScissorTest);
  Renderer::disable(Renderer::Feature::Blending);

  swapBuffers();
  redraw();
}

void SceneBase3D::update() {
  std::cout << "You should override update function!" << std::endl;
}

void SceneBase3D::drawGui() {
  _gui.newFrame();

  ImGui::Text("You should override drawGui function!");

  _gui.drawFrame();
}

void SceneBase3D::showNavigationHelp() const {

  if (ImGui::CollapsingHeader("Navigation Help")) {
    ImGui::Text("Rotate:   Alt + LMB");
    ImGui::Text("Zoom:     Alt + RMB");
    ImGui::Text("Pan:      Alt + MMB");
    ImGui::Text("Recenter: f");
  }
}

void SceneBase3D::viewportEvent(const Vector2i &size) {
  GL::defaultFramebuffer.setViewport({{}, size});

  _camera.setProjectionMatrix(Matrix4::perspectiveProjection(
      Deg{60.0}, Vector2{size}.aspectRatio(), 0.1f, 1000.0f));
}

void SceneBase3D::keyPressEvent(KeyEvent &event) {
  if (_gui.handleKeyPressEvent(event)) {
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

void SceneBase3D::keyReleaseEvent(KeyEvent &event) {
  if (_gui.handleKeyReleaseEvent(event)) {
    redraw();
    return;
  }
}

void SceneBase3D::mousePressEvent(MouseEvent &event) {
  if (_gui.handleMousePressEvent(event)) {
    std::cout.flush();
    redraw();
    return;
  }

  if (event.button() == MouseEvent::Button::Left) {
    _previousMousePosition = event.position();
    event.setAccepted();
  }
}

void SceneBase3D::mouseReleaseEvent(MouseEvent &event) {
  if (_gui.handleMouseReleaseEvent(event)) {
    redraw();
    return;
  }

  event.setAccepted();
  redraw();
}

void SceneBase3D::mouseMoveEvent(MouseMoveEvent &event) {
  if (_gui.handleMouseMoveEvent(event)) {
    redraw();
    return;
  }

  Vector2i mousePosNew = event.position();
  Vector2i mousePosOld = _previousMousePosition;
  Vector2i screenSize  = GL::defaultFramebuffer.viewport().size();

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

void SceneBase3D::mouseScrollEvent(MouseScrollEvent &event) {
  if (_gui.handleMouseScrollEvent(event)) {
    redraw();
    return;
  }
}

void SceneBase3D::textInputEvent(TextInputEvent &event) {
  if (_gui.handleTextInputEvent(event)) {
    redraw();
    return;
  }
}

} // namespace Magnum
