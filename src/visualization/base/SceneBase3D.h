#pragma once


#include <Magnum/GL/Buffer.h>
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
#include <Magnum/GL/Renderer.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/GL/Shader.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Trade/MeshData3D.h>

// #include <MagnumImGui.h>
#include <Magnum/ImGuiIntegration/Context.hpp>
#include <imgui.h>

#include "Camera.h"

#include <iostream>

namespace Magnum {

class SceneBase3D : public Platform::Application {
public:
  using Object3D = SceneGraph::Object<SceneGraph::MatrixTransformation3D>;
  using Scene3D  = SceneGraph::Scene<SceneGraph::MatrixTransformation3D>;

public:
  explicit SceneBase3D(const Arguments &arguments);

  virtual ~SceneBase3D(){};

public:
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

  // MagnumImGui _gui;
  ImGuiIntegration::Context _gui{NoCreate};
};

} // namespace Magnum
