#pragma once

#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Scene.h>

#include <algorithm>
#include <cmath>
#include <tuple>

namespace Magnum {

template <typename Object3D, typename Scene3D> class Camera {
public:
  Camera(Scene3D &scene) {
    _object  = new Object3D{&scene};
    _feature = new SceneGraph::Camera3D{*_object};

    recomputeTransformation();
  }

  void setProjectionMatrix(Matrix4 matrix) {
    _feature->setProjectionMatrix(matrix);
  }

  void centerToOrigin() {
    target = Vector3{0, 0, 0};
    recomputeTransformation();
  }

  void draw(SceneGraph::DrawableGroup3D &drawables) {
    _feature->draw(drawables);
  }

  void rotate(Vector2i screenPosNew, Vector2i screenPosOld,
              Vector2i screenSize) {

    const Vector2 delta =
        Vector2{screenPosNew - screenPosOld} / Vector2{screenSize};

    rotate(-3.0f * delta.x(), 3.0f * delta.y());
  }

  void rotate(float delta_longitude, float delta_latitude) {
    longitude += delta_longitude;
    latitude =
        std::clamp(1.0 * (latitude + delta_latitude), -M_PI / 2, M_PI / 2);
    recomputeTransformation();
  }

  void zoom(Vector2i screenPosNew, Vector2i screenPosOld, Vector2i screenSize) {
    float mult = 1.0 + 2.0 * (screenPosOld.y() - screenPosNew.y()) /
                           (0.5 * screenSize.y());
    zoom(mult);
  }

  void zoom(float mult) {
    targetDistance *= mult;
    recomputeTransformation();
  }

  void pan(Vector2i screenPosNew, Vector2i screenPosOld, Vector2i screenSize) {

    auto point_from_camera = [&, this](Vector2i screenPos) -> Vector3 {
      auto[dir, cam_pos] = cameraRayCast(screenPos, screenSize);
      return cam_pos + targetDistance * dir;
    };

    Vector3 targetShift =
        point_from_camera(screenPosOld) - point_from_camera(screenPosNew);

    pan(targetShift);
  }

  void pan(Vector3 targetShift) {
    target += targetShift;

    recomputeTransformation();
  }

  std::tuple<Vector3, Vector3> cameraRayCast(Vector2i pixel,
                                             Vector2i screenSize) const {
    Vector2 screenPos =
        2.0f * (Vector2{pixel} / Vector2{screenSize} - Vector2{.5f, 0.5f});
    screenPos[1] *= -1.f;

    return cameraRayCast(screenPos);
  }

  std::tuple<Vector3, Vector3> cameraRayCast(Vector2 screenPos) const {
    Matrix4 camTrans = _object->transformation();
    Matrix4 camProj  = _feature->projectionMatrix();
    Matrix4 trans    = camTrans * camProj.inverted();

    Vector3 point = Vector3{screenPos[0], screenPos[1], 0} +
                    camProj.transformPoint(Vector3{0, 0, -1});

    point = trans.transformPoint(point);

    Vector3 camPos =
        _object->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
    Vector3 dir = (point - camPos).normalized();
    return {dir, camPos};
  }

private:
  void recomputeTransformation() {
    Matrix4 trans = Matrix4::translation(target) *
                    Matrix4::rotationZ(Rad{longitude}) *
                    Matrix4::rotationX(-Rad{latitude}) *
                    Matrix4::translation(Vector3{0.f, -targetDistance, 0.f}) *
                    Matrix4::rotationX(Rad{M_PI / 2});
    _object->setTransformation(trans);
  }

private:
  Object3D *            _object;
  SceneGraph::Camera3D *_feature;

  Vector3 target         = {0.f, 0.f, 0.f};
  float   longitude      = M_PI / 4;
  float   latitude       = M_PI / 4;
  float   targetDistance = 10.0f;
};

} // namespace Magnum
