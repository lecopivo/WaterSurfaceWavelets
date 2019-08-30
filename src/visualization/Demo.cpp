#include "base/SceneBase3D.h"
#include "drawables/Primitives3D.h"

#include "../WaveGrid.h"
#include "waterSurface/WaterSurfaceMesh.h"

#include <iostream>

using namespace Magnum;
using namespace WaterWavelets;

// Simulation Settings
auto settings = []() { // this use of lambda is basically emulating Designated
                       // initializers from C++20

  auto s = WaveGrid::Settings{};

  s.size     = 50;
  s.min_zeta = log2(0.03);
  s.max_zeta = log2(10);

  s.n_x     = 100;
  s.n_theta = DIR_NUM;
  s.n_zeta  = 1;

  s.initial_time = 100;

  return s;
}();

// Visualization Settings
int   visGridResolution  = 100;
float amplitudeMult      = 4.0;
bool  update_screen_grid = true;
float logdt              = -0.9;
int   directionToShow    = -1;

class Scene3D : public SceneBase3D {
public:
  explicit Scene3D(const Arguments &arguments)
      : SceneBase3D(arguments), _waveGrid(settings) {

    _water_surface =
        new WaterSurfaceMesh(&_scene, &_drawables, visGridResolution);

    _plane = new Drawables::Plane(&_scene, &_drawables, 100, 100);

    _plane->setVertices([&](int i, Drawables::Plane::VertexData &v) {
      auto &pos = v.position;
      pos *= 0.99 * settings.size;
      pos.z() =
          -10 * tanh(0.1 * _waveGrid.m_enviroment.levelset({pos.x(), pos.y()}));
    });
    _plane->_shader.setDiffuseColor(Color4{0.4, 0.4, 0.4, 1});
    
    _camera.zoom(10);
  }

  void update() override {

    // Load amplitude data from simulation to visualization grid
    if (update_screen_grid) {
      _water_surface->setVertices([&](int i, WaterSurfaceMesh::VertexData &v) {

        int     ix        = i / (visGridResolution + 1);
        int     iy        = i % (visGridResolution + 1);
        Vector2 screenPos = Vector2{(2.0f * ix) / visGridResolution - 1.0f,
                                    (2.0f * iy) / visGridResolution - 1.0f};

        auto[dir, camPos] = _camera.cameraRayCast(screenPos);
        dir               = dir.normalized();
        float t           = -camPos.z() / dir.z();
        t                 = t < 0 ? 1000 : t;
        v.position        = camPos + t * dir;
        v.position.z()    = 0;

        for (int itheta = 0; itheta < DIR_NUM; itheta++) {
          float theta = _waveGrid.idxToPos(itheta, WaveGrid::Theta);
          Vec4  pos4{v.position.x(), v.position.y(), theta,
                    _waveGrid.idxToPos(0, WaveGrid::Zeta)};
          if (directionToShow == -1 || directionToShow == itheta)
            v.amplitude[itheta] = amplitudeMult * _waveGrid.amplitude(pos4);
          else
            v.amplitude[itheta] = 0;
        }
      });
    }

    // Load profile data from simulation to visualization grid
    _water_surface->loadProfile(_waveGrid.m_profileBuffers[0]);

    // Time step of simulation
    _waveGrid.timeStep(_waveGrid.cflTimeStep() * pow(10, logdt),
                       update_screen_grid);

    auto[dir, camPos] = _camera.cameraRayCast(
					      _previousMousePosition, GL::defaultFramebuffer.viewport().size());

    // float t = -camPos.z()/dir.z();
    // Vector3 pos = camPos + t * dir;
    // _waveGrid.addPointDisturbance({pos[0],pos[1]}, 0.1);
  }

  void drawGui() override {
    _gui.newFrame();

    ImGui::Text("Water Surface Wavelts Demo!");

    ImGui::SliderInt("direction", &directionToShow, -1, DIR_NUM - 1);
    ImGui::SliderFloat("amplitude", &amplitudeMult, 0, 20);
    ImGui::SliderFloat("log(dt)", &logdt, -2, 2);
    ImGui::DragFloat("time", &_waveGrid.m_time);
    ImGui::Checkbox("update screen grid", &update_screen_grid);
    if(ImGui::Button("Reset")){
      _waveGrid = WaveGrid(settings);
    }
    static bool triangulation = false;
    if (ImGui::Checkbox("triangulation", &triangulation)) {
      _water_surface->_showTriangulation = triangulation;
    }

    showNavigationHelp();
    _gui.drawFrame();
  }

private:
  WaveGrid _waveGrid;

  Drawables::Plane *_plane;
  WaterSurfaceMesh *_water_surface;
};

MAGNUM_APPLICATION_MAIN(Scene3D)
