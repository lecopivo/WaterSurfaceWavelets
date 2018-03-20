#pragma once

#include "global.h"

namespace WaterWavelets {

class Enviroment {
public:
  bool inDomain(Vec2 pos) const;
  Real levelset(Vec2  pos) const;
  Vec2 levelsetGrad(Vec2 pos) const;

public:
  Vec2 islandCenter = {0.0, 0.0};
  Real islandRadius = 20.0;
};

} // namespace WaterWavelets
