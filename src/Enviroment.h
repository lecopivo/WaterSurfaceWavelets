#pragma once

#include "Global.h"

namespace WaterWavelets {

class Environment {
public:
  bool inDomain(Vec2 pos) const;
  Real levelset(Vec2  pos) const;
  Vec2 levelsetGrad(Vec2 pos) const;

public:
  Vec2 islandCenter = {0.0, 0.0};
  Real islandRadius = 10.0;
};

} // namespace WaterWavelets

