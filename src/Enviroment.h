#pragma once

#include "Global.h"

namespace WaterWavelets {

class Environment {
public:
  Environment(float size);

  bool inDomain(Vec2 pos) const;
  Real levelset(Vec2 pos) const;
  Vec2 levelsetGrad(Vec2 pos) const;

public:
  float _dx;
};

} // namespace WaterWavelets
