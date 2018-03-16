#include "Enviroment.h"

#include "math/ArrayAlgebra.h"

namespace WaterWavelets {

bool Enviroment::inDomain(Vec2 pos) const { return levelset(pos) >= 0; }

Real Enviroment::levelset(Vec2 pos) const {
  return norm(pos - islandCenter) - islandRadius;
}

Vec2 Enviroment::levelsetGrad(Vec2 pos) const {
  Vec2 x = pos - islandCenter;
  return normalized(x);
}

} // namespace WaterWavelets
