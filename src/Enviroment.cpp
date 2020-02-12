#include "Enviroment.h"

#include "data/harbor_data.cpp"
#include "data/island_data.cpp"
#include "math/ArrayAlgebra.h"
#include "math/interpolation/Interpolation.h"

#include <cmath>
#include <iostream>

namespace WaterWavelets {

auto &raw_data = harbor_data;

const int N = sqrt(sizeof(raw_data) / sizeof(float));

auto data_grid = [](int i, int j) -> float {
  // outside of the data grid just return some high arbitrary number
  if (i < 0 || i >= N || j < 0 || j >= N)
    return 100;

  return raw_data[j + i * N];
};

auto dx_data_grid = [](int i, int j) -> float {
  if (i < 0 || i >= (N - 1) || j < 0 || j >= N)
    return 0;

  return raw_data[j + (i + 1) * N] - raw_data[j + i * N];
};

auto dy_data_grid = [](int i, int j) -> float {
  if (i < 0 || i >= N || j < 0 || j >= (N - 1))
    return 0;

  return raw_data[(j + 1) + i * N] - raw_data[j + i * N];
};

auto igrid =
    InterpolationDimWise(LinearInterpolation, LinearInterpolation)(data_grid);

auto igrid_dx = InterpolationDimWise(LinearInterpolation,
                                     LinearInterpolation)(dx_data_grid);

auto igrid_dy = InterpolationDimWise(LinearInterpolation,
                                     LinearInterpolation)(dy_data_grid);

auto grid = [](Vec2 pos, float dx) -> float {
  pos *= 1 / dx;
  pos += Vec2{N / 2 - 0.5f, N / 2 - 0.5f};
  return igrid(pos[0], pos[1]) * dx;
};

auto grid_dx = [](Vec2 pos, float dx) -> float {
  pos *= 1 / dx;
  pos += Vec2{N / 2 - 1.0f, N / 2 - 0.5f};
  return igrid_dx(pos[0], pos[1]);
};

auto grid_dy = [](Vec2 pos, float dx) -> float {
  pos *= 1 / dx;
  pos += Vec2{N / 2 - 0.5f, N / 2 - 1.0f};
  return igrid_dy(pos[0], pos[1]);
};

Environment::Environment(float size) : _dx((2 * size) / N) {}

bool Environment::inDomain(Vec2 pos) const { return levelset(pos) >= 0; }

Real Environment::levelset(Vec2 pos) const { return grid(pos, _dx); }

Vec2 Environment::levelsetGrad(Vec2 pos) const {
  Vec2 grad = Vec2{grid_dx(pos, _dx), grid_dy(pos, _dx)};
  return normalized(grad);
}

} // namespace WaterWavelets
