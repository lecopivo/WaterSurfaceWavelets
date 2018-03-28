#include "WaveGrid.h"
#include "math/Math.h"

#include <algorithm>

namespace WaterWavelets {

constexpr int pos_modulo(int n, int d) { return (n % d + d) % d; }

constexpr Real tau = 6.28318530718; // https://tauday.com/tau-manifesto

  WaveGrid::WaveGrid(Settings s) : m_spectrum(10),m_enviroment(s.size) {

  m_amplitude.resize(s.n_x, s.n_x, s.n_theta, s.n_zeta);
  m_newAmplitude.resize(s.n_x, s.n_x, s.n_theta, s.n_zeta);

  Real zeta_min = m_spectrum.minZeta();
  Real zeta_max = m_spectrum.maxZeta();

  m_xmin = {-s.size, -s.size, 0.0, zeta_min};
  m_xmax = {s.size, s.size, tau, zeta_max};

  for (int i = 0; i < 4; i++) {
    m_dx[i]  = (m_xmax[i] - m_xmin[i]) / m_amplitude.dimension(i);
    m_idx[i] = 1.0 / m_dx[i];
  }

  m_time = s.initial_time;

  m_profileBuffers.resize(s.n_zeta);
  precomputeGroupSpeeds();
}

void WaveGrid::timeStep(const Real dt, bool fullUpdate) {
  if (fullUpdate) {
    advectionStep(dt);
    diffusionStep(dt);
  }
  precomputeProfileBuffers();
  m_time += dt;
}

Real WaveGrid::cflTimeStep() const {
  return std::min(m_dx[X], m_dx[Y]) / groupSpeed(gridDim(Zeta) - 1);
}

std::pair<Vec3, Vec3> WaveGrid::waterSurface(Vec2 pos) const {

  Vec3 surface = {0, 0, 0};
  Vec3 tx      = {0, 0, 0};
  Vec3 ty      = {0, 0, 0};

  for (int izeta = 0; izeta < gridDim(Zeta); izeta++) {
    Real  zeta    = idxToPos(izeta, Zeta);
    auto &profile = m_profileBuffers[izeta];

    int  DIR_NUM = gridDim(Theta);
    int  N       = 4 * DIR_NUM;
    Real da      = 1.0 / N;
    Real dx      = DIR_NUM * tau / N;
    for (Real a = 0; a < 1; a += da) {

      Real angle  = a * tau;
      Vec2 kdir   = Vec2{cosf(angle), sinf(angle)};
      Real kdir_x = kdir * pos;

      Vec4 wave_data =
          dx * amplitude({pos[X], pos[Y], angle, zeta}) * profile(kdir_x);

      surface +=
          Vec3{kdir[0] * wave_data[0], kdir[1] * wave_data[0], wave_data[1]};

      tx += kdir[0] * Vec3{wave_data[2], 0, wave_data[3]};
      ty += kdir[1] * Vec3{0, wave_data[2], wave_data[3]};
    }
  }

  Vec3 normal = normalized(cross(tx, ty));

  return {surface, normal};
}

auto WaveGrid::extendedGrid() const {
  return [this](int ix, int iy, int itheta, int izeta) {
    // wrap arround for angle
    itheta = pos_modulo(itheta, gridDim(Theta));

    // return zero for wavenumber outside of a domain
    if (izeta < 0 || izeta >= gridDim(Zeta)) {
      return 0.0f;
    }

    // return a default value for points outside of the simulation box
    if (ix < 0 || ix >= gridDim(X) || iy < 0 || iy >= gridDim(Y)) {
      return defaultAmplitude(itheta, izeta);
    }

    // if the point is in the domain the return the actual value of the grid
    return m_amplitude(ix, iy, itheta, izeta);
  };
}

auto WaveGrid::interpolatedAmplitude() const {

  auto extended_grid = extendedGrid();

  // This function indicated which grid points are in domain and which are not
  auto domain = [this](int ix, int iy, int itheta, int izeta) -> bool {
    return m_enviroment.inDomain(nodePosition(ix, iy));
  };

  auto interpolation = InterpolationDimWise(
      // CubicInterpolation, CubicInterpolation,
      LinearInterpolation, LinearInterpolation, LinearInterpolation,
      ConstantInterpolation);

  auto interpolated_grid =
      DomainInterpolation(interpolation, domain)(extended_grid);

  return [interpolated_grid, this](Vec4 pos4) mutable {
    Vec4 ipos4 = posToGrid(pos4);
    return interpolated_grid(ipos4[X], ipos4[Y], ipos4[Theta], ipos4[Zeta]);
  };
}

Real WaveGrid::amplitude(Vec4 pos4) const {
  return interpolatedAmplitude()(pos4);
}

Real WaveGrid::gridValue(Idx idx4) const {
  return extendedGrid()(idx4[X], idx4[Y], idx4[Theta], idx4[Zeta]);
}

std::vector<Vec4> WaveGrid::trajectory(Vec4 pos4, Real length) const {

  std::vector<Vec4> trajectory;
  Real              dist = 0;

  for (Real dist = 0; dist <= length;) {

    trajectory.push_back(pos4);

    Vec2 vel = groupVelocity(pos4);
    Real dt  = dx(X) / norm(vel);

    pos4[X] += dt * vel[X];
    pos4[Y] += dt * vel[Y];

    pos4 = boundaryReflection(pos4);

    dist += dt * norm(vel);
  }
  trajectory.push_back(pos4);
  return trajectory;
}

void WaveGrid::addPointDisturbance(const Vec2 pos, const Real val) {
  // Find the closest point on the grid to the point `pos`
  int ix = posToIdx(pos[X], X);
  int iy = posToIdx(pos[Y], Y);
  if (ix >= 0 && ix < gridDim(X) && iy >= 0 && iy < gridDim(Y)) {

    for (int itheta = 0; itheta < gridDim(Theta); itheta++) {
      m_amplitude(ix, iy, itheta, 0) += val;
    }
  }
}

void WaveGrid::advectionStep(const Real dt) {

  auto amplitude = interpolatedAmplitude();

#pragma omp parallel for collapse(2)
  for (int ix = 0; ix < gridDim(X); ix++) {
    for (int iy = 0; iy < gridDim(Y); iy++) {

      Vec2 pos = nodePosition(ix, iy);

      // update only points in the domain
      if (m_enviroment.inDomain(pos)) {

        for (int itheta = 0; itheta < gridDim(Theta); itheta++) {
          for (int izeta = 0; izeta < gridDim(Zeta); izeta++) {

            Vec4 pos4 = idxToPos({ix, iy, itheta, izeta});
            Vec2 vel  = groupVelocity(pos4);

            // Tracing back in Semi-Lagrangian
            Vec4 trace_back_pos4 = pos4;
            trace_back_pos4[X] -= dt * vel[X];
            trace_back_pos4[Y] -= dt * vel[Y];

            // Take care of boundaries
            trace_back_pos4 = boundaryReflection(trace_back_pos4);

            m_newAmplitude(ix, iy, itheta, izeta) = amplitude(trace_back_pos4);
          }
        }
      }
    }
  }

  std::swap(m_newAmplitude, m_amplitude);
}

Vec4 WaveGrid::boundaryReflection(const Vec4 pos4) const {
  Vec2 pos = Vec2{pos4[X], pos4[Y]};
  Real ls  = m_enviroment.levelset(pos);
  if (ls >= 0) // no reflection is needed if point is in the domain
    return pos4;

  // Boundary normal is approximatex by the levelset gradient
  Vec2 n = m_enviroment.levelsetGrad(pos);

  Real theta = pos4[Theta];
  Vec2 kdir  = Vec2{cosf(theta), sinf(theta)};

  // Reflect point and wave-vector direction around boundary
  // Here we rely that `ls` is equal to the signed distance from the boundary
  pos  = pos - 2.0 * ls * n;
  kdir = kdir - 2.0 * (kdir * n) * n;

  Real reflected_theta = atan2(kdir[Y], kdir[X]);

  // We are assuming that after one reflection you are back in the domain. This
  // assumption is valid if you boundary is not so crazy.
  // This assert tests this assumption.
  assert(m_enviroment.inDomain(pos));

  return Vec4{pos[X], pos[Y], reflected_theta, pos4[Zeta]};
}

void WaveGrid::diffusionStep(const Real dt) {

  auto grid = extendedGrid();

#pragma omp parallel for collapse(2)
  for (int ix = 0; ix < gridDim(X); ix++) {
    for (int iy = 0; iy < gridDim(Y); iy++) {

      float ls = m_enviroment.levelset(nodePosition(ix, iy));

      for (int itheta = 0; itheta < gridDim(Theta); itheta++) {
        for (int izeta = 0; izeta < gridDim(Zeta); izeta++) {

          Vec4 pos4  = idxToPos({ix, iy, itheta, izeta});
          Real gamma = 0.025 * groupSpeed(izeta) * dt * m_idx[X];

          // do diffusion only if you are 2 grid nodes away from boudnary
          if (ls >= 4 * dx(X)) {
            m_newAmplitude(ix, iy, itheta, izeta) =
                (1 - gamma) * grid(ix, iy, itheta, izeta) +
                gamma * 0.5 *
                    (grid(ix, iy, itheta + 1, izeta) +
                     grid(ix, iy, itheta - 1, izeta));
          } else {
            m_newAmplitude(ix, iy, itheta, izeta) = grid(ix, iy, itheta, izeta);
          }
          // auto dispersion = [](int i) { return 1.0; };
          // Real delta =
          //     1e-5 * dt * pow(m_dx[3], 2) * dispersion(waveNumber(izeta));
          // 0.5 * delta *
          //     (m_amplitude(ix, iy, itheta, izeta + 1) +
          //      m_amplitude(ix, iy, itheta, izeta + 1));
        }
      }
    }
  }
  std::swap(m_newAmplitude, m_amplitude);
}

void WaveGrid::precomputeProfileBuffers() {

  for (int izeta = 0; izeta < gridDim(Zeta); izeta++) {

    Real zeta_min = idxToPos(izeta, Zeta) - 0.5 * dx(Zeta);
    Real zeta_max = idxToPos(izeta, Zeta) + 0.5 * dx(Zeta);

    // define spectrum

    m_profileBuffers[izeta].precompute(m_spectrum, m_time, zeta_min, zeta_max);
  }
}

void WaveGrid::precomputeGroupSpeeds() {
  m_groupSpeeds.resize(gridDim(Zeta));
  for (int izeta = 0; izeta < gridDim(Zeta); izeta++) {

    Real zeta_min = idxToPos(izeta, Zeta) - 0.5 * dx(Zeta);
    Real zeta_max = idxToPos(izeta, Zeta) + 0.5 * dx(Zeta);

    auto result = integrate(100, zeta_min, zeta_max, [&](Real zeta) -> Vec2 {
      Real waveLength = pow(2, zeta);
      Real waveNumber = tau / waveLength;
      Real cg         = 0.5 * sqrt(9.81 / waveNumber);
      Real density    = m_spectrum(zeta);
      return {cg * density, density};
    });

    m_groupSpeeds[izeta] =
        3 /*the 3 should not be here !!!*/ * result[0] / result[1];
  }
}

Real WaveGrid::idxToPos(const int idx, const int dim) const {
  return m_xmin[dim] + (idx + 0.5) * m_dx[dim];
}

Vec4 WaveGrid::idxToPos(const Idx idx) const {
  return Vec4{idxToPos(idx[X], X), idxToPos(idx[Y], Y),
              idxToPos(idx[Theta], Theta), idxToPos(idx[Zeta], Zeta)};
}

Real WaveGrid::posToGrid(const Real pos, const int dim) const {
  return (pos - m_xmin[dim]) * m_idx[dim] - 0.5;
}

Vec4 WaveGrid::posToGrid(const Vec4 pos4) const {
  return Vec4{posToGrid(pos4[X], X), posToGrid(pos4[Y], Y),
              posToGrid(pos4[Theta], Theta), posToGrid(pos4[Zeta], Zeta)};
}

int WaveGrid::posToIdx(const Real pos, const int dim) const {
  return round(posToGrid(pos, dim));
}

WaveGrid::Idx WaveGrid::posToIdx(const Vec4 pos4) const {
  return Idx{posToIdx(pos4[X], X), posToIdx(pos4[Y], Y),
             posToIdx(pos4[Theta], Theta), posToIdx(pos4[Zeta], Zeta)};
}

Vec2 WaveGrid::nodePosition(int ix, int iy) const {
  return Vec2{idxToPos(ix, 0), idxToPos(iy, 1)};
}

Real WaveGrid::waveLength(int izeta) const {
  Real zeta = idxToPos(izeta, Zeta);
  return pow(2, zeta);
}

Real WaveGrid::waveNumber(int izeta) const { return tau / waveLength(izeta); }

Real WaveGrid::dispersionRelation(Real knum) const {
  const Real g = 9.81;
  return sqrt(knum * g);
}

Real WaveGrid::dispersionRelation(Vec4 pos4) const {
  Real knum = waveNumber(pos4[Zeta]);
  return dispersionRelation(knum);
}

// Real WaveGrid::groupSpeed(Vec4 pos4) const {
//   Real       knum = waveNumber(pos4[Zeta]);
//   const Real g    = 9.81;
//   return 0.5 * sqrt(g / knum);
// }

Real WaveGrid::groupSpeed(int izeta) const { return m_groupSpeeds[izeta]; }

Vec2 WaveGrid::groupVelocity(Vec4 pos4) const {
  int  izeta = posToIdx(pos4[Zeta], Zeta);
  Real cg    = groupSpeed(izeta);
  Real theta = pos4[Theta];
  return cg * Vec2{cosf(theta), sinf(theta)};
}

Real WaveGrid::defaultAmplitude(const int itheta, const int izeta) const {
  if (itheta == 5 * gridDim(Theta) / 16)
    return 0.1;
  return 0.0;
}

int WaveGrid::gridDim(const int dim) const {
  return m_amplitude.dimension(dim);
}

Real WaveGrid::dx(int dim) const { return m_dx[dim]; }

} // namespace WaterWavelets
