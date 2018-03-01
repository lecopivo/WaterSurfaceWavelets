#include "WaveGrid.hpp"

#include <algorithm>

constexpr int pos_modulo(int n, int d) { return (n % d + d) % d; }

WaveGrid::WaveGrid(Settings s) {

  m_amplitude.resize(s.n_x, s.n_y, s.n_theta, s.n_k);
  m_newAmplitude.resize(s.n_x, s.n_y, s.n_theta, s.n_k);

  Real kmin = 8.0;
  Real kmax = 16.0;

  m_xmin = {s.xmin, s.ymin, 0.0, kmin};
  m_xmax = {s.xmax, s.ymax, 2.0 * M_PI, kmax};

  for (int i = 0; i < 4; i++) {
    m_dx[i] = (m_xmax[i] - m_xmin[i]) / m_amplitude.dimension(i);
    m_idx[i] = 1.0 / m_dx[i];
  }

  // boundary info
  islandCenter << 0.5 * (m_xmin[0] + m_xmax[0]), 0.5 * (m_xmin[1] + m_xmax[1]);
  islandRadius = 0.3 * std::min(m_xmax[0] - m_xmin[0], m_xmax[1] - m_xmin[1]);
}

void WaveGrid::timeStep(const double dt, const double t) {
  advectionStep(dt);
  // diffusionStep(dt);
  precompute1DBuffers(t);
}

Real WaveGrid::cflTimeStep() const{
  return std::min(m_dx[0],m_dx[1])/groupSpeed(waveNumber(0));
}

auto WaveGrid::interpolateAmplitude(Grid const &grid) const {

  // Make acess to grid possible even outside of the domain
  auto extended_grid = [&grid, this](int a_x, int a_y, int b_theta, int c_k) {
    // wrap arround for angle
    b_theta = pos_modulo(b_theta, grid.dimension(2));

    // return zero for wavenumber outside of a domain
    if (c_k < 0 || c_k > grid.dimension(3)) {
      return 0.0;
    }

    // return a default value for points outside of the simulation box
    if (a_x < 0 || a_x >= grid.dimension(0) || a_y < 0 ||
        a_y >= grid.dimension(1)) {
      if (b_theta == 0)
        return defaultAmplitude(b_theta, c_k);
      else
        return 0.0;
    }

    // if the point is in the domain the return the actual value of the grid
    return grid(a_x, a_y, b_theta, c_k);
  };

  auto domain = [this](int a_x, int a_y, int b_theta, int c_k) -> bool {
    return inDomain(nodePosition(a_x, a_y));
  };

  auto interpolation =
      InterpolationDimWise(LinearInterpolation, LinearInterpolation,
                           LinearInterpolation, ConstantInterpolation);

  auto igrid = DomainInterpolation(interpolation, domain)(extended_grid);

  return [igrid, this](Vec4 pos4) mutable {
    Vec4 gridPos4 = posToGrid(pos4);
    return igrid(gridPos4(0), gridPos4(1), gridPos4(2), gridPos4(3));
  };
}

void WaveGrid::advectionStep(const double dt) {

  auto amplitude = interpolateAmplitude(m_amplitude);

  for (int a_x = 0; a_x < gridDim(0); a_x++) {
    for (int a_y = 0; a_y < gridDim(1); a_y++) {

      Vec2 pos = nodePosition(a_x, a_y);

      if (inDomain(pos)) {

        for (int b_theta = 0; b_theta < gridDim(2); b_theta++) {
          for (int c_k = 0; c_k < gridDim(3); c_k++) {

            Real theta = waveVectorAngle(b_theta);
            Real k = waveNumber(c_k);
            Vec2 cg = groupVelocity(theta, k);

	    Vec2 newPos = pos - dt*cg;
            Vec4 newPos4;
            newPos4 << newPos(0), newPos(1), theta, k;

            newPos4 = boundaryReflection(newPos4);

            m_newAmplitude(a_x, a_y, b_theta, c_k) = amplitude(newPos4);
          }
        }
      }
    }
  }

  std::swap(m_newAmplitude, m_amplitude);
}

/*
 *
 * \param gridPos position in 4D grid. In grid space coordinates!
 * \return position in 4D grid which corresponds to a reflected position
 * around boundary if necessary
 */
Vec4 WaveGrid::boundaryReflection(const Vec4 pos4) const {
  Vec2 pos = pos4.segment<2>(0);
  Real ls = levelset(pos);
  if (ls >= 0) // no reflection is needed if point is in the domain
    return pos4;

  // Boundary normal is approximatex by the levelset gradient
  Vec2 n = levelsetGrad(pos);

  pos = pos - 2.0*ls * n;
  assert(inDomain(pos));

  // reflect wave vector around boundary normal
  Real theta = pos4(2); // This is in grid space!
  Real normalAngle = atan2(n(1), n(0));
  Real reflected_theta = -theta + normalAngle + M_PI;

  Vec4 newPos4;
  newPos4 << pos(0), pos(1), reflected_theta, pos4(3);

  return newPos4;
}

bool WaveGrid::inDomain(const Vec2 pos) const { return levelset(pos) >= 0; }

Real WaveGrid::levelset(const Vec2 pos) const {
  return (pos - islandCenter).norm() - islandRadius;
}

Vec2 WaveGrid::levelsetGrad(const Vec2 pos) const {
  Vec2 x = pos - islandCenter;
  return x.normalized();
}

// void diffusionStep(const double dt) {
//   for (int a_x = 0; a_x < gridDim(0); a_x++) {
//     for (int a_y = 0; a_x < gridDim(1); a_y++) {
//       for (int b_theta = 0; a_x < gridDim(2); b_theta++) {
//         for (int c_k = 0; a_x < gridDim(3); c_k++) {

//           // this is
//           double gamma = 0.025 * groupSpeed(waveNumber(c_k)) * dt / m_dx(0);
//           double delta =
//               1e-5 * dt * pow(m_dx(3), 2) * dispersion(waveNumber(c_k));

//           m_newAmplitude({a_x, a_y, b_theta, c_k}) =
//               (1 - gamma_k - gamma_alpha) *
//                   m_amplitude({a_x, a_y, b_theta, c_k}) +
//               0.5 * gamma_k *
//                   (m_amplitude({a_x, a_y, b_theta + 1, c_k}) +
//                    m_amplitude({a_x, a_y, b_theta - 1, c_k})) +
//               0.5 * gamma_alpha *
//                   (m_amplitude({a_x, a_y, b_theta, c_k + 1}) +
//                    m_amplitude({a_x, a_y, b_theta, c_k + 1}));
//         }
//       }
//     }
//   }
// }

void WaveGrid::precompute1DBuffers(const double dt) {
  for (int c_k = 0; c_k < gridDim(3); c_k++) {
    
  }
}

Real WaveGrid::waterHeight(const Vec2 x) { return 0.0; }

Real WaveGrid::idxToPos(const int idx, const int dim) const {
  return m_xmin[dim] + (idx + 0.5) * m_dx[dim];
}

Vec4 WaveGrid::idxToPos(const Idx idx) const {
  Vec4 pos4;
  pos4 << idxToPos(idx[0], 0), idxToPos(idx[1], 1), idxToPos(idx[2], 2),
      idxToPos(idx[3], 3);
  return pos4;
}

Real WaveGrid::posToGrid(const Real pos, const int dim) const {
  return (pos - m_xmin[dim]) * m_idx[dim] - 0.5;
}

Vec4 WaveGrid::posToGrid(const Vec4 pos4) const {
  Vec4 grid4;
  grid4 << posToGrid(pos4(0), 0), posToGrid(pos4(1), 1), posToGrid(pos4(2), 2),
      posToGrid(pos4(3), 3);
  return grid4;
}

int WaveGrid::posToIdx(const Real pos, const int dim) const {
  return round(posToGrid(pos, dim));
}

Idx WaveGrid::posToIdx(const Vec4 pos4) const {
  return {posToIdx(pos4(0), 0), posToIdx(pos4(1), 1), posToIdx(pos4(2), 2),
          posToIdx(pos4(3), 3)};
}

Vec2 WaveGrid::nodePosition(const int a_x, const int a_y) const {
  Vec2 pos;
  pos << idxToPos(a_x, 0), idxToPos(a_y, 1);
  return pos;
}

Real WaveGrid::waveVectorAngle(const int angleIdx) const {
  return idxToPos(angleIdx, 2);
}

Real WaveGrid::waveNumber(const int waveNumberIdx) const {
  return idxToPos(waveNumberIdx, 3);
}

Vec2 WaveGrid::waveDirection(const int angleIdx) const {
  Real angle = waveVectorAngle(angleIdx);
  Vec2 d;
  d << cos(angle), sin(angle);
  return d;
}

Vec2 WaveGrid::waveVector(const int b_theta, const int c_k) const {
  return waveDirection(b_theta) * waveNumber(c_k);
}

Real WaveGrid::dispersionRelation(const Real k) const {
  const Real g = 9.81;
  return sqrt(k * g);
}

Real WaveGrid::groupSpeed(const Real k) const {
  const Real g = 9.81;
  return 0.5 * sqrt(g / k);
}

Vec2 WaveGrid::groupVelocity(const Real theta, const Real k) const {
  Real cg = groupSpeed(k);
  Vec2 vel;
  vel << cg * cos(theta), cg * sin(theta);
  return vel;
}

Real WaveGrid::defaultAmplitude(const int b_theta, const int c_k) const {
  if(c_k==0)
    return 1.0;
  return 0.0;
}

int WaveGrid::gridDim(const int dim) const {
  return m_amplitude.dimension(dim);
}
