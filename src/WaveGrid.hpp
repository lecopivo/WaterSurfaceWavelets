#pragma once

#include <array>
#include <iostream>
#include <vector>

#include <Eigen/../unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

#include <Interpolation/DomainInterpolation.h>
#include <Interpolation/Interpolation.h>

constexpr int pos_modulo(int n, int d);

using Real = double;

class BoundaryLevelSet {
public:
  // evaluate Levelset 
  Real operator()(Real x, Real y) { return 0; }

  
};

/*! \typedef Vec4 Location in 4D grid!
 *
 * \brief Location (x,y,theta,zeta) in 4D grid.
 *
 *  The location is determined by four numbers: x,y,theta,zeta
 *  $x \in [-size,size]$ - the first spatial coordinate
 *  $y \in [-size,size]$ - the second spatial coordinate
 *  $theta \in [0,2 \pi)$  - direction of wavevector, theta==0 corresponds to
 * wavevector in +x direction
 *  $zeta \in [\log_2(minWavelength),\log_2(maxWavelength)]$ - zeta is log2 of
 * wavelength
 *
 *  The reason behind using zeta instead of wavenumber or wavelength is the
 * nature of wavelengths and it has better discretization properties. We want
 * to have a nice cascade of waves with exponentially increasing wavelengths.
 */

class WaveGrid {
public:
  using Vec2 = std::array<Real, 2>;
  using Vec3 = std::array<Real, 3>;
  using Vec4 = std::array<Real, 4>;

  using Idx = std::array<int, 4>;

  using Grid = Eigen::Tensor<Real, 4>;

public:
  struct Settings {
    // domain size
    Real size;
    Real max_wavelength;
    Real min_wavelength;

    // discretization nodes
    int n_x     = 100;
    int n_theta = 16;
    int n_zeta  = 1;

    // Additional settings
    Real initial_time = 0;

    // spectrum type
    enum SpectrumType {
      LinearBasis,
      PiersonMoskowitz
    } spectrumType = PiersonMoskowitz;
  };

public:
  WaveGrid(Settings s);

  void timeStep(const double dt);

  /** \brief Position and normal of water surface
   *
   * For initial position `pos` you get
   * \param pos Position where you
   */
  std::pair<Vec3, Vec3> waterSurface(Real x, Real y) const;

  Real cflTimeStep() const;

  Real amplitude(Real x, Real y, Real theta, Real zeta) const;

  std::vector<Vec4> trajectory(Vec4 pos4, Real length);

  void addPointDisturbance(const Vec2 pos, const double val);

public:
  // private:
  void advectionStep(const double dt);

  /*
   * \param gridPos position in 4D grid. In grid space coordinates!
   * \return position in 4D grid which corresponds to a reflected position
   * around boundary if necessary
   */
  Vec4 boundaryReflection(const Vec4 gridPos) const;

  auto interpolateAmplitude(Grid const &grid) const;

  bool inDomain(const Vec2 pos) const;
  Real levelset(const Vec2 pos) const;
  Vec2 levelsetGrad(const Vec2 pos) const;

  void diffusionStep(const double dt);

  void precompute1DBuffers(const double t);

public:
  // private:
  Real idxToPos(const int idx, const int dim) const;
  Vec4 idxToPos(const Idx idx) const;

  Real posToGrid(const Real pos, const int dim) const;
  Vec4 posToGrid(const Vec4 pos4) const;

  int posToIdx(const Real pos, const int dim) const;
  Idx posToIdx(const Vec4 pos4) const;

  Vec2 nodePosition(const int a_x, const int a_y) const;
  Real waveVectorAngle(const int angleIdx) const;
  Real waveNumber(const int waveNumberIdx) const;
  Vec2 waveDirection(const int angleIdx) const;
  Vec2 waveVector(const int a_x, const int a_y) const;

  Real dispersionRelation(const Real k) const;
  Real groupSpeed(const Real k) const;
  Real groupSpeed(const Vec4 pos4) const;
  Vec2 groupVelocity(const Real theta, const Real k) const;
  Vec2 groupVelocity(const Vec4 pos4) const;
  Real defaultAmplitude(const int b_theta, const int c_k) const;

  int  gridDim(const int dim) const;
  Real dx() const;

public:
  Grid m_amplitude, m_newAmplitude;

  std::vector<std::vector<Vec4>> m_profileBuffers;

  std::array<Real, 4> m_xmin;
  std::array<Real, 4> m_xmax;
  std::array<Real, 4> m_dx;
  std::array<Real, 4> m_idx;

  BoundaryLevelset m_boundary;

  // boundary info
  Vec2 islandCenter;
  Real islandRadius;
};
