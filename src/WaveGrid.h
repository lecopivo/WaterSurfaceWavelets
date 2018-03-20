#pragma once

// #include <Eigen/../unsupported/Eigen/CXX11/Tensor>
// #include <Eigen/Dense>

#include "Enviroment.h"
#include "global.h"

namespace WaterWavelets {


  // using Grid = Eigen::Tensor<Real, 4>;

  class Grid{
  public:

    Grid();

    void resize(int n0, int n1, int n2, int n3);

    Real& operator()(int i0, int i1, int i2, int i3);

    Real const& operator()(int i0, int i1, int i2, int i3) const;

    int dimension(int dim) const;

  private:
    std::vector<Real> m_data;
    std::array<int,4> m_dimensions;
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

class WaveGrid{
public:
  using Idx = std::array<int, 4>;

  enum Coord { X = 0, Y = 1, Theta = 2, Zeta = 3 };

public:
  struct Settings {
    // domain size
    Real size = 50;
    Real max_wavelength = 0.01;
    Real min_wavelength = 10;

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

  void timeStep(const Real dt);

  /** \brief Position and normal of water surface
   *
   * For initial position `pos` you get
   * \param pos Position where you
   */
  std::pair<Vec3, Vec3> waterSurface(Vec2 pos) const;

  Real cflTimeStep() const;

  // Real amplitude(Real x, Real y, Real theta, Real zeta) const;
  Real amplitude(Vec4 pos4) const;

  // Real gridValue(int ix, int iy, int itheta, int izeta) const;
  Real gridValue(Idx idx4) const;

  std::vector<Vec4> trajectory(Vec4 pos4, Real length) const;

  void addPointDisturbance(Vec2 pos, Real val);

public:
  // private:
  void advectionStep(Real dt);

  /*
   * \param gridPos position in 4D grid. In grid space coordinates!
   * \return position in 4D grid which corresponds to a reflected position
   * around boundary if necessary
   */
  Vec4 boundaryReflection(Vec4 gridPos) const;

  /*
   * Returns: Int × Int × Int × Int -> Real
   */
  auto extendedGrid() const;

  /*
   * Returns interpolated amplitude
   * The return value is a function, `Vec4 -> Real`, taking a position in the
   * physical space and returning amplitude
   *
   */
  auto interpolatedAmplitude() const;

  void diffusionStep(Real dt);

  Real bufferPeriod(int izeta) const;
  void precomputeProfileBuffers(Real time);

  /*
   * Returns $\Psi_i$, Defined by equation (21).
   * The return value is a lambda representing function $\Psi_i$.
   *
   * The lambda accepts a single number, a dot product between position and wave
   direction, and returns four numbers:
   * 1. horizontal offset
   * 2. vertical offset
   * 3. derivative of horizontal offset
   * 4. derivative of vertical offset
   *
   * Assumptions: The function `precomputeProfileBuffers` has been called. That
   function actually preforms the numerical integration and stores the result
   into internal buffers. If you step in time with `timeStep` then
   `precomputeProfileBuffers` gets called automatically and you don't have to
   worry about it.
   */
  auto profileBuffer(int i) const;

public:
  // private:
  Real idxToPos(int idx, int dim) const;
  Vec4 idxToPos(Idx idx) const;

  Real posToGrid(Real pos, int dim) const;
  Vec4 posToGrid(Vec4 pos4) const;

  int posToIdx(Real pos, int dim) const;
  Idx posToIdx(Vec4 pos4) const;

  Vec2 nodePosition(int ix, int iy) const;

  Real waveLength(int izeta) const;
  Real waveNumber(int izeta) const;

  Real dispersionRelation(Real k) const;
  Real dispersionRelation(Vec4 pos4) const;
  Real groupSpeed(Vec4 pos4) const;

  Vec2 groupVelocity(Vec4 pos4) const;
  Real defaultAmplitude(int itheta, int izeta) const;

  int  gridDim(int dim) const;
  Real dx(int dim) const;

public:
  Grid m_amplitude, m_newAmplitude;

  std::vector<std::vector<Vec4>> m_profileBuffers;

  std::array<Real, 4> m_xmin;
  std::array<Real, 4> m_xmax;
  std::array<Real, 4> m_dx;
  std::array<Real, 4> m_idx;

  Real m_time;

  Enviroment m_enviroment;
};

} // namespace WaterWavelets
