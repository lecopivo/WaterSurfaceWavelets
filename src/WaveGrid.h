#pragma once

#include "Enviroment.h"
#include "Global.h"
#include "Grid.h"
#include "ProfileBuffer.h"
#include "Spectrum.h"

namespace WaterWavelets {

/**
@ WaveGrid main class representing water surface.

@section Usage


@section Discretization

explain that \theta_0 = 0.5*dtheta

  The location is determined by four numbers: x,y,theta,zeta
  $x \in [-size,size]$ - the first spatial coordinate
  $y \in [-size,size]$ - the second spatial coordinate
  $theta \in [0,2 \pi)$  - direction of wavevector, theta==0 corresponds to
 wavevector in +x direction
  $zeta \in [\log_2(minWavelength),\log_2(maxWavelength)]$ - zeta is log2 of
 wavelength

  The reason behind using zeta instead of wavenumber or wavelength is the
 nature of wavelengths and it has better discretization properties. We want
 to have a nice cascade of waves with exponentially increasing wavelengths.
*/

class WaveGrid {
public:
  using Idx = std::array<int, 4>;

  enum Coord { X = 0, Y = 1, Theta = 2, Zeta = 3 };

public:
  /**
     @brief Settings to set up @ref WaveGrid

     The physical size of the resulting domain is:
     [-size,size]x[-size,size]x[0,2pi)x[min_zeta,max_zeta]

     The final grid resolution is: n_x*n_x*n_theta*n_zeta
   */
  struct Settings {
    /** Spatial size of the domain will be [-size,size]x[-size,size] */
    Real size = 50;
    /** Maximal zeta to simulate. */
    Real max_zeta = 0.01;
    /** Minimal zeta to simulate. */
    Real min_zeta = 10;

    /** Number of nodes per spatial dimension. */
    int n_x = 100;
    /** Number of discrete wave directions. */
    int n_theta = 16;
    /** Number of nodes in zeta. This determines resolution in wave number @see
     * @ref Discretization. */
    int n_zeta = 1;

    /** Set up initial time. Default value is 100 because at time zero you get
     * wierd coherency patterns */
    Real initial_time = 100;

    /** Select spectrum type. Currently only PiersonMoskowitz is supported. */
    enum SpectrumType {
      LinearBasis,
      PiersonMoskowitz
    } spectrumType = PiersonMoskowitz;
  };

public:
  /**
   * @brief Construct WaveGrid based on supplied @ref Settings
   * @param s settings to initialize WaveGrid
   */
  WaveGrid(Settings s);

  /** @brief Preform one time step.
   * @param dt time step
   * @param fullUpdate If true preform standard time step. If false, update only profile buffers.
   *
   * One time step consists of doing an advection step, diffusion step and
   * computation of profile buffers.
   *
   * To choose a reasonable dt we provide function @ref cflTimeStep()
   */
  void timeStep(const Real dt, bool fullUpdate =true);

  /**
   * @brief Position and normal of water surface
   * @param pos Position where you want to know position and normal
   * @return A pair: position, normal
   *
   * Returned position is not only a vertical displacement but because we use
   * Gerstner waves we get also a horizontal displacement.
   */
  std::pair<Vec3, Vec3> waterSurface(Vec2 pos) const;

  /**
   * @brief Time step based on CFL conditions
   *
   * Returns time in which the fastest waves move one grid cell across.
   * It is usefull for setting reasonable time step in @ref timeStep()
   */
  Real cflTimeStep() const;

  /**
   * @brief Amplitude at a point
   * @param pos4 Position in physical coordinates
   * @return Returns interpolated amplitude at a given point
   *
   * The point pos4 has to be physical coordinates: i.e. in box
   * [-size,size]x[-size,size]x[0,2pi)x[min_zeta,max_zeta] @ref Settings
   * Default value(@ref defaultAmplitude()) is returned is point is outside of
   * this box.
   */
  Real amplitude(Vec4 pos4) const;

  /**
   * @brief Amplitude value at a grid not
   * @param idx4 Node index to get amplitude at
   * @return Amplitude value at the give node.
   *
   * If idx4 is outside of discretized grid
   * {0,...,n_x-1}x{0,...,n_x-1}x{0,...,n_theta}x{0,...,n_zeta} then the default
   * value is returned,@ref defaultAmplitude().
   */
  Real gridValue(Idx idx4) const;

  /**
   * @brief Wave trajectory
   * @param pos4 Starting location
   * @param length Trajectory length
   * @return Trajectory of a wave starting at position pos4.
   *
   * This method was used for debugin purposes, mainly checking that boudary
   * reflection works correctly.
   */
  std::vector<Vec4> trajectory(Vec4 pos4, Real length) const;

  /**
   * @brief Adds point disturbance to a point
   * @param pos Position of disturbance, in physical coordinates.
   * @val Strength of disturbance
   *
   * This function basically increases amplitude for all directions at a certain
   * point.
   */
  void addPointDisturbance(Vec2 pos, Real val);

public:
  /**
   * @brief Preforms advection step
   * @param dt Time for one step.
   */
  void advectionStep(Real dt);

  /**
   * @brief Preforms diffusion step
   * @param dt Time for one step.
   */
  void diffusionStep(Real dt);

  /**
   * @brief Precomputes profile buffers
   *
   * The "parameter" to this function is the internal time(@ref m_time) at which
   * the profile buffers are precomputed.
   */
  void precomputeProfileBuffers();

  /**
   * @brief Precomputed group speed.
   *
   * This basically computes the "expected group speed" for the currently chosen
   * wave spectrum.
   */
  void precomputeGroupSpeeds();

  /**
   * @brief Boundary checking.
   * @param pos4 Point to be checked
   * @return If the intput point was inside of a boudary then this function
   * returns reflected point.
   *
   * If the point pos4 is not inside of the boundary then it is returned.
   *
   * If the point is inside of the boundary then it is reflected. This means
   * that the spatial position and also the wave direction is reflected w.t.r.
   * to the boundary.
   */
  Vec4 boundaryReflection(Vec4 pos4) const;

  /**
   * @brief Extends discrete grid with default values
   * @return Returns a function with signature(Int × Int × Int × Int -> Real).
   *
   * We store amplitudes on a 4D grid of the size (n_x*n_x*n_theta*n_zeta)(@ref
   * Settings). Sometimes it is usefull not to worry about grid bounds and get
   * default value for points outside of the grid or it wrapps arround for the
   * theta-coordinate.
   *
   * This function is doing exactly this, it returns a function which accepts
   * four integers and returns an amplitude. You do not have to worry about
   * the bounds.
   */
  auto extendedGrid() const;

  /**
   * @brief Preforms linear interpolation on the grid
   * @return Returns a function(Vec4 -> Real): accepting a point in physical
   * coordinates and returns interpolated amplitude.
   *
   * This function preforms a linear interpolation on the computational grid in
   * all four coordinates. It returns a function accepting @ref Vec4 and returns
   * interpolated amplitude. Any input point is valid because the interpolation
   * is done on @ref extendedGrid().
   */
  auto interpolatedAmplitude() const;

public:
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
  Real groupSpeed(int izeta) const;
  // Real groupSpeed(Vec4 pos4) const;

  // Vec2 groupVelocity(int izeta) const;
  Vec2 groupVelocity(Vec4 pos4) const;
  Real defaultAmplitude(int itheta, int izeta) const;

  int  gridDim(int dim) const;
  Real dx(int dim) const;

public:
  Grid     m_amplitude, m_newAmplitude;
  Spectrum m_spectrum;

  std::vector<ProfileBuffer> m_profileBuffers;

  std::array<Real, 4> m_xmin;
  std::array<Real, 4> m_xmax;
  std::array<Real, 4> m_dx;
  std::array<Real, 4> m_idx;

  std::vector<Real> m_groupSpeeds;

  Real m_time;

  Environment m_enviroment;
};

} // namespace WaterWavelets
