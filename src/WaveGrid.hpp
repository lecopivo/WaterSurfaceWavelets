#pragma once

#include <array>
#include <iostream>
#include <vector>

#include <Eigen/../unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

#include <Interpolation/DomainInterpolation.h>
#include <Interpolation/Interpolation.h>

using Real = double;

using Vec2 = Eigen::Matrix<Real, 2, 1>;
using Vec4 = Eigen::Matrix<Real, 4, 1>;

using PositionIdx = std::array<int, 2>;
using WaveVectorIdx = std::array<int, 2>;
using Idx = std::array<int, 4>;

using Grid = Eigen::Tensor<Real, 4>;

constexpr int pos_modulo(int n, int d);

class WaveGrid {
public:
  struct Settings {
    // domain size
    Real xmin;
    Real xmax;
    Real ymin;
    Real ymax;

    // discretization nodes
    int n_x;
    int n_y;
    int n_theta;
    int n_k;

    // spectrum type
    enum { HAHA_SPECTRUM } spectrumType;
  };

public:
  WaveGrid(Settings s);

  void timeStep(const double dt, const double t);

  Real waterHeight(const Vec2 pos);

  Real cflTimeStep() const;

  Real amplitude(const Vec4 pos) const;

  std::vector<Vec4> trajectory(Vec4 pos4, Real length);

  void addPointDisturbance(const Vec2 pos,const double val);
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

  int gridDim(const int dim) const;
  Real dx() const;

public:
  Grid m_amplitude, m_newAmplitude;

  std::array<Real, 4> m_xmin;
  std::array<Real, 4> m_xmax;
  std::array<Real, 4> m_dx;
  std::array<Real, 4> m_idx;

  // boundary info
  Vec2 islandCenter;
  Real islandRadius;
};
