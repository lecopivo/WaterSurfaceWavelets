#include <array>
#include <iostream>
#include <vector>

template <int DIM, typename T>
class Grid {
public:
  static constexpr int Dim = DIM;
  using Idx                = std::array<int, Dim>;

public:
  Grid(const Idx sizes)
      : m_sizes(sizes)
      , m_offsets([sizes]() {
        Idx offsets;
        offsets[Dim - 1] = 1;
        for (int i = Dim - 2; i >= 0; i--) {
          offsets[i] = offsets[i + 1] * sizes[i + 1];
        }
        return offsets;
      }())
      , m_data(m_offsets[Dim - 1]) {}

  T const &operator()(const Idx idx) const { return m_data[linearIndex(idx)]; }

  T &operator()(const Idx idx) { return m_data[linearIndex(idx)]; }

  int gridSize(const int dim) const { return m_sizes[dim]; }

private:
  int linearIndex(const Idx idx) const {
    int index = 0;
    for (int i = 0; i < Dim; i++) {
      index += idx[i] * m_offsets[i];
    }
    return index;
  }

private:
  const Idx      m_sizes;
  const Idx      m_offsets;
  std::vector<T> m_data;
};

class WaveGrid {
public:
  using Pos  = std::array<double, 2>;
  using Dir  = std::array<double, 2>;
  using Vec4 = std::array<double, 4>;

  using PosIdx = std::array<int, 2>;
  using DirIdx = std::array<int, 2>;
  using Idx    = std::array<int, 4>;

public:
  void TimeStep(const double dt, const double t) {
    AdvectionStep(dt);
    DiffusionStep(dt);
    Precompute1DBuffers(t);
  }

  void AdvectionStep(const double dt) {
    for (int i = 0; i < gridSize(0); i++) {
      for (int j = 0; i < gridSize(1); j++) {
        for (int k = 0; i < gridSize(2); k++) {
          for (int l = 0; i < gridSize(3); l++) {

            Pos x     = gridPosition({i, j});
            Dir v     = groupVelocity({k, l});
            Pos new_x = Pos{x[0] - dt * v[0], x[1] - dt * v[1]};

            HandleBoundaryConditions(x, new_x, v);

            m_newAmplitude({i, j, k, l}) =
                interpolate(m_amplitude, {new_x[0], new_x[1], v[0], v[1]});
          }
        }
      }
    }

    std::swap(m_newAmplitude, m_amplitude);
  }

  void DiffusionStep(const double dt) {
    for (int i = 0; i < gridSize(0); i++) {
      for (int j = 0; i < gridSize(1); j++) {
        for (int k = 0; i < gridSize(2); k++) {
          for (int l = 0; i < gridSize(3); l++) {

            double gamma_k     = 0.1;
            double gamma_alpha = 0.1;

            m_newAmplitude({i, j, k, l}) =
                (1 - gamma_k - gamma_alpha) * m_amplitude({i, j, k, l}) +
                0.5 * gamma_k *
                    (m_amplitude({i, j, k + 1, l}) +
                     m_amplitude({i, j, k - 1, l})) +
                0.5 * gamma_alpha *
                    (m_amplitude({i, j, k, l + 1}) +
                     m_amplitude({i, j, k, l + 1}));
          }
        }
      }
    }
  }

  void Precompute1DBuffers(const double dt) {
    for (int k = 0; k < gridSize(2); k++) {
    }
  }

  void HandleBoundaryConditions(const Pos x, Pos &new_x, Dir &v) {

    if (CrossedBoundary(x, new_x)) {

      Pos intersection;
      Dir normal;

      BoundaryIntersection(x, new_x, intersection, normal);

      new_x += -2.0 * normal * normal.dot(new_x - intersection);
      v += -2.0 * normal * normal.dot(v);
    }
  }

  double WaterHeight(const Pos x) { return 0.0; }

public:
  Pos gridPosition(const PosIdx idx) {}

  Dir groupVelocity(const DirIdx idx) {}

  int gridSize(const int dim) const { return m_amplitude.gridSize(dim); }

private:
  Grid<4, double> m_amplitude, m_newAmplitude;
};

int main() { return 0; }
