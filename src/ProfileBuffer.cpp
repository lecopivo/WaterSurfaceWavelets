#include "ProfileBuffer.h"

namespace WaterWavelets {

std::array<float, 4> ProfileBuffer::operator()(float p) const {

  const int N = m_data.size();

  // Guard from acessing outside of the buffer by wrapping
  auto extended_buffer = [=](int i) -> std::array<float, 4> {
    return m_data[pos_modulo(i, N)];
  };

  // Preform linear interpolation
  auto interpolated_buffer = LinearInterpolation(extended_buffer);

  // rescale `p` to interval [0,1)
  return interpolated_buffer(N * p / m_period);
}

float ProfileBuffer::dispersionRelation(float k) const {
  constexpr float g = 9.81;
  return sqrt(k * g);
};

std::array<float, 4> ProfileBuffer::gerstner_wave(float phase /*=knum*x*/,
                                                  float knum) const {
  float s = sin(phase);
  float c = cos(phase);
  return std::array<float, 4>{-s, c, -knum * c, -knum * s};
};

float ProfileBuffer::cubic_bump(float x) const {
  if (abs(x) >= 1)
    return 0.0f;
  else
    return x * x * (2 * abs(x) - 3) + 1;
};

}; // namespace WaterWavelets
