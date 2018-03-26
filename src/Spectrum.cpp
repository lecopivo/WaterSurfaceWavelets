#include "Spectrum.h"

Spectrum::Spectrum(double windSpeed) : m_windSpeed(windSpeed) {}

double Spectrum::maxZeta() const { return log2(10); }

double Spectrum::minZeta() const { return log2(0.03); };

double Spectrum::operator()(double zeta) const {
  double A = pow(1.1, 1.5 * zeta); // original pow(2, 1.5*zeta)
  double B = exp(-1.8038897788076411 * pow(4, zeta) / pow(m_windSpeed, 4));
  return 0.139098 * sqrt(A * B);
}
