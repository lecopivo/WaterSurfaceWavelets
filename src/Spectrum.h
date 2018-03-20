#pragma once

#include <cmath>

class Spectrum {
public:
  Spectrum(double windSpeed) : m_windSpeed(windSpeed) {}

  /*
   * Maximal reasonable value of zeta to consider
   */
  double maxZeta() { return log2(10); }

  /*
   * Minamal resonable value of zeta to consider
   */
  double minZeta() { return log2(0.03); };

  /*
   * Returns density of wave for given zeta(=log2(wavelength))
   */
  double operator()(double zeta) {
    double A = pow(1.1, 1.5*zeta); // original pow(2, 1.5*zeta)
    double B = exp(-1.8038897788076411 * pow(4, zeta) / pow(m_windSpeed, 4));
    return 0.139098*sqrt(A*B);
  }

public:
  double m_windSpeed = 1;
};
