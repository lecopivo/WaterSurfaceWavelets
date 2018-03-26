#pragma once

#include <cmath>

class Spectrum {
public:
  Spectrum(double windSpeed);

  /*
   * Maximal reasonable value of zeta to consider
   */
  double maxZeta() const;

  /*
   * Minamal resonable value of zeta to consider
   */
  double minZeta() const;

  /*
   * Returns density of wave for given zeta(=log2(wavelength))
   */
  double operator()(double zeta) const;

public:
  double m_windSpeed = 1;
};
