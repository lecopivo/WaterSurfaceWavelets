#pragma once

#include <cmath>

double sine_wave(double a, double k, double kx, double t) {

  constexpr double g = 9.81;

  double omega = sqrt(g / k);
  double theta = kx - omega * t;
  double nu    = a * cos(theta);
  
  return nu;
}

double stokes_wave(double a, double k, double kx, double t) {

  constexpr double g = 9.81;

  double omega = k * (1 + 0.5 * pow(k * a, 2)) * sqrt(g / k);
  double theta = kx - omega * t;
  double nu    = a * (cos(theta) + 0.5 * (k * a) * cos(2 * theta) +
                   3.0 / 8.0 * pow(k * a, 2) * cos(3 * theta));

  return nu;
}

std::pair<double, double> gerstner_wave(double a, double k, double kx,
                                        double t) {

  constexpr double g = 9.81;

  double omega = sqrt(g * k);
  double theta = kx - omega * t;
  double X     = a / k * sin(theta);
  double Y     = -a / k * cos(theta);
  
  return {X, Y};
}
