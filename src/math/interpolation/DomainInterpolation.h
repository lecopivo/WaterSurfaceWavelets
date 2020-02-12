#pragma once

#include <iostream>

#include <Eigen/Dense>

/*
 * Domain inteprolation - Interpolates function values only for pointS in the
 * domain
 *
 * \tparam Interpolation This is any function satisfying concept Interpolatiobn
 * \tparam Domain This is bool valued function on integers returning true for
 points inside of the domain and false otherwise
 * \param interpolation Interpolation to use.
 * \param domain Function indicating domain.
 */
template <class Interpolation, class Domain>
auto DomainInterpolation(Interpolation interpolation, Domain domain) {
  return [=](auto fun) mutable {

    // The function `dom_fun` collects function values and weights inside of
    // the domain
    auto dom_fun = [=](auto... x) mutable -> Eigen::Vector2d {
      Eigen::Vector2d val;
      if (domain(x...) == true) {
        val << fun(x...), 1.0;
      } else {
        val << 0.0, 0.0;
      }
      return val;
    };

    // Interpolates `dom_fun`
    auto int_fun = interpolation(dom_fun);

    return [=](auto... x) mutable {
      auto val = int_fun(x...);
      double f = val[0];
      double w = val[1];
      return w != 0.0 ? f / w : 0.0;
    };
  };
}
