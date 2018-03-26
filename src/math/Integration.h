#pragma once

#include "ArrayAlgebra.h"

template <typename Fun>
auto integrate(int integration_nodes, double x_min, double x_max, Fun const &fun) {

  double dx = (x_max - x_min) / integration_nodes;
  double x  = x_min + 0.5 * dx;

  auto result = dx * fun(x); // the first integration node
  for (int i = 1; i < integration_nodes; i++) { // proceed with other nodes, notice `int i= 1`
    x += dx;
    result += dx * fun(x);
  }

  return result;
}
