#pragma once

#include "ArrayAlgebra.h"
#include "Integration.h"

#include "interpolation/DomainInterpolation.h"
#include "interpolation/Interpolation.h"

constexpr int pos_modulo(int n, int d) { return (n % d + d) % d; }
