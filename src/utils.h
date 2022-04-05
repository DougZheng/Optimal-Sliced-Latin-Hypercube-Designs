#pragma once

#include <cassert>

namespace SLHD {
namespace Utils {
static inline double QuickPow15(double x) {
  double y = x;
  y *= y, y *= y, y *= y, y *= y;
  return y / x;
}

static inline double QuickPow(double x, int p) {
  assert(p >= 0);
  if (p == 15) return QuickPow15(x);
  double y = 1;
  while (p > 0) {
    if (p & 1) y *= x;
    x *= x;
    p >>= 1;
  }
  return y;
}
} // namespace Utils
} // namespace SLHD