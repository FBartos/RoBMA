#include "StepRowNormalizer.h"

#include <JRmath.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace jags {
namespace RoBMA {
namespace {

constexpr double open_bound = 1e300;
constexpr double log_sqrt_two_pi = 0.91893853320467274178;
const double infinity = std::numeric_limits<double>::infinity();

// log(1 - exp(x)) for x <= 0 (Maechler 2012).
double log1mexp(double x)
{
  return x > -0.69314718055994530942 ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x));
}

// log P(a < Z < c) for a standard normal Z, on the log scale in either tail.
double log_interval(double a, double c)
{
  if (!(a < c)) return -infinity;
  if (a >= 0) {
    const double qa = pnorm(a, 0.0, 1.0, 0, 1);
    if (c == infinity) return qa;
    return qa + log1mexp(pnorm(c, 0.0, 1.0, 0, 1) - qa);
  }
  if (c <= 0) {
    const double pc = pnorm(c, 0.0, 1.0, 1, 1);
    if (a == -infinity) return pc;
    return pc + log1mexp(pnorm(a, 0.0, 1.0, 1, 1) - pc);
  }
  const double outside = (a == -infinity ? 0.0 : pnorm(a, 0.0, 1.0, 1, 0)) +
    (c == infinity ? 0.0 : pnorm(c, 0.0, 1.0, 0, 0));
  return std::log1p(-outside);
}

}

bool step_row_log_normalizer(StepRowLaw const &law, double mean, double *log_a, double *d1, double *d2)
{
  if (!(law.sd > 0) || !std::isfinite(law.sd) || !(law.se > 0) || !std::isfinite(law.se) ||
      !std::isfinite(mean) || law.bins < 1) {
    return false;
  }
  const double shift = law.effect_sign * mean;
  // the standardized bounds of each bin with weight
  double stack[3 * 16];
  std::vector<double> heap;
  double *score = stack;
  if (law.bins > 16) {
    heap.resize(3 * static_cast<std::size_t>(law.bins));
    score = heap.data();
  }
  double top = -infinity;
  for (int b = 0; b < law.bins; ++b) {
    double *lower = score + 3 * b, *upper = lower + 1, *log_term = lower + 2;
    *lower = law.lower[b] <= -open_bound ? -infinity : (law.lower[b] * law.se - shift) / law.sd;
    *upper = law.upper[b] >= open_bound ? infinity : (law.upper[b] * law.se - shift) / law.sd;
    *log_term = law.omega[b] > 0 ? std::log(law.omega[b]) + log_interval(*lower, *upper) : -infinity;
    top = std::max(top, *log_term);
  }
  if (!std::isfinite(top)) return false;
  long double total = 0;
  for (int b = 0; b < law.bins; ++b) total += std::exp(score[3 * b + 2] - top);
  *log_a = top + std::log(static_cast<double>(total));
  // ratios of the boundary densities to A: A'/A = -sign / sd * sum omega
  // (phi(upper) - phi(lower)) / A and A''/A = -sum omega (upper phi(upper) -
  // lower phi(lower)) / (A sd^2)
  long double s1 = 0, s2 = 0;
  for (int b = 0; b < law.bins; ++b) {
    if (!(law.omega[b] > 0)) continue;
    const double log_weight = std::log(law.omega[b]) - log_sqrt_two_pi - *log_a;
    const double lower = score[3 * b], upper = score[3 * b + 1];
    if (std::isfinite(upper)) {
      const double e = std::exp(log_weight - 0.5 * upper * upper);
      s1 += e;
      s2 += static_cast<long double>(upper) * e;
    }
    if (std::isfinite(lower)) {
      const double e = std::exp(log_weight - 0.5 * lower * lower);
      s1 -= e;
      s2 -= static_cast<long double>(lower) * e;
    }
  }
  const double ratio1 = -law.effect_sign * static_cast<double>(s1) / law.sd;
  const double ratio2 = -static_cast<double>(s2) / (law.sd * law.sd);
  *d1 = ratio1;
  *d2 = ratio2 - ratio1 * ratio1;
  return std::isfinite(*log_a) && std::isfinite(*d1) && std::isfinite(*d2);
}

}
}
