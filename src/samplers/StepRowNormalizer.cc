#include "StepRowNormalizer.h"
#include "../selnorm/selnorm.h"
#include "../selnorm/selnorm-tail.h"

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

namespace {

// Q(x) = P(Z > x) by the package's fast upper tail, infinite arguments included.
double upper(double x)
{
  if (x == infinity) return 0.0;
  if (x == -infinity) return 1.0;
  return selnorm_tail::scalar_kernel::upper_tail(x);
}

// A bin boundary on the standardized scale: its score, the tail beyond its
// magnitude Q(|score|) and the normal density there (zero at infinite scores).
struct Boundary {
  double score = 0, tail = 0, density = 0;
};

Boundary boundary(double score)
{
  constexpr double inverse_sqrt_two_pi = 0.39894228040143267794;
  Boundary out;
  out.score = score;
  if (std::isfinite(score)) {
    out.tail = upper(std::fabs(score));
    out.density = inverse_sqrt_two_pi * std::exp(-0.5 * score * score);
  }
  return out;
}

// The fast path: interval probabilities from the package's upper tail and the
// boundary densities in linear space, while A is far from underflow; false
// sends the row to the log-scale path. Adjacent bins of a step partition share
// their boundary, which is evaluated once.
bool linear_normalizer(StepRowLaw const &law, double shift, double *log_a, double *d1, double *d2)
{
  long double total = 0, s1 = 0, s2 = 0;
  Boundary last;
  bool have_last = false;
  auto at = [&](double score) {
    if (!have_last || last.score != score) {
      last = boundary(score);
      have_last = true;
    }
    return last;
  };
  for (int b = 0; b < law.bins; ++b) {
    if (!(law.omega[b] > 0)) continue;
    const Boundary up = at(law.upper[b] >= open_bound ? infinity : (law.upper[b] * law.se - shift) / law.sd);
    const Boundary low = at(law.lower[b] <= -open_bound ? -infinity : (law.lower[b] * law.se - shift) / law.sd);
    // P(low < Z < up) without cancellation: in the upper tail Q(low) - Q(up),
    // in the lower tail Phi(up) - Phi(low), across zero 1 - Q(up) - Phi(low)
    double probability;
    if (low.score >= 0) probability = low.tail - up.tail;
    else if (up.score <= 0) probability = up.tail - low.tail;
    else probability = 1.0 - up.tail - low.tail;
    const long double omega = law.omega[b];
    total += omega * probability;
    s1 += omega * (static_cast<long double>(up.density) - low.density);
    if (up.density > 0) s2 += omega * up.score * up.density;
    if (low.density > 0) s2 -= omega * low.score * low.density;
  }
  if (!(total > 1e-280L) || !std::isfinite(static_cast<double>(total))) return false;
  *log_a = std::log(static_cast<double>(total));
  const double ratio1 = -law.effect_sign * static_cast<double>(s1 / total) / law.sd;
  const double ratio2 = -static_cast<double>(s2 / total) / (law.sd * law.sd);
  *d1 = ratio1;
  *d2 = ratio2 - ratio1 * ratio1;
  return std::isfinite(*d1) && std::isfinite(*d2);
}

}

bool step_row_log_normalizer(StepRowLaw const &law, double mean, double *log_a, double *d1, double *d2)
{
  if (!(law.sd > 0) || !std::isfinite(law.sd) || !(law.se > 0) || !std::isfinite(law.se) ||
      !std::isfinite(mean) || law.bins < 1) {
    return false;
  }
  const double shift = law.effect_sign * mean;
  if (linear_normalizer(law, shift, log_a, d1, d2)) return true;
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

double step_row_log_density(StepRowLaw const &law, double power, int bin, int mode, bool telescope,
                            double x, double mean)
{
  if (mode == SELKERNEL_NORMAL || bin == 0) {
    const double log_lik = dnorm(x, mean, law.sd, true);
    return power == 1 ? log_lik : power * log_lik;
  }
  double phack_z_zero[2] = {0, 0};
  double segment_bounds_zero[1] = {0};
  int segment_zero[1] = {0};
  SelNormKernelData data;
  data.n_bins = law.bins;
  data.n_segments = 0;
  data.effect_sign = law.effect_sign;
  data.q = 0;
  data.z_lower = law.lower;
  data.z_upper = law.upper;
  data.phack_z_source = phack_z_zero;
  data.phack_z_dest = phack_z_zero;
  data.segment_bounds = segment_bounds_zero;
  data.segment_step_bin = segment_zero;
  data.segment_phack_region = segment_zero;
  data.segment_step_bin_real = 0;
  data.segment_phack_region_real = 0;
  data.trusted_step_partition = true;
  data.telescope_probabilities = telescope;
  return cpp_selnorm_kernel_lpdf(x, mean, law.sd, mean, law.sd, law.se, power, law.omega, bin, 0, 0,
                                 SELKERNEL_STEP, data, 1, false);
}

}
}
