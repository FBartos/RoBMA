#include "DSELNORMMVSTEP.h"

#include <JRmath.h>
#include <rng/RNG.h>
#include <util/nainf.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "../selnorm/selnorm.h"
#include "../selnorm/selnorm-mv.h"
#include "selnorm-jags-bounds.h"

namespace jags {
namespace RoBMA {

DSELNORMMVSTEP::DSELNORMMVSTEP() : VectorDist("dselnorm_mnorm_step", 14) {}

bool DSELNORMMVSTEP::checkParameterLength(
    std::vector<unsigned int> const &len) const
{
  if (len[0] == 0 || len[2] != len[0] || len[6] != len[0] ||
      len[1] != len[0] * (len[0] + 1) / 2 ||
      len[3] == 0 || len[4] != len[3] || len[5] != len[3]) {
    return false;
  }
  for (unsigned int i = 7; i < 10; ++i) {
    if (len[i] != 1) return false;
  }
  for (unsigned int i = 11; i < 14; ++i) {
    if (len[i] != 1) return false;
  }
  return len[10] > 0;
}

bool DSELNORMMVSTEP::checkParameterValue(
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  const int k = static_cast<int>(len[0]);
  const int n_bins = static_cast<int>(len[3]);
  const int sign = static_cast<int>(*par[7]);
  const int telescope = static_cast<int>(*par[8]);
  const int kernel_mode = static_cast<int>(*par[9]);
  const int points = static_cast<int>(*par[11]);
  const int scrambles = static_cast<int>(*par[12]);
  const double tolerance = *par[13];
  const SelNormJagsBounds z_lower(par[4], len[4]);
  const SelNormJagsBounds z_upper(par[5], len[5]);

  if (!(sign == 1 || sign == -1) ||
      !(telescope == 0 || telescope == 1) ||
      !(kernel_mode == SELKERNEL_NORMAL || kernel_mode == SELKERNEL_STEP) ||
      points < 1 || scrambles < 2 || !(tolerance > 0) ||
      !std::isfinite(tolerance) ||
      len[10] != static_cast<unsigned int>(2 * k * points * scrambles)) {
    return false;
  }
  for (int i = 0; i < k; ++i) {
    const int obs_bin = static_cast<int>(par[6][i]);
    if (!std::isfinite(par[0][i]) || !std::isfinite(par[2][i]) ||
        !(par[2][i] > 0) || obs_bin < 1 || obs_bin > n_bins) {
      return false;
    }
  }
  for (unsigned int i = 0; i < len[1]; ++i) {
    if (!std::isfinite(par[1][i])) return false;
  }
  for (int bin = 0; bin < n_bins; ++bin) {
    if (!std::isfinite(par[3][bin]) || par[3][bin] < 0 ||
        !(z_lower[bin] < z_upper[bin])) {
      return false;
    }
  }
  for (unsigned int i = 0; i < len[10]; ++i) {
    if (!std::isfinite(par[10][i]) || !(par[10][i] > 0) ||
        !(par[10][i] < 1)) {
      return false;
    }
  }
  return selnorm_is_descending_step_partition(
    z_lower.data(), z_upper.data(), n_bins
  );
}

double DSELNORMMVSTEP::logDensity(
    double const *x, unsigned int length, PDFType type,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len,
    double const *lower, double const *upper) const
{
  const int k = static_cast<int>(length);
  const int n_bins = static_cast<int>(len[3]);
  const int kernel_mode = static_cast<int>(*par[9]);
  const int points = static_cast<int>(*par[11]);
  const int scrambles = static_cast<int>(*par[12]);
  const SelNormJagsBounds z_lower(par[4], len[4]);
  const SelNormJagsBounds z_upper(par[5], len[5]);
  std::vector<int> obs_bin(static_cast<std::size_t>(k));
  for (int i = 0; i < k; ++i) {
    obs_bin[static_cast<std::size_t>(i)] = static_cast<int>(par[6][i]);
  }

  double relative_mcse = 0.0;
  const double log_density = cpp_selnorm_mnorm_step_lpdf(
    x,
    par[0],
    par[1],
    k,
    par[2],
    par[3],
    n_bins,
    z_lower.data(),
    z_upper.data(),
    obs_bin.data(),
    static_cast<int>(*par[7]),
    static_cast<int>(*par[8]) == 1,
    kernel_mode,
    par[10],
    points,
    scrambles,
    &relative_mcse
  );
  if (!std::isfinite(relative_mcse) || relative_mcse > *par[13]) {
    throw std::runtime_error(
      "Exact selection-normalizer integration failed its requested relative-error tolerance. Increase 'points_per_scramble' or 'scrambles' in 'selection_control'."
    );
  }
  if (!std::isfinite(log_density)) return JAGS_NEGINF;
  return log_density;
}

void DSELNORMMVSTEP::randomSample(
    double *x, unsigned int length,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len,
    double const *lower, double const *upper, RNG *rng) const
{
  for (unsigned int i = 0; i < length; ++i) {
    x[i] = std::numeric_limits<double>::quiet_NaN();
  }
}

void DSELNORMMVSTEP::support(
    double *lower, double *upper, unsigned int length,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  for (unsigned int i = 0; i < length; ++i) {
    lower[i] = JAGS_NEGINF;
    upper[i] = JAGS_POSINF;
  }
}

unsigned int DSELNORMMVSTEP::length(
    std::vector<unsigned int> const &len) const
{
  return len[0];
}

void DSELNORMMVSTEP::typicalValue(
    double *x, unsigned int length,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len,
    double const *lower, double const *upper) const
{
  for (unsigned int i = 0; i < length; ++i) x[i] = par[0][i];
}

bool DSELNORMMVSTEP::isSupportFixed(
    std::vector<bool> const &fixmask) const
{
  return true;
}

}
}
