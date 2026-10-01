#include "DSELNORMMVSTEP.h"

#include <JRmath.h>
#include <rng/RNG.h>
#include <util/nainf.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "../selnorm/selnorm.h"
#include "../selnorm/selnorm-mv.h"
#include "selnorm-jags-bounds.h"

namespace jags {
namespace RoBMA {
DSELNORMMVSTEP::DSELNORMMVSTEP() : VectorDist("dselnorm_mnorm_step", 18) {}

bool DSELNORMMVSTEP::checkParameterLength(
    std::vector<unsigned int> const &len) const
{
  const unsigned int maximum = std::numeric_limits<int>::max();
  // Dense matrices, QMC offsets and quadrature offsets use signed int indices.
  if (len.size() != 18 || len[0] == 0 || len[0] > maximum / len[0] ||
      len[3] == 0 || len[3] > maximum || len[10] == 0 || len[10] > maximum ||
      len[15] == 0 || len[15] > maximum || len[17] < 3 || len[17] > maximum) {
    return false;
  }
  const std::uint64_t packed = static_cast<std::uint64_t>(len[0]) * (len[0] + 1) / 2;
  if (len[2] != len[0] || len[6] != len[0] || len[1] != packed ||
      len[4] != len[3] || len[5] != len[3]) {
    return false;
  }
  for (unsigned int i = 7; i < 10; ++i) {
    if (len[i] != 1) return false;
  }
  for (unsigned int i = 11; i < 15; ++i) {
    if (len[i] != 1) return false;
  }
  return len[15] == len[16];
}

bool DSELNORMMVSTEP::checkParameterValue(
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  return checkDynamicParameterValue(par, len) && checkControlParameterValue(par, len);
}

bool DSELNORMMVSTEP::checkDynamicParameterValue(
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  const int k = static_cast<int>(len[0]);
  const int n_bins = static_cast<int>(len[3]);
  if (!(*par[9] == SELKERNEL_NORMAL || *par[9] == SELKERNEL_STEP)) {
    return false;
  }
  for (int i = 0; i < k; ++i) {
    if (!std::isfinite(par[0][i])) return false;
  }
  for (unsigned int i = 0; i < len[1]; ++i) {
    if (!std::isfinite(par[1][i])) return false;
  }
  for (int bin = 0; bin < n_bins; ++bin) {
    if (!std::isfinite(par[3][bin]) || par[3][bin] < 0) return false;
  }
  return true;
}

bool DSELNORMMVSTEP::checkControlParameterValue(
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  const int k = static_cast<int>(len[0]);
  const int n_bins = static_cast<int>(len[3]);
  const double tolerance = *par[13];
  const int maximum = std::numeric_limits<int>::max();
  if (!(*par[7] == 1 || *par[7] == -1) ||
      !(*par[8] == 0 || *par[8] == 1) ||
      !selnorm_jags_integer_in_range(*par[11], 1, maximum) ||
      !selnorm_jags_integer_in_range(*par[12], 2, maximum) ||
      !selnorm_jags_integer_in_range(*par[14], SELVECTOR_PRODUCT, SELVECTOR_BEST_TWO_SIDED) ||
      !(tolerance > 0) || !std::isfinite(tolerance)) {
    return false;
  }
  // Bound every multiplication by the already validated QMC array length.
  // Malformed large controls must fail before narrowing or signed overflow.
  if (!selnorm_jags_qmc_length(len[10], 2U * len[0],
      static_cast<unsigned int>(*par[11]), static_cast<unsigned int>(*par[12]))) return false;
  const SelNormJagsBounds z_lower(par[4], len[4]);
  const SelNormJagsBounds z_upper(par[5], len[5]);
  for (int i = 0; i < k; ++i) {
    if (!std::isfinite(par[2][i]) || !(par[2][i] > 0) ||
        !selnorm_jags_integer_in_range(par[6][i], 1, n_bins)) {
      return false;
    }
  }
  for (int bin = 0; bin < n_bins; ++bin) {
    if (!(z_lower[bin] < z_upper[bin])) {
      return false;
    }
  }
  for (unsigned int i = 0; i < len[10]; ++i) {
    if (!std::isfinite(par[10][i]) || !(par[10][i] > 0) ||
        !(par[10][i] < 1)) {
      return false;
    }
  }
  int previous_order = 0;
  unsigned int total_order = 0;
  for (unsigned int i = 0; i < len[17]; ++i) {
    const double raw_order = par[17][i];
    if (!std::isfinite(raw_order) || raw_order != std::floor(raw_order) ||
        raw_order <= previous_order || raw_order > len[15] - total_order) return false;
    previous_order = static_cast<int>(raw_order);
    total_order += previous_order;
  }
  if (total_order != len[15]) return false;
  for (unsigned int i = 0; i < len[15]; ++i) {
    if (!std::isfinite(par[15][i]) || std::isnan(par[16][i]) ||
        par[16][i] == std::numeric_limits<double>::infinity()) return false;
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
  SelNormDenseIntegration integration;
  integration.nodes = par[15];
  integration.log_weights = par[16];
  integration.orders = par[17];
  integration.rule_count = len[17];
  integration.tolerance = *par[13];
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
    &relative_mcse, nullptr, nullptr, static_cast<int>(*par[14]), &integration
  );
  if (!std::isfinite(relative_mcse) || relative_mcse > *par[13]) {
    throw std::runtime_error(
      "Selection normalizer was rejected by diagnostics: relative "
      "Monte Carlo standard error was " + std::to_string(relative_mcse) +
      ". Increase 'points_per_scramble' or 'scrambles' in "
      "'selection_control'."
    );
  }
  if (!std::isfinite(log_density)) return JAGS_NEGINF;
  return log_density;
}

double DSELNORMMVSTEP::surrogateLogDensity(
    double const *x, unsigned int length, PDFType type,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len,
    double const *lower, double const *upper,
    SelNormCoarseSettings const &settings) const
{
  const int k = static_cast<int>(length);
  const int n_bins = static_cast<int>(len[3]);
  const SelNormJagsBounds z_lower(par[4], len[4]);
  const SelNormJagsBounds z_upper(par[5], len[5]);
  std::vector<int> obs_bin(k);
  for (int i = 0; i < k; ++i) obs_bin[i] = static_cast<int>(par[6][i]);
  SelNormDenseIntegration integration;
  integration.nodes = par[15];
  integration.log_weights = par[16];
  integration.orders = par[17];
  integration.rule_count = len[17];
  integration.tolerance = *par[13];
  return cpp_selnorm_mnorm_step_surrogate_lpdf(x, par[0], par[1], k, par[2], par[3],
    n_bins, z_lower.data(), z_upper.data(), obs_bin.data(), static_cast<int>(*par[7]),
    static_cast<int>(*par[8]) == 1, static_cast<int>(*par[9]), static_cast<int>(*par[14]),
    integration, settings);
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
