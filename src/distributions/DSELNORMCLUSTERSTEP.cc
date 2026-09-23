#include "DSELNORMCLUSTERSTEP.h"

#include <rng/RNG.h>
#include <util/nainf.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "../selnorm/selnorm.h"
#include "../selnorm/selnorm-mv.h"
#include "selnorm-jags-bounds.h"

namespace jags {
namespace RoBMA {

DSELNORMCLUSTERSTEP::DSELNORMCLUSTERSTEP() :
  VectorDist("dselnorm_cluster_step", 20) {}

bool DSELNORMCLUSTERSTEP::checkParameterLength(
    std::vector<unsigned int> const &len) const
{
  if (len.size() != 20) return false;
  const unsigned int maximum = std::numeric_limits<int>::max();
  for (unsigned int length : len) if (length > maximum) return false;
  if (len[0] == 0 || len[0] > maximum / len[0]) return false;
  if (len[0] == 0 || len[1] != len[0] || len[2] != len[0] ||
      len[3] != len[0] || len[7] != len[0] || len[4] == 0 ||
      len[5] != len[4] || len[6] != len[4] || len[11] == 0 ||
      len[12] != len[11] || len[13] < 2 || len[14] == 0) return false;
  for (unsigned int i = 8; i <= 10; ++i) {
    if (len[i] != 1) return false;
  }
  for (unsigned int i = 15; i <= 19; ++i) {
    if (len[i] != 1) return false;
  }
  return true;
}

bool DSELNORMCLUSTERSTEP::checkParameterValue(
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  const int k = static_cast<int>(len[0]);
  const int n_bins = static_cast<int>(len[4]);
  const int limit = std::numeric_limits<int>::max();
  const double tolerance = *par[18];
  if (!(*par[8] == 1 || *par[8] == -1) || !(*par[9] == 0 || *par[9] == 1) ||
      !(*par[10] == SELKERNEL_NORMAL || *par[10] == SELKERNEL_STEP) ||
      !selnorm_jags_integer_in_range(*par[15], 4, limit) ||
      !selnorm_jags_integer_in_range(*par[16], 4, limit) || *par[16] < *par[15] ||
      !selnorm_jags_integer_in_range(*par[17], 2, limit) ||
      !selnorm_jags_integer_in_range(*par[19], SELVECTOR_PRODUCT, SELVECTOR_BEST_TWO_SIDED) ||
      !(tolerance > 0.0) || !std::isfinite(tolerance)) return false;
  if (!selnorm_jags_qmc_length(len[14], 2,
      static_cast<unsigned int>(*par[16]), static_cast<unsigned int>(*par[17]))) return false;
  const SelNormJagsBounds z_lower(par[5], len[5]);
  const SelNormJagsBounds z_upper(par[6], len[6]);
  for (int i = 0; i < k; ++i) {
    if (!std::isfinite(par[0][i]) || !std::isfinite(par[1][i]) ||
        !(par[1][i] > 0.0) || !std::isfinite(par[2][i]) ||
        !std::isfinite(par[3][i]) || !(par[3][i] > 0.0) ||
        !selnorm_jags_integer_in_range(par[7][i], 1, n_bins)) return false;
  }
  for (int bin = 0; bin < n_bins; ++bin) {
    if (!std::isfinite(par[4][bin]) || par[4][bin] < 0.0 ||
        !(z_lower[bin] < z_upper[bin])) return false;
  }
  for (unsigned int i = 0; i < len[11]; ++i) {
    if (!std::isfinite(par[11][i]) || !std::isfinite(par[12][i])) return false;
  }
  unsigned int quadrature_length = 0;
  int previous_order = 0;
  for (unsigned int i = 0; i < len[13]; ++i) {
    const double raw_order = par[13][i];
    if (!selnorm_jags_integer_in_range(raw_order, 1, limit) ||
        raw_order <= previous_order || raw_order > len[11] - quadrature_length) return false;
    const int order = static_cast<int>(raw_order);
    quadrature_length += static_cast<unsigned int>(order);
    previous_order = order;
  }
  if (quadrature_length != len[11]) return false;
  // As in the higher-rank factor route, the fixed QMC array is validated at
  // construction and each coordinate is checked only if fallback consumes it.
  return selnorm_is_descending_step_partition(
    z_lower.data(), z_upper.data(), n_bins
  );
}

double DSELNORMCLUSTERSTEP::logDensity(
    double const *x, unsigned int length, PDFType type,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len,
    double const *lower, double const *upper) const
{
  const int k = static_cast<int>(length);
  const int n_bins = static_cast<int>(len[4]);
  const SelNormJagsBounds z_lower(par[5], len[5]);
  const SelNormJagsBounds z_upper(par[6], len[6]);
  std::vector<int> obs_bin(static_cast<std::size_t>(k));
  for (int i = 0; i < k; ++i) {
    obs_bin[static_cast<std::size_t>(i)] = static_cast<int>(par[7][i]);
  }

  double relative_mcse = 0.0;
  double relative_change = 0.0;
  const double log_density = cpp_selnorm_cluster_step_lpdf(
    x, par[0], par[1], par[2], k, par[3], par[4], n_bins,
    z_lower.data(), z_upper.data(), obs_bin.data(),
    static_cast<int>(*par[8]), static_cast<int>(*par[9]) == 1,
    static_cast<int>(*par[10]), par[11], par[12], par[13],
    static_cast<int>(len[13]), par[14], static_cast<int>(*par[15]),
    static_cast<int>(*par[16]), static_cast<int>(*par[17]), *par[18],
    &relative_mcse, &relative_change, nullptr, nullptr, static_cast<int>(*par[19])
  );
  if (!std::isfinite(relative_mcse) || !std::isfinite(relative_change) ||
      std::max(relative_mcse, relative_change) > *par[18]) {
    throw std::runtime_error(
      "Selection cluster normalizer was rejected by diagnostics: "
      "relative Monte Carlo standard error was " +
      std::to_string(relative_mcse) +
      " and nested-design relative change was " +
      std::to_string(relative_change) + ". Increase "
      "'max_points_per_scramble' or 'scrambles' in 'selection_control'."
    );
  }
  if (!std::isfinite(log_density)) return JAGS_NEGINF;
  return log_density;
}

void DSELNORMCLUSTERSTEP::randomSample(
    double *x, unsigned int length,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len,
    double const *lower, double const *upper, RNG *rng) const
{
  for (unsigned int i = 0; i < length; ++i) {
    x[i] = std::numeric_limits<double>::quiet_NaN();
  }
}

void DSELNORMCLUSTERSTEP::support(
    double *lower, double *upper, unsigned int length,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  for (unsigned int i = 0; i < length; ++i) {
    lower[i] = JAGS_NEGINF;
    upper[i] = JAGS_POSINF;
  }
}

unsigned int DSELNORMCLUSTERSTEP::length(
    std::vector<unsigned int> const &len) const
{
  return len[0];
}

void DSELNORMCLUSTERSTEP::typicalValue(
    double *x, unsigned int length,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len,
    double const *lower, double const *upper) const
{
  for (unsigned int i = 0; i < length; ++i) x[i] = par[0][i];
}

bool DSELNORMCLUSTERSTEP::isSupportFixed(
    std::vector<bool> const &fixmask) const
{
  return true;
}

}
}
