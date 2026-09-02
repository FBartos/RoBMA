#include "DSELNORMCLUSTERSTEP.h"

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

DSELNORMCLUSTERSTEP::DSELNORMCLUSTERSTEP() :
  VectorDist("dselnorm_cluster_step", 15) {}

bool DSELNORMCLUSTERSTEP::checkParameterLength(
    std::vector<unsigned int> const &len) const
{
  if (len[0] == 0 || len[1] != len[0] || len[2] != len[0] ||
      len[3] != len[0] || len[7] != len[0] || len[4] == 0 ||
      len[5] != len[4] || len[6] != len[4] || len[11] == 0 ||
      len[12] != len[11] || len[13] < 2) return false;
  for (unsigned int i = 8; i <= 10; ++i) {
    if (len[i] != 1) return false;
  }
  return len[14] == 1;
}

bool DSELNORMCLUSTERSTEP::checkParameterValue(
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  const int k = static_cast<int>(len[0]);
  const int n_bins = static_cast<int>(len[4]);
  const int sign = static_cast<int>(*par[8]);
  const int telescope = static_cast<int>(*par[9]);
  const int kernel_mode = static_cast<int>(*par[10]);
  const double tolerance = *par[14];
  const SelNormJagsBounds z_lower(par[5], len[5]);
  const SelNormJagsBounds z_upper(par[6], len[6]);

  if (!(sign == 1 || sign == -1) || !(telescope == 0 || telescope == 1) ||
      !(kernel_mode == SELKERNEL_NORMAL || kernel_mode == SELKERNEL_STEP) ||
      !(tolerance > 0.0) || !std::isfinite(tolerance)) return false;
  for (int i = 0; i < k; ++i) {
    const int obs_bin = static_cast<int>(par[7][i]);
    if (!std::isfinite(par[0][i]) || !std::isfinite(par[1][i]) ||
        !(par[1][i] > 0.0) || !std::isfinite(par[2][i]) ||
        !std::isfinite(par[3][i]) || !(par[3][i] > 0.0) ||
        obs_bin < 1 || obs_bin > n_bins) return false;
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
    const int order = static_cast<int>(raw_order);
    if (!std::isfinite(raw_order) || raw_order != static_cast<double>(order) ||
        order <= previous_order) return false;
    quadrature_length += static_cast<unsigned int>(order);
    previous_order = order;
  }
  if (quadrature_length != len[11]) return false;
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

  double relative_change = 0.0;
  const double log_density = cpp_selnorm_cluster_step_lpdf(
    x, par[0], par[1], par[2], k, par[3], par[4], n_bins,
    z_lower.data(), z_upper.data(), obs_bin.data(),
    static_cast<int>(*par[8]), static_cast<int>(*par[9]) == 1,
    static_cast<int>(*par[10]), par[11], par[12], par[13],
    static_cast<int>(len[13]), *par[14],
    &relative_change
  );
  if (!std::isfinite(relative_change) || relative_change > *par[14]) {
    throw std::runtime_error(
      "Exact selection cluster normalizer was rejected by diagnostics: "
      "successive quadrature relative change was " +
      std::to_string(relative_change) + "."
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
