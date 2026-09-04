#include "DSELNORMFACTORSTEP.h"

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

DSELNORMFACTORSTEP::DSELNORMFACTORSTEP() :
  VectorDist("dselnorm_factor_step", 19) {}

bool DSELNORMFACTORSTEP::checkParameterLength(
    std::vector<unsigned int> const &len) const
{
  if (len[0] == 0 || len[1] != len[0] || len[2] % len[0] != 0 ||
      len[3] != len[0] || len[7] != len[0] || len[4] == 0 ||
      len[5] != len[4] || len[6] != len[4]) return false;
  const unsigned int rank = len[2] / len[0];
  if (rank < 2 || rank > 4) return false;
  for (unsigned int i = 8; i <= 10; ++i) {
    if (len[i] != 1) return false;
  }
  if (len[11] == 0 || len[12] != len[11] || len[13] < 3 || len[14] == 0) {
    return false;
  }
  for (unsigned int i = 15; i <= 18; ++i) {
    if (len[i] != 1) return false;
  }
  return true;
}

bool DSELNORMFACTORSTEP::checkParameterValue(
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  const int k = static_cast<int>(len[0]);
  const int rank = static_cast<int>(len[2] / len[0]);
  const int n_bins = static_cast<int>(len[4]);
  const int sign = static_cast<int>(*par[8]);
  const int telescope = static_cast<int>(*par[9]);
  const int kernel_mode = static_cast<int>(*par[10]);
  const int initial_points = static_cast<int>(*par[15]);
  const int max_points = static_cast<int>(*par[16]);
  const int scrambles = static_cast<int>(*par[17]);
  const double tolerance = *par[18];
  const SelNormJagsBounds z_lower(par[5], len[5]);
  const SelNormJagsBounds z_upper(par[6], len[6]);

  if (!std::isfinite(*par[8]) || *par[8] != static_cast<double>(sign) ||
      !std::isfinite(*par[9]) || *par[9] != static_cast<double>(telescope) ||
      !std::isfinite(*par[10]) ||
      *par[10] != static_cast<double>(kernel_mode) ||
      !std::isfinite(*par[15]) ||
      *par[15] != static_cast<double>(initial_points) ||
      !std::isfinite(*par[16]) ||
      *par[16] != static_cast<double>(max_points) ||
      !std::isfinite(*par[17]) ||
      *par[17] != static_cast<double>(scrambles) ||
      !(sign == 1 || sign == -1) || !(telescope == 0 || telescope == 1) ||
      !(kernel_mode == SELKERNEL_NORMAL || kernel_mode == SELKERNEL_STEP) ||
      initial_points < 4 || max_points < initial_points || scrambles < 2 ||
      !(tolerance > 0.0) ||
      !std::isfinite(tolerance) ||
      len[14] !=
        static_cast<unsigned int>(2 * rank * max_points * scrambles)) {
    return false;
  }
  for (int i = 0; i < k; ++i) {
    const int obs_bin = static_cast<int>(par[7][i]);
    if (!std::isfinite(par[0][i]) || !std::isfinite(par[1][i]) ||
        !(par[1][i] > 0.0) || !std::isfinite(par[3][i]) ||
        !(par[3][i] > 0.0) || !std::isfinite(par[7][i]) ||
        par[7][i] != static_cast<double>(obs_bin) ||
        obs_bin < 1 || obs_bin > n_bins) return false;
  }
  for (unsigned int i = 0; i < len[2]; ++i) {
    if (!std::isfinite(par[2][i])) return false;
  }
  for (int bin = 0; bin < n_bins; ++bin) {
    if (!std::isfinite(par[4][bin]) || par[4][bin] < 0.0 ||
        !(z_lower[bin] < z_upper[bin])) return false;
  }
  int quadrature_length = 0;
  int previous_order = 0;
  for (unsigned int rule = 0; rule < len[13]; ++rule) {
    const int order = static_cast<int>(par[13][rule]);
    if (!std::isfinite(par[13][rule]) ||
        par[13][rule] != static_cast<double>(order) ||
        order <= previous_order) return false;
    quadrature_length += order;
    previous_order = order;
  }
  if (quadrature_length != static_cast<int>(len[11])) return false;
  for (unsigned int i = 0; i < len[11]; ++i) {
    if (!std::isfinite(par[11][i]) || !std::isfinite(par[12][i])) {
      return false;
    }
  }
  // The QMC design is a large fixed JAGS data array. Scanning it here would
  // repeat invariant O(max_points * scrambles * rank) work for every density
  // update. The factor evaluator validates each coordinate lazily when a QMC
  // fallback actually consumes it; package-generated designs are validated at
  // construction.
  return selnorm_is_descending_step_partition(
    z_lower.data(), z_upper.data(), n_bins
  );
}

double DSELNORMFACTORSTEP::logDensity(
    double const *x, unsigned int length, PDFType type,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len,
    double const *lower, double const *upper) const
{
  const int k = static_cast<int>(length);
  const int rank = static_cast<int>(len[2] / len[0]);
  const int n_bins = static_cast<int>(len[4]);
  const int initial_points = static_cast<int>(*par[15]);
  const int max_points = static_cast<int>(*par[16]);
  const int scrambles = static_cast<int>(*par[17]);
  const SelNormJagsBounds z_lower(par[5], len[5]);
  const SelNormJagsBounds z_upper(par[6], len[6]);
  std::vector<int> obs_bin(static_cast<std::size_t>(k));
  for (int i = 0; i < k; ++i) {
    obs_bin[i] = static_cast<int>(par[7][i]);
  }

  double relative_mcse = 0.0;
  double relative_change = 0.0;
  const double log_density = cpp_selnorm_factor_step_lpdf(
    x, par[0], par[1], par[2], k, rank, par[3], par[4], n_bins,
    z_lower.data(), z_upper.data(), obs_bin.data(),
    static_cast<int>(*par[8]), static_cast<int>(*par[9]) == 1,
    static_cast<int>(*par[10]), par[11], par[12], par[13],
    static_cast<int>(len[13]), par[14], initial_points, max_points,
    scrambles, *par[18], &relative_mcse, &relative_change
  );
  const double diagnostic = std::max(relative_mcse, relative_change);
  if (!std::isfinite(diagnostic) || diagnostic > *par[18]) {
    throw std::runtime_error(
      "Exact selection factor normalizer was rejected by diagnostics: "
      "relative Monte Carlo standard error was " +
      std::to_string(relative_mcse) +
      " and nested-design relative change was " +
      std::to_string(relative_change) + ". Increase "
      "'max_points_per_scramble' or 'scrambles' in "
      "'selection_control'."
    );
  }
  if (!std::isfinite(log_density)) return JAGS_NEGINF;
  return log_density;
}

void DSELNORMFACTORSTEP::randomSample(
    double *x, unsigned int length,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len,
    double const *lower, double const *upper, RNG *rng) const
{
  for (unsigned int i = 0; i < length; ++i) {
    x[i] = std::numeric_limits<double>::quiet_NaN();
  }
}

void DSELNORMFACTORSTEP::support(
    double *lower, double *upper, unsigned int length,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  for (unsigned int i = 0; i < length; ++i) {
    lower[i] = JAGS_NEGINF;
    upper[i] = JAGS_POSINF;
  }
}

unsigned int DSELNORMFACTORSTEP::length(
    std::vector<unsigned int> const &len) const
{
  return len[0];
}

void DSELNORMFACTORSTEP::typicalValue(
    double *x, unsigned int length,
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len,
    double const *lower, double const *upper) const
{
  for (unsigned int i = 0; i < length; ++i) x[i] = par[0][i];
}

bool DSELNORMFACTORSTEP::isSupportFixed(
    std::vector<bool> const &fixmask) const
{
  return true;
}

}
}
