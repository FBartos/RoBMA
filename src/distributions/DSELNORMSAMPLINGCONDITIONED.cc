#include "DSELNORMSAMPLINGCONDITIONED.h"
#include "../selnorm/selnorm.h"
#include "../selnorm/selnorm-mv.h"
#include "selnorm-jags-bounds.h"

#include <util/nainf.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace jags {
namespace RoBMA {

DSELNORMSAMPLINGCONDITIONED::DSELNORMSAMPLINGCONDITIONED() :
  VectorDist("dselnorm_sampling_conditioned", 29) {}

bool DSELNORMSAMPLINGCONDITIONED::checkParameterLength(
    std::vector<unsigned int> const &len) const
{
  if (len.size() != 29) return false;
  const unsigned int maximum = std::numeric_limits<int>::max();
  for (unsigned int length : len) if (length > maximum) return false;
  if (len[0] == 0 || len[0] > maximum / len[0]) return false;
  if (len[0] == 0 || len[1] != static_cast<std::uint64_t>(len[0]) * (len[0] + 1) / 2 ||
      len[2] != len[0] || len[3] == 0 || len[4] != 1 ||
      len[5] != len[0] || len[6] != len[0] || len[7] != len[0] ||
      len[8] == 0 || len[9] != len[8] || len[10] != len[8] ||
      len[11] != len[0] || len[16] != len[0] || len[17] == 0 ||
      len[22] == 0 || len[23] != len[22] || len[24] < 3 ||
      len[25] == 0 || len[26] != len[25] || len[27] < 3 || len[28] != 3) return false;
  for (int i = 12; i <= 15; ++i) if (len[i] != 1) return false;
  for (int i = 18; i <= 21; ++i) if (len[i] != 1) return false;
  return true;
}

bool DSELNORMSAMPLINGCONDITIONED::checkParameterValue(
    std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  const int k = static_cast<int>(len[0]);
  const int limit = std::numeric_limits<int>::max();
  if (!selnorm_jags_integer_in_range(*par[4], 0, limit / k) ||
      !(*par[12] == 1 || *par[12] == -1) ||
      !(*par[13] == SELKERNEL_NORMAL || *par[13] == SELKERNEL_STEP) ||
      !(*par[14] == 0 || *par[14] == 1) ||
      !selnorm_jags_integer_in_range(*par[15], SELVECTOR_PRODUCT, SELVECTOR_BEST_TWO_SIDED) ||
      !std::isfinite(*par[21]) || !(*par[21] > 0.0)) return false;
  const int rank = static_cast<int>(*par[4]);
  if (rank > 0 && rank > limit / rank) return false;
  if (len[3] != static_cast<unsigned int>(std::max(1, k * rank))) return false;
  if (!selnorm_jags_integer_in_range(*par[18], 1, limit) ||
      !selnorm_jags_integer_in_range(*par[19], 1, limit) || *par[19] < *par[18] ||
      !selnorm_jags_integer_in_range(*par[20], 2, limit)) return false;
  if (!(rank == 0 && len[17] == 1) &&
      !selnorm_jags_qmc_length(len[17], std::max(2 * k, rank),
        static_cast<unsigned int>(*par[19]), static_cast<unsigned int>(*par[20]))) return false;
  for (int row = 0; row < k; ++row) {
    if (!std::isfinite(par[0][row]) || !std::isfinite(par[2][row]) ||
        par[2][row] < 0.0 || !std::isfinite(par[5][row]) ||
        !std::isfinite(par[6][row]) || !std::isfinite(par[7][row]) ||
        !(par[7][row] > 0.0) ||
        !selnorm_jags_integer_in_range(par[11][row], 1, len[8]) ||
        !selnorm_jags_integer_in_range(par[16][row], 1, limit)) return false;
  }
  for (int argument : {1, 3, 8, 17}) {
    for (unsigned int i = 0; i < len[argument]; ++i) {
      if (!std::isfinite(par[argument][i])) return false;
      if (argument == 8 && par[argument][i] < 0.0) return false;
      if (argument == 17 && !(par[argument][i] > 0.0 && par[argument][i] < 1.0)) return false;
    }
  }
  for (int argument = 22; argument < 29; ++argument) {
    for (unsigned int i = 0; i < len[argument]; ++i) {
      if (!std::isfinite(par[argument][i])) return false;
    }
  }
  for (int offset : {22, 25}) {
    int total = 0, previous = 0;
    for (unsigned int i = 0; i < len[offset + 2]; ++i) {
      const double order = par[offset + 2][i];
      if (!selnorm_jags_integer_in_range(order, 1, limit) ||
          order <= previous || order > len[offset] - total) return false;
      previous = static_cast<int>(order);
      total += previous;
    }
    if (len[offset] != static_cast<unsigned int>(total)) return false;
  }
  for (int rank_index = 0; rank_index < 3; ++rank_index) {
    const double count = par[28][rank_index];
    if (count != std::floor(count) || count < 3 || count > len[27]) return false;
  }
  const SelNormJagsBounds lower(par[9], len[9]), upper(par[10], len[10]);
  return selnorm_is_descending_step_partition(lower.data(), upper.data(), len[8]);
}

double DSELNORMSAMPLINGCONDITIONED::logDensity(double const *x,
    unsigned int length, PDFType type, std::vector<double const *> const &par,
    std::vector<unsigned int> const &len, double const *lower,
    double const *upper) const
{
  const int k = static_cast<int>(length);
  const SelNormJagsBounds z_lower(par[9], len[9]), z_upper(par[10], len[10]);
  std::vector<int> bins(k), groups(k);
  std::vector<double> e(k), delta(k);
  for (int row = 0; row < k; ++row) {
    bins[row] = static_cast<int>(par[11][row]);
    groups[row] = static_cast<int>(par[16][row]);
  }
  double mcse = 0.0, change = 0.0, normalizer = 0.0;
  const SelNormConditionedQuadrature quadrature = {par[22], par[23], par[24],
    static_cast<int>(len[24]), par[25], par[26], par[27],
    static_cast<int>(len[27]), par[28]};
  const double value = cpp_selnorm_sampling_conditioned_lpdf(x, par[0], par[1],
    par[2], par[3], k, static_cast<int>(*par[4]), par[5], par[6], par[7], par[8],
    len[8], z_lower.data(), z_upper.data(), bins.data(), static_cast<int>(*par[12]),
    static_cast<int>(*par[13]), *par[14] != 0, static_cast<int>(*par[15]),
    groups.data(), par[17], static_cast<int>(*par[18]), static_cast<int>(*par[19]),
    static_cast<int>(*par[20]), *par[21], &mcse, &change, e.data(), delta.data(),
    &normalizer, quadrature);
  if (!std::isfinite(mcse) || mcse > *par[21]) {
    throw std::runtime_error("Selection normalizer was rejected by diagnostics: relative Monte Carlo standard error was " +
      std::to_string(mcse) + ". Increase 'max_points_per_scramble' or 'scrambles' in 'selection_control'.");
  }
  if (!std::isfinite(change) || change > *par[21]) {
    throw std::runtime_error("Selection normalizer was rejected by diagnostics: relative nested-design change was " +
      std::to_string(change) + ". Increase 'max_points_per_scramble' in 'selection_control'.");
  }
  return std::isfinite(value) ? value : JAGS_NEGINF;
}

void DSELNORMSAMPLINGCONDITIONED::randomSample(double *x, unsigned int length,
    std::vector<double const *> const &par, std::vector<unsigned int> const &len,
    double const *lower, double const *upper, RNG *rng) const
{
  std::fill_n(x, length, std::numeric_limits<double>::quiet_NaN());
}

void DSELNORMSAMPLINGCONDITIONED::typicalValue(double *x, unsigned int length,
    std::vector<double const *> const &par, std::vector<unsigned int> const &len,
    double const *lower, double const *upper) const
{
  std::copy_n(par[0], length, x);
}

unsigned int DSELNORMSAMPLINGCONDITIONED::length(std::vector<unsigned int> const &len) const
{
  return len[0];
}

void DSELNORMSAMPLINGCONDITIONED::support(double *lower, double *upper,
    unsigned int length, std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const
{
  std::fill_n(lower, length, JAGS_NEGINF);
  std::fill_n(upper, length, JAGS_POSINF);
}

bool DSELNORMSAMPLINGCONDITIONED::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}

}
}
