#include "DSELNORMSTEPSWITCH.h"

#include <rng/RNG.h>
#include <util/nainf.h>
#include <stdexcept>

#include <JRmath.h>
#include <cmath>

#include "../selnorm/selnorm.h"
#include "selnorm-jags-bounds.h"
#include "../samplers/StepRowNormalizer.h"

namespace jags {
namespace RoBMA {

DSELNORMSTEPSWITCH::DSELNORMSTEPSWITCH() : VectorDist("dselnorm_step_switch", 11) {}

bool DSELNORMSTEPSWITCH::checkParameterLength(std::vector<unsigned long> const &len) const
{
  if (len.size() != 11) return false;
  const unsigned int maximum = std::numeric_limits<int>::max();
  for (unsigned long length : len) if (length > maximum) return false;
  for (unsigned int i = 0; i < 4; ++i) if (len[i] != 1) return false;
  return len[4] > 0 &&
    len[4] == len[5] &&
    len[5] == len[6] &&
    len[7] == 1 &&
    len[8] == 1 &&
    len[9] == 1 &&
    len[10] == 1;
}

bool DSELNORMSTEPSWITCH::checkParameterValue(std::vector<double const *> const &par,
                                             std::vector<unsigned long> const &len) const
{
  const int n_bins      = static_cast<int>(len[4]);
  const SelNormJagsBounds z_lower(par[5], len[5]);
  const SelNormJagsBounds z_upper(par[6], len[6]);

  if (!std::isfinite(*par[0]) ||
      !(*par[1] > 0 && *par[2] > 0 && *par[3] > 0) ||
      !(*par[8] == 1 || *par[8] == -1) ||
      !(*par[10] == 0 || *par[10] == 1) ||
      !(*par[9] == SELKERNEL_NORMAL || *par[9] == SELKERNEL_STEP)) {
    return false;
  }

  if (*par[9] == SELKERNEL_NORMAL) {
    return true;
  }

  if (!selnorm_jags_integer_in_range(*par[7], 0, n_bins)) {
    return false;
  }

  for (int b = 0; b < n_bins; ++b) {
    if (!std::isfinite(par[4][b]) || par[4][b] < 0 ||
        !(z_lower[b] < z_upper[b])) {
      return false;
    }
  }

  return selnorm_is_descending_step_partition(
    z_lower.data(), z_upper.data(), n_bins
  );
}

double DSELNORMSTEPSWITCH::logDensity(double const *x,
                                      PDFType type,
                                      std::vector<double const *> const &par,
                                      std::vector<unsigned long> const &len) const
{
  // The one implementation of the row density (StepRowNormalizer.h), which
  // the location block also evaluates; the bounds are decoded only for a row
  // in step mode with an observed bin.
  StepRowLaw law;
  law.sd          = *par[1];
  law.se          = *par[2];
  law.effect_sign = static_cast<int>(*par[8]);
  law.bins        = static_cast<int>(len[4]);
  law.omega       = par[4];
  const int kernel_mode = static_cast<int>(*par[9]);
  const int bin         = static_cast<int>(*par[7]);
  const bool telescope  = static_cast<int>(*par[10]) == 1;
  if (kernel_mode == SELKERNEL_NORMAL || bin == 0) {
    return step_row_log_density(law, *par[3], bin, kernel_mode, telescope, *x, *par[0]);
  }
  const SelNormJagsBounds z_lower(par[5], len[5]);
  const SelNormJagsBounds z_upper(par[6], len[6]);
  law.lower = z_lower.data();
  law.upper = z_upper.data();
  return step_row_log_density(law, *par[3], bin, kernel_mode, telescope, *x, *par[0]);
}

void DSELNORMSTEPSWITCH::randomSample(double *x,
                                      std::vector<double const *> const &par,
                                      std::vector<unsigned long> const &len,
                                      RNG *rng) const
{
  throw std::runtime_error("dselnorm_step_switch supports observed likelihood nodes only.");
}

void DSELNORMSTEPSWITCH::support(double *lower, double *upper,
                                 std::vector<double const *> const &par,
                                 std::vector<unsigned long> const &len) const
{
  *lower = JAGS_NEGINF;
  *upper = JAGS_POSINF;
}

unsigned long DSELNORMSTEPSWITCH::length(std::vector<unsigned long> const &len) const
{
  return 1;
}

bool DSELNORMSTEPSWITCH::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}

}
}
