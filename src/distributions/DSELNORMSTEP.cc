#include "DSELNORMSTEP.h"

#include <rng/RNG.h>
#include <util/nainf.h>

#include <cmath>

#include "../selnorm/selnorm.h"

namespace jags {
namespace RoBMA {

DSELNORMSTEP::DSELNORMSTEP() : VectorDist("dselnorm_step", 10) {}

bool DSELNORMSTEP::checkParameterLength(std::vector<unsigned int> const &len) const
{
  return len[4] > 0 &&
    len[4] == len[5] &&
    len[5] == len[6] &&
    len[7] == 1 &&
    len[8] == 1 &&
    len[9] == 1;
}

bool DSELNORMSTEP::checkParameterValue(std::vector<double const *> const &par,
                                       std::vector<unsigned int> const &len) const
{
  const int n_bins = static_cast<int>(len[4]);
  const int obs_bin = static_cast<int>(*par[7]);
  const int sign = static_cast<int>(*par[8]);
  const int telescope_probabilities = static_cast<int>(*par[9]);

  if (!(*par[1] > 0 && *par[2] > 0 && *par[3] > 0) ||
      !(sign == 1 || sign == -1) ||
      !(telescope_probabilities == 0 || telescope_probabilities == 1) ||
      obs_bin < 1 || obs_bin > n_bins) {
    return false;
  }

  for (int b = 0; b < n_bins; ++b) {
    if (!std::isfinite(par[4][b]) || par[4][b] < 0 ||
        !(par[5][b] < par[6][b])) {
      return false;
    }
  }

  return selnorm_is_descending_step_partition(par[5], par[6], n_bins);
}

double DSELNORMSTEP::logDensity(double const *x, unsigned int length,
                                PDFType type,
                                std::vector<double const *> const &par,
                                std::vector<unsigned int> const &len,
                                double const *lower, double const *upper) const
{
  double phack_z_zero[2] = {0, 0};
  double segment_bounds_zero[1] = {0};
  int segment_zero[1] = {0};

  SelNormKernelData data;
  data.n_bins               = static_cast<int>(len[4]);
  data.n_segments           = 0;
  data.effect_sign          = static_cast<int>(*par[8]);
  data.q                    = 0;
  data.z_lower              = par[5];
  data.z_upper              = par[6];
  data.phack_z_source       = phack_z_zero;
  data.phack_z_dest         = phack_z_zero;
  data.segment_bounds       = segment_bounds_zero;
  data.segment_step_bin     = segment_zero;
  data.segment_phack_region = segment_zero;
  data.segment_step_bin_real = 0;
  data.segment_phack_region_real = 0;
  data.trusted_step_partition = true;
  data.telescope_probabilities = static_cast<int>(*par[9]) == 1;

  return cpp_selnorm_kernel_lpdf(
    *x,
    *par[0],
    *par[1],
    *par[0],
    *par[1],
    *par[2],
    *par[3],
    par[4],
    static_cast<int>(*par[7]),
    0,
    0,
    SELKERNEL_STEP,
    data,
    1,
    false
  );
}

void DSELNORMSTEP::randomSample(double *x, unsigned int length,
                                std::vector<double const *> const &par,
                                std::vector<unsigned int> const &len,
                                double const *lower, double const *upper,
                                RNG *rng) const
{
}

void DSELNORMSTEP::support(double *lower, double *upper, unsigned int length,
                           std::vector<double const *> const &par,
                           std::vector<unsigned int> const &len) const
{
  for (unsigned int i = 0; i < length; ++i) {
    lower[i] = JAGS_NEGINF;
    upper[i] = JAGS_POSINF;
  }
}

unsigned int DSELNORMSTEP::length(std::vector<unsigned int> const &len) const
{
  return 1;
}

void DSELNORMSTEP::typicalValue(double *x, unsigned int length,
                                std::vector<double const *> const &par,
                                std::vector<unsigned int> const &len,
                                double const *lower, double const *upper) const
{
}

bool DSELNORMSTEP::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}

}
}
