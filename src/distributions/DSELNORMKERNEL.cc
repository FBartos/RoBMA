#include "DSELNORMKERNEL.h"

#include <rng/RNG.h>
#include <util/nainf.h>

#include <cmath>
#include "../selnorm/selnorm.h"
#include "selnorm-jags-bounds.h"

namespace jags {
namespace RoBMA {

DSELNORMKERNEL::DSELNORMKERNEL() : VectorDist("dselnorm_kernel", 19) {}

bool DSELNORMKERNEL::checkParameterLength(std::vector<unsigned int> const &len) const
{
  return len[6] == len[7] &&
    len[7] == len[8] &&
    len[13] == 2 &&
    len[14] == 2 &&
    len[16] == len[17] &&
    len[15] == len[16] + 1;
}

bool DSELNORMKERNEL::checkParameterValue(std::vector<double const *> const &par,
                                         std::vector<unsigned int> const &len) const
{
  const int phack_kind = static_cast<int>(*par[12]);
  const int kernel_mode = static_cast<int>(*par[18]);

  return *par[1] > 0 && *par[3] > 0 && *par[4] > 0 && *par[5] > 0 &&
    (*par[10] == 1 || *par[10] == -1) &&
    (*par[11] >= 0 && *par[11] < 1) &&
    (phack_kind == 0 || phack_kind == 1 || phack_kind == 2) &&
    (kernel_mode >= SELKERNEL_NORMAL && kernel_mode <= SELKERNEL_STEP_PHACK_POWER);
}

double DSELNORMKERNEL::logDensity(double const *x, unsigned int length,
                                  PDFType type,
                                  std::vector<double const *> const &par,
                                  std::vector<unsigned int> const &len,
                                  double const *lower, double const *upper) const
{
  const int n_segments = static_cast<int>(len[16]);
  const SelNormJagsBounds z_lower(par[7], len[7]);
  const SelNormJagsBounds z_upper(par[8], len[8]);
  const SelNormJagsBounds segment_bounds(par[15], len[15]);

  SelNormKernelData data;
  data.n_bins               = static_cast<int>(len[7]);
  data.n_segments           = n_segments;
  data.effect_sign          = static_cast<int>(*par[10]);
  data.q                    = 1;
  data.z_lower              = z_lower.data();
  data.z_upper              = z_upper.data();
  data.phack_z_source       = par[13];
  data.phack_z_dest         = par[14];
  data.segment_bounds       = segment_bounds.data();
  data.segment_step_bin     = 0;
  data.segment_phack_region = 0;
  data.segment_step_bin_real = par[16];
  data.segment_phack_region_real = par[17];
  data.trusted_step_partition = false;
  data.telescope_probabilities = false;

  return cpp_selnorm_kernel_lpdf(
    *x,
    *par[0],
    *par[1],
    *par[2],
    *par[3],
    *par[4],
    *par[5],
    par[6],
    static_cast<int>(*par[9]),
    *par[11],
    static_cast<int>(*par[12]),
    static_cast<int>(*par[18]),
    data
  );
}

void DSELNORMKERNEL::randomSample(double *x, unsigned int length,
                                  std::vector<double const *> const &par,
                                  std::vector<unsigned int> const &len,
                                  double const *lower, double const *upper,
                                  RNG *rng) const
{
}

void DSELNORMKERNEL::support(double *lower, double *upper, unsigned int length,
                             std::vector<double const *> const &par,
                             std::vector<unsigned int> const &len) const
{
  for (unsigned int i = 0; i < length; ++i) {
    lower[i] = JAGS_NEGINF;
    upper[i] = JAGS_POSINF;
  }
}

unsigned int DSELNORMKERNEL::length(std::vector<unsigned int> const &len) const
{
  return 1;
}

void DSELNORMKERNEL::typicalValue(double *x, unsigned int length,
                                  std::vector<double const *> const &par,
                                  std::vector<unsigned int> const &len,
                                  double const *lower, double const *upper) const
{
}

bool DSELNORMKERNEL::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}

}
}
