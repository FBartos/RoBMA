#include "DWWNMIX.h"

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <JRmath.h>
#include <numeric>

#include "../source/mnorm.h"
#include "../source/wmnorm.h"
#include "../source/tools.h"


namespace jags {
namespace RoBMA {

DWWNMIX::DWWNMIX() : VectorDist("dwwnorm_mix", 7) {}

bool DWWNMIX::checkParameterLength(std::vector<unsigned int> const &len) const
{
  // there is one less cut-point then weights
  return true;
}

bool DWWNMIX::checkParameterValue(std::vector<double const *> const &par,
			    std::vector<unsigned int> const &len) const
{
  return true;
}

// Log Density
double DWWNMIX::logDensity(double const *x, unsigned int length, PDFType type,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
    // extract parameters
    const double *mu     = par[0];
    const double *sigma  = par[1];
    const double *crit_y_all = par[2];
    const double *omega_all  = par[3];

    const double *crit_y_mapping = par[4];
    const double crit_y_mapping_max = static_cast<int>(*par[5]);

    const double weight  = *par[6];

    double log_lik;

    log_lik = cpp_wnorm_mix_lpdf(
      *x, *mu, *sigma, crit_y_all, omega_all, crit_y_mapping,
      static_cast<int>(crit_y_mapping_max)
    ) * weight;

    return log_lik;
}

void DWWNMIX::randomSample(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DWWNMIX::support(double *lower, double *upper, unsigned int length,
	     std::vector<double const *> const &par,
	     std::vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DWWNMIX::length(std::vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DWWNMIX::typicalValue(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DWWNMIX::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}


}}
