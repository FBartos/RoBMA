
#include "DWB.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>


// define parameters
// p       = par[0]
// N       = par[1]
// weights = par[2]
// and their dimensions

namespace jags {
namespace RoBMA {

DWB::DWB() : VectorDist("dwbinom", 3) {}


bool DWB::checkParameterLength(std::vector<unsigned long> const &len) const
{
  // there is one less cut-point then weights
  return true;
}

bool DWB::checkParameterValue(std::vector<double const *> const &par,
			    std::vector<unsigned long> const &len) const
{
  // p is between 0-1, N is a nonnegative integer, and weight is positive
  bool p_OK      = std::isfinite(*par[0]) && *par[0] >= 0.0 && *par[0] <= 1.0;
  bool N_OK      = std::isfinite(*par[1]) && *par[1] >= 0.0 && std::floor(*par[1]) == *par[1];
  bool weight_OK = std::isfinite(*par[2]) && *par[2] > 0.0;

  return p_OK && N_OK && weight_OK;
}

bool DWB::checkParameterDiscrete(std::vector<bool> const &mask) const
{
  return mask[1];
}

double DWB::logDensity(double const *x, PDFType type,
			  std::vector<double const *> const &par,
			  std::vector<unsigned long> const &len) const
{
  double p = *par[0];
  double N = *par[1];
  double weight  = *par[2];

  if(!std::isfinite(*x) || std::floor(*x) != *x || *x < 0.0 || *x > N){
    return JAGS_NEGINF;
  }

  // compute the nominator
  double log_lik = dbinom(*x, N, p, true) * weight;

  return log_lik;
}

void DWB::randomSample(double *x,
			  std::vector<double const *> const &par,
			  std::vector<unsigned long> const &len,
			  RNG *rng) const
{
  double p = *par[0];
  double N = *par[1];
  x[0] = qbinom(rng->uniform(), N, p, true, false);
}

void DWB::support(double *lower, double *upper,
	     std::vector<double const *> const &par,
	     std::vector<unsigned long> const &len) const
{
  *lower = 0.0;
  *upper = *par[1];
}

unsigned long DWB::length(std::vector<unsigned long> const &len) const
{
  return 1;
}


bool DWB::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return fixmask[1];
}

bool DWB::isDiscreteValued(std::vector<bool> const &mask) const
{
  return true;
}


}
}

