
#include "DWN.h"
#include <util/nainf.h>
#include <stdexcept>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>


// define parameters
// mu  = par[0]
// var = 1/par[1]
// weight = par[2]


namespace jags {
namespace RoBMA {

DWN::DWN() : VectorDist("dwnorm", 3) {}


bool DWN::checkParameterLength(std::vector<unsigned long> const &len) const
{
  return true;
}

bool DWN::checkParameterValue(std::vector<double const *> const &par,
			    std::vector<unsigned long> const &len) const
{
  // var and weight is positive
  bool var_OK = std::isfinite(*par[1]) && *par[1] > 0.0;
  bool weight_OK = std::isfinite(*par[2]) && *par[2] > 0.0;

  return std::isfinite(*par[0]) && var_OK && weight_OK;
}

double DWN::logDensity(double const *x, PDFType type,
			  std::vector<double const *> const &par,
			  std::vector<unsigned long> const &len) const
{
  double mu  = *par[0];
  double var = 1/ *par[1];
  double weight  = *par[2];

  // compute the nominator
  double log_lik = dnorm(*x, mu, std::sqrt(var), true) * weight;

  return log_lik;
}

void DWN::randomSample(double *x,
			  std::vector<double const *> const &par,
			  std::vector<unsigned long> const &len,
			  RNG *rng) const
{
  throw std::runtime_error("dwnorm supports observed likelihood nodes only.");
}

void DWN::support(double *lower, double *upper,
	     std::vector<double const *> const &par,
	     std::vector<unsigned long> const &len) const
{
  *lower = JAGS_NEGINF;
  *upper = JAGS_POSINF;
}

unsigned long DWN::length(std::vector<unsigned long> const &len) const
{
  return 1;
}


bool DWN::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}


}
}

