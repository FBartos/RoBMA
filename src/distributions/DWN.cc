
#include "DWN.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>


// define parameters
// mu  = par[0]
// var = 1/par[1]
// weights = par[2]
// and their dimensions
#define n_crit_x(len) (len[2])
#define n_omega(len) (len[3])


namespace jags {
namespace RoBMA {

DWN::DWN() : VectorDist("dwnorm", 3) {}


bool DWN::checkParameterLength(std::vector<unsigned int> const &len) const
{
  // there is one less cut-point then weights
  return true;
}

bool DWN::checkParameterValue(std::vector<double const *> const &par,
			    std::vector<unsigned int> const &len) const
{
  // var and weight is positive
  bool var_OK = *par[1] > 0.0;
  bool weight_OK = *par[2] > 0.0;

  return var_OK && weight_OK;
}

double DWN::logDensity(double const *x, unsigned int length, PDFType type,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  double mu  = *par[0];
  double var = 1/ *par[1];
  double weight  = *par[2];

  // compute the nominator
  double log_lik = dnorm(*x, mu, std::sqrt(var), true) + std::log(weight);

  return log_lik;
}

void DWN::randomSample(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DWN::support(double *lower, double *upper, unsigned int length,
	     std::vector<double const *> const &par,
	     std::vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DWN::length(std::vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DWN::typicalValue(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DWN::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}


}
}

