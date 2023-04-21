
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


bool DWB::checkParameterLength(std::vector<unsigned int> const &len) const
{
  // there is one less cut-point then weights
  return true;
}

bool DWB::checkParameterValue(std::vector<double const *> const &par,
			    std::vector<unsigned int> const &len) const
{
  // p is between 0-1, n and weight are positive
  bool p_OK = *par[0] >= 0.0 && *par[0] <= 1.0 ; 
  bool N_OK = *par[1] >= 0.0; 
  bool weight_OK = *par[2] > 0.0;

  return p_OK && N_OK && weight_OK;
}

double DWB::logDensity(double const *x, unsigned int length, PDFType type,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  double p = *par[0];
  double N = *par[1];
  double weight  = *par[2];

  // compute the nominator
  double log_lik = dbinom(*x, N, p, true) * weight;

  return log_lik;
}

void DWB::randomSample(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DWB::support(double *lower, double *upper, unsigned int length,
	     std::vector<double const *> const &par,
	     std::vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DWB::length(std::vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DWB::typicalValue(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DWB::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}


}
}

