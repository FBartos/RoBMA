#include "DWN1.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>

#include "../source/mnorm.h"
#include "../source/wmnorm.h"
#include "../source/tools.h"

// define parameters
// mu  = par[0]
// var = 1/par[1]
#define crit_x(par) (par[2])
#define omega(par) (par[3])
// and their dimensions
#define n_crit_x(len) (len[2])
#define n_omega(len) (len[3])


namespace jags {
namespace RoBMA {

DWN1::DWN1() : VectorDist("dwnorm_1s", 4) {}


bool DWN1::checkParameterLength(std::vector<unsigned int> const &len) const
{
  // there is one less cut-point then weights
  return n_crit_x(len) == n_omega(len) - 1;
}

bool DWN1::checkParameterValue(std::vector<double const *> const &par,
			    std::vector<unsigned int> const &len) const
{
  bool omega_OK  = true;

  // all omegas are within [0, 1]
  for(unsigned j = 0; j < (n_omega(len)-1); ++j){
    omega_OK = omega_OK && ( omega(par)[j] >= 0.0 ) && ( omega(par)[j] <= 1.0 );
  }

  // var is positive
  bool var_OK = *par[1] > 0.0;

  return omega_OK && var_OK;
}

double DWN1::logDensity(double const *x, unsigned int length, PDFType type,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // reassign the addresses to pointers
  const double *mu     = par[0];
  const double var     = 1/ *par[1];
  const double sigma   = std::sqrt(var);
  const double *crit_x = par[2];
  const double *omega  = par[3];

  // information about the dimensions
  const int J = len[3]; // of the weights

  // the log weighted normal likelihood
  double log_lik = cpp_wnorm_1s_lpdf(x, mu, &sigma, crit_x, omega, J);

  return log_lik;
}

void DWN1::randomSample(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DWN1::support(double *lower, double *upper, unsigned int length,
	     std::vector<double const *> const &par,
	     std::vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DWN1::length(std::vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DWN1::typicalValue(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DWN1::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}


}
}

