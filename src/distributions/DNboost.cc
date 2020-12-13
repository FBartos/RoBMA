#include "DNboost.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <boost/math/distributions/normal.hpp>

using std::vector;
using std::log;
using boost::math::normal;
using boost::math::pdf;

// define parameters
// mu  = 1/par[0]
// var = par[1]

namespace jags {
namespace weightedt { 

DNboost::DNboost() : VectorDist("dnorm_boost", 3) {}


bool DNboost::checkParameterLength(vector<unsigned int> const &len) const
{
  return true;
}

bool DNboost::checkParameterValue(vector<double const *> const &par,
			    vector<unsigned int> const &len) const
{
  // df are positive
  bool var_OK = *par[1] > 0.0;

  return var_OK;
}

double DNboost::logDensity(double const *x, unsigned int length, PDFType type,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // the sigma parameter is ignored
  double mu  = *par[0];
  double var = 1/ *par[1];

  // create the boost distribution
  normal n_dist(mu, sqrt(var));

  // compute the log likelihood
  double log_lik = log(pdf(n_dist, *x));

  return log_lik;
}

void DNboost::randomSample(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DNboost::support(double *lower, double *upper, unsigned int length,
	     vector<double const *> const &par,
	     vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DNboost::length(vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DNboost::typicalValue(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DNboost::isSupportFixed(vector<bool> const &fixmask) const
{
  return true;
}


}
}

