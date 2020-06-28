#include "DTboost.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <boost/math/distributions/non_central_t.hpp>

using std::vector;
using std::log;
using boost::math::non_central_t;
using boost::math::pdf;

namespace jags {
namespace weightedt { 

DTboost::DTboost() : VectorDist("dnt_boost", 3) {}


bool DTboost::checkParameterLength(vector<unsigned int> const &len) const
{
  return true;
}

bool DTboost::checkParameterValue(vector<double const *> const &par,
			    vector<unsigned int> const &len) const
{
  // df are positive
  bool df_OK = *par[2] > 0.0;

  return df_OK;
}

double DTboost::logDensity(double const *x, unsigned int length, PDFType type,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // the sigma parameter is ignored
  double ncp   = *par[0];
  double df    = *par[2];

  // create the boost distribution
  non_central_t t_dist(df, ncp);

  // compute the log likelihood
  double log_lik = log(pdf(t_dist, *x));

  return log_lik;
}

void DTboost::randomSample(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DTboost::support(double *lower, double *upper, unsigned int length,
	     vector<double const *> const &par,
	     vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DTboost::length(vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DTboost::typicalValue(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DTboost::isSupportFixed(vector<bool> const &fixmask) const
{
  return true;
}


}
}

