
#include "DWT1boost.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include <boost/math/distributions/non_central_t.hpp>

using std::vector;
using std::log;
using std::exp;
using std::fabs;
using boost::math::non_central_t;
using boost::math::pdf;
using boost::math::cdf;

// define parameters
// df = par[0]
// ncp = par[1]
#define crit_t(par) (par[2])
#define omega(par) (par[3])
// and their dimensions
#define n_crit_t(len) (len[2])
#define n_omega(len) (len[3])


namespace jags {
namespace weightedt { 

DWT1boost::DWT1boost() : VectorDist("dwt_1s_boost", 4) {}


bool DWT1boost::checkParameterLength(vector<unsigned int> const &len) const
{
  // there is one less cut-point then weights
  return n_crit_t(len) == n_omega(len) - 1;
}

bool DWT1boost::checkParameterValue(vector<double const *> const &par,
			    vector<unsigned int> const &len) const
{
  bool omega_OK  = true;
  bool df_OK;

  // all omegas are within [0, 1] and the last omega == 1
  for(unsigned j = 0; j < (n_omega(len)-1); ++j){
    omega_OK = omega_OK && ( omega(par)[j] >= 0.0 ) && ( omega(par)[j] <= 1.0 );
  }
  // numerical imprecission for last omega is not a problem since it's assumed to be 1 later on
  omega_OK = omega_OK && ( fabs(omega(par)[n_omega(len)-1] - 1.0) < 0.001 ); 

  // df are positive
  df_OK = *par[0] > 0.0;

  return omega_OK && df_OK;
}

double DWT1boost::logDensity(double const *x, unsigned int length, PDFType type,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  double df    = *par[0];
  double ncp   = *par[1];

  double w;
  double nom;
  double denom = 0;
  double denom_sum;
  vector<double> denoms;
  double log_lik;

  // select weight to correspond to the current cut-off
  if(*x >= crit_t(par)[n_crit_t(len)-1]){
    // using 1 instead of omega(par)[n_omega(len)-1] because of possible numerical imprecission
    w = log(1.0);
  }else if(*x < crit_t(par)[0]){
    w = log(omega(par)[0]);
  }else{
    for(unsigned i = 1; i < n_omega(len); ++i){
      if( (*x < crit_t(par)[i] ) && (*x >= crit_t(par)[i-1]) ){
        w = log(omega(par)[i]);
        break;
      }
    }
  }

  // create the boost distribution object
  non_central_t t_dist(df, ncp);

  // compute the nominator
  nom = log(pdf(t_dist, *x)) + w;

  // compute the probabilities between cutpoints
  // the first one
  denoms.push_back(cdf(t_dist, crit_t(par)[0]));
  if(denoms[0] < 0.0){ // check and correct for possibly negative numbers due to numerical imprecission
    denoms[0] = 0.0;
  }
  denom_sum = denoms[0];
  // the ones in the middle
  if(n_omega(len) > 1){
    for(unsigned j = 1; j < n_omega(len) - 1; ++j){
      denoms.push_back(cdf(t_dist, crit_t(par)[j]) - denom_sum);
      if(denoms[j] < 0.0){ // check and correct for possibly negative numbers due to numerical imprecission
        denoms[j] = 0.0;
      }
      denom_sum = denom_sum + denoms[j];
    }
  }
  // the last one
  denoms.push_back(1.0 - denom_sum);
  if(denoms[n_omega(len)-1] < 0.0){ // check and correct for possibly negative numbers due to numerical imprecission
    denoms[n_omega(len)-1] = 0.0;
  }

  // weight and sum the denominators
  for(unsigned k = 0; k < n_omega(len); ++k){
    denom = denom + exp(log(denoms[k]) + log(omega(par)[k]));
  }

  // compute the log likelihood
  log_lik = nom-log(denom);

  return log_lik;
}

void DWT1boost::randomSample(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DWT1boost::support(double *lower, double *upper, unsigned int length,
	     vector<double const *> const &par,
	     vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DWT1boost::length(vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DWT1boost::typicalValue(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DWT1boost::isSupportFixed(vector<bool> const &fixmask) const
{
  return true;
}


}
}

