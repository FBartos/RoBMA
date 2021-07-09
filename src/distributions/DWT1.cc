
#include "DWT1.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>

using std::vector;
using std::log;
using std::exp;
using std::fabs;


// define parameters
// df = par[0]
// ncp = par[1]
#define crit_t(par) (par[2])
#define omega(par) (par[3])
// and their dimensions
#define n_crit_t(len) (len[2])
#define n_omega(len) (len[3])


namespace jags {
namespace RoBMA { 

DWT1::DWT1() : VectorDist("dwt_1s", 4) {}


bool DWT1::checkParameterLength(vector<unsigned int> const &len) const
{
  // there is one less cut-point then weights
  return n_crit_t(len) == n_omega(len) - 1;
}

bool DWT1::checkParameterValue(vector<double const *> const &par,
			    vector<unsigned int> const &len) const
{
  bool omega_OK  = true;
  bool df_OK;

  // all omegas are within [0, 1]
  for(unsigned j = 0; j < (n_omega(len)-1); ++j){
    omega_OK = omega_OK && ( omega(par)[j] >= 0.0 ) && ( omega(par)[j] <= 1.0 );
  }
  
  // df are positive
  df_OK = *par[0] > 0.0;

  return omega_OK && df_OK;
}

double DWT1::logDensity(double const *x, unsigned int length, PDFType type,
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
    w = log(omega(par)[n_omega(len)-1]);
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

  // compute the nominator
  nom = dnt(*x, df, ncp, true) + w;

  // compute the probabilities between cutpoints
  // the first one
  denoms.push_back(pnt(crit_t(par)[0], df, ncp, true, false));
  if(denoms[0] < 0.0){ // check and correct for possibly negative numbers due to numerical imprecission
    denoms[0] = 0.0;
  }
  denom_sum = denoms[0];
  // the ones in the middle
  if(n_omega(len) > 1){
    for(unsigned j = 1; j < n_omega(len) - 1; ++j){
      denoms.push_back(pnt(crit_t(par)[j], df, ncp, true, false) - denom_sum);
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

void DWT1::randomSample(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DWT1::support(double *lower, double *upper, unsigned int length,
	     vector<double const *> const &par,
	     vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DWT1::length(vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DWT1::typicalValue(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DWT1::isSupportFixed(vector<bool> const &fixmask) const
{
  return true;
}


}
}

