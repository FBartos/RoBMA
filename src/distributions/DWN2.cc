#include "DWN2.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>

using std::vector;
using std::log;
using std::exp;
using std::fabs;


// define parameters
// mu  = par[0]
// var = 1/par[1]
#define crit_x(par) (par[2])
#define omega(par) (par[3])
// and their dimensions
#define n_crit_x(len) (len[2])
#define n_omega(len) (len[3])


namespace jags {
namespace weightedt {

DWN2::DWN2() : VectorDist("dwnorm_2s", 4) {}


bool DWN2::checkParameterLength(vector<unsigned int> const &len) const
{
  // there is one less cut-point then weights
  return n_crit_x(len) == n_omega(len) - 1;
}

bool DWN2::checkParameterValue(vector<double const *> const &par,
			    vector<unsigned int> const &len) const
{
  bool crit_x_OK = true;
  bool omega_OK  = true;
  bool var_OK;

  // all crit_x must be non-negative
  for(unsigned i = 1; i < n_crit_x(len); ++i){
    crit_x_OK = crit_x_OK && ( crit_x(par)[i] >= 0.0 );
  }

  // all omegas are within [0, 1] and the last omega == 1
  for(unsigned j = 0; j < (n_omega(len)-1); ++j){
    omega_OK = omega_OK && ( omega(par)[j] >= 0.0 ) && ( omega(par)[j] <= 1.0 );
  }
  // numerical imprecission for last omega is not a problem since it's assumed to be 1 later on
  omega_OK = omega_OK && ( fabs(omega(par)[n_omega(len)-1] - 1.0) < 0.001 );

  // var are positive
  var_OK = *par[1] > 0.0;

  return crit_x_OK && omega_OK && var_OK;
}

double DWN2::logDensity(double const *x, unsigned int length, PDFType type,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{

  double abs_x = fabs(*x);
  double mu    = *par[0];
  double var   = 1/ *par[1];

  double w;
  double nom;
  double denom = 0;
  double denom_sum;
  vector<double> denoms;
  double log_lik;



  // select weight to correspond to the current cut-off
  if(abs_x >= crit_x(par)[n_crit_x(len)-1]){
    // using 1 instead of omega(par)[n_omega(len)-1] because of the numerical imprecision
    w = log(1.0);
  }else if(abs_x < crit_x(par)[0]){
    w = log(omega(par)[0]);
  }else{
    for(unsigned i = 1; i < n_omega(len); ++i){
      if( ( abs_x < crit_x(par)[i] ) && ( abs_x >= crit_x(par)[i-1]) ){
        w = log(omega(par)[i]);
        break;
      }
    }
  }

  // compute the nominator
  nom = dnorm(*x, mu, sqrt(var), true) + w;

  // compute the probabilities between cutpoints
  // the first one
  denoms.push_back(pnorm(crit_x(par)[0], mu, sqrt(var), true, false) -  pnorm(-crit_x(par)[0], mu, sqrt(var), true, false));
  if(denoms[0] < 0.0){ // check and correct for possibly negative numbers due to numerical imprecision
    denoms[0] = 0.0;
  }
  denom_sum = denoms[0];
  // the ones in the middle
  if(n_omega(len) > 1){
    for(unsigned j = 1; j < n_omega(len) - 1; ++j){
      denoms.push_back(pnorm(crit_x(par)[j], mu, sqrt(var), true, false) -  pnorm(-crit_x(par)[j], mu, sqrt(var), true, false) - denom_sum);
      if(denoms[j] < 0.0){ // check and correct for possibly negative numbers due to numerical imprecision
        denoms[j] = 0.0;
      }
      denom_sum = denom_sum + denoms[j];
    }
  }
  // the last one
  denoms.push_back(1.0 - denom_sum);
  if(denoms[n_omega(len)-1] < 0.0){ // check and correct for possibly negative numbers due to numerical imprecision
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

void DWN2::randomSample(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DWN2::support(double *lower, double *upper, unsigned int length,
	     vector<double const *> const &par,
	     vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DWN2::length(vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DWN2::typicalValue(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DWN2::isSupportFixed(vector<bool> const &fixmask) const
{
  return true;
}


}
}

