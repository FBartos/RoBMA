#include "DWWN2.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>







// define parameters
// mu  = par[0]
// var = 1/par[1]
#define crit_x(par) (par[2])
#define omega(par) (par[3])
// weights = par[4]
// and their dimensions
#define n_crit_x(len) (len[2])
#define n_omega(len) (len[3])


namespace jags {
namespace RoBMA {

DWWN2::DWWN2() : VectorDist("dwwnorm_2s", 5) {}


bool DWWN2::checkParameterLength(std::vector<unsigned int> const &len) const
{
  // there is one less cut-point then weights
  return n_crit_x(len) == n_omega(len) - 1;
}

bool DWWN2::checkParameterValue(std::vector<double const *> const &par,
			    std::vector<unsigned int> const &len) const
{
  bool crit_x_OK = true;
  bool omega_OK  = true;

  // all crit_x must be non-negative
  for(unsigned i = 1; i < n_crit_x(len); ++i){
    crit_x_OK = crit_x_OK && ( crit_x(par)[i] >= 0.0 );
  }

  // all omegas are within [0, 1]
  for(unsigned j = 0; j < (n_omega(len)-1); ++j){
    omega_OK = omega_OK && ( omega(par)[j] >= 0.0 ) && ( omega(par)[j] <= 1.0 );
  }

  // var and weight is positive
  bool var_OK = *par[1] > 0.0;
  bool weight_OK = *par[4] > 0.0;

  return crit_x_OK && omega_OK && var_OK && weight_OK;
}

double DWWN2::logDensity(double const *x, unsigned int length, PDFType type,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{

  double abs_x = std::fabs(*x);
  double mu    = *par[0];
  double var   = 1/ *par[1];
  double weight = *par[4];

  double denom = 0;
  double w = -68;
  std::vector<double> denoms;

  // select weight to correspond to the current cut-off
  if(abs_x >= crit_x(par)[n_crit_x(len)-1]){
    w = std::log(omega(par)[n_omega(len)-1]);
  }else if(abs_x < crit_x(par)[0]){
    w = std::log(omega(par)[0]);
  }else{
    for(unsigned i = 1; i < n_omega(len); ++i){
      if( ( abs_x < crit_x(par)[i] ) && ( abs_x >= crit_x(par)[i-1]) ){
        w = std::log(omega(par)[i]);
        break;
      }
    }
  }

  // compute the nominator
  double nom = dnorm(*x, mu, std::sqrt(var), true) + w;

  // compute the probabilities between cutpoints
  // the first one
  denoms.push_back(pnorm(crit_x(par)[0], mu, std::sqrt(var), true, false) -  pnorm(-crit_x(par)[0], mu, std::sqrt(var), true, false));
  // check and correct for possibly negative numbers due to numerical imprecision
  if(denoms[0] < 0.0){
    denoms[0] = 0.0;
  }
  double denom_sum = denoms[0];
  // the ones in the middle
  if(n_omega(len) > 1){
    for(unsigned j = 1; j < n_omega(len) - 1; ++j){
      denoms.push_back(pnorm(crit_x(par)[j], mu, std::sqrt(var), true, false) -  pnorm(-crit_x(par)[j], mu, std::sqrt(var), true, false) - denom_sum);
      // check and correct for possibly negative numbers due to numerical imprecision
      if(denoms[j] < 0.0){
        denoms[j] = 0.0;
      }
      denom_sum = denom_sum + denoms[j];
    }
  }
  // the last one
  denoms.push_back(1.0 - denom_sum);
  // check and correct for possibly negative numbers due to numerical imprecision
  if(denoms[n_omega(len)-1] < 0.0){
    denoms[n_omega(len)-1] = 0.0;
  }

  // weight and sum the denominators
  for(unsigned k = 0; k < n_omega(len); ++k){
    denom = denom + std::exp(std::log(denoms[k]) + std::log(omega(par)[k]));
  }

  // compute the log likelihood
  double log_lik = (nom-std::log(denom)) * weight;

  return log_lik;
}

void DWWN2::randomSample(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  // not implemented
}

void DWWN2::support(double *lower, double *upper, unsigned int length,
	     std::vector<double const *> const &par,
	     std::vector<unsigned int> const &len) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

unsigned int DWWN2::length(std::vector<unsigned int> const &len) const
{
  // no idea how this works
  return 1;
}


void DWWN2::typicalValue(double *x, unsigned int length,
			  std::vector<double const *> const &par,
			  std::vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
  // not implemented
}


bool DWWN2::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}


}
}

