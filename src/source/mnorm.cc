#include "mnorm.h"
#include <mvtnormAPI.h>
#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cfloat>
#include <cmath>
#include <vector>
#include <array>
#include <JRmath.h>
#include "../matrix/matrix.h"

using namespace std;

// wrapper around the mvtnorm package
double cpp_mnorm_cdf(double *lower, double *upper, int *infin, double *mu, double *sigma_stdev, double *sigma_corr, int K)
{
  // create dynamically allocated arrays for the standardized locations
  double * lower_std;
  double * upper_std;
  double * delta;

  lower_std = new double [K];
  upper_std = new double [K];
  delta     = new double [K];

  // standardized boundary points
  for(int k = 0; k < K; k++){
    *(lower_std + k) = ( *(lower + k) - *(mu + k)) / *(sigma_stdev + k);
    *(upper_std + k) = ( *(upper + k) - *(mu + k)) / *(sigma_stdev + k);
    *(delta + k)     = 0;
  }

  // mvtnorm settings
  double releps = 0;      // default in mvtnorm: 0
  int    maxpts = 25000;  // default in mvtnorm: 25000
  double abseps = 1e-3;   // default in mvtnorm: 0.001, absolute precision
  int    rnd    = 1;      // Get/PutRNGstate
  int    nu     = 0;      // degrees of freedom, 0 = normal

  // return values
  double error  = 0;
  double value  = 0;
  int    inform = 0;

  mvtnorm_C_mvtdst(&K,
                   &nu,
                   lower_std,
                   upper_std,
                   infin,
                   sigma_corr,
                   delta,
                   &maxpts, &abseps, &releps,
                   &error, &value, &inform, &rnd);

  // clean the memory
  delete[] lower_std;
  delete[] upper_std;
  delete[] delta;

  return value;
}

// adapted from the BUGS module from JAGS
double cpp_mnorm_lpdf(double const *x, double const *mu, double const *sigma, const int K)
{

  vector<double> T(K * K);
  jags::RoBMA::inverse_spd (&T[0], sigma, K);

  double loglik = 0;
  vector<double> delta(K);
  for (int i = 0; i < K; ++i) {
    delta[i] = x[i] - mu[i];
    loglik -= (delta[i] * T[i + i * K] * delta[i])/2;
    for (int j = 0; j < i; ++j) {
      loglik -= delta[i] * T[i + j * K] * delta[j];
    }
  }

  loglik -= jags::RoBMA::logdet(sigma, K)/2 + K * M_LN_SQRT_2PI;

  return loglik;
}
