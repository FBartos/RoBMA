#include "mnorm.h"
#include <mvtnormAPI.h>

// wrapper around the mvtnorm package
double pmnorm(double const *lower, double const *upper, double const *mu, double const *sigma_stdev, double const *sigma_corr, const int K);
{

  // create dynamically allocated arrays for the standardized locations
  double * lower_std;
  double * upper_std;
  double * delta;
  double * infin;

  lower_std = new double [K];
  upper_std = new double [K];
  delta     = new double [K];
  infin     = new double [K];

  // standardized boundary points
  for(int i = 0; i < K; i++){
    lower_std[i] = (lower[i] - mean[i]) / sigma_stdev[i];
    upper_std[i] = (upper[i] - mean[i]) / sigma_stdev[i];
  }

  // return 0 on the same upper and lower bounds
  for(int i = 0; i < K; i++){
    if(fabs(lower_std[i] - upper_std[i]) < EPSILON)
      return 0;
  }

  //


  // mvtnorm settings
  double releps = 0;      // default in mvtnorm: 0
  int maxpts    = 25000;  // default in mvtnorm: 25000
  double abseps = 1e-3;   // default in mvtnorm: 0.001, absolute precision
  int rnd       = 1;      // Get/PutRNGstate
  int nu        = 0;      // degrees of freedom, 0 = normal

  // return values
  double error;
  double value;
  int    inform;

  mvtnorm_C_mvtdst(K,
                   nu,
                   lower_std,
                   upper_std,
                   infin, 
                   sigma_corr,
                   delta,
                   maxpts, abseps, releps,
                   error, value, inform, rnd);

  return value;

                   

  double y = 0;
  return y;
}


