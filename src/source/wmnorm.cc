#include "wmnorm.h"
#include "mnorm.h"
#include "tools.h"
#include <iostream>
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::cout;
using std::endl;

double cpp_wmnorm_1s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int K, const int J)
{

  // obtain product of the weights (on log scale)
  double log_w = 0;
  for(int k = 0; k < K; k++){
    log_w += log_weight_onesided(&x[k], &crit_x[k * (J - 1)], &omega[0], J);
  }

  // the log weighted normal likelihood
  double log_lik = cpp_mnorm_lpdf(&x[0], &mu[0], &sigma[0], K) + log_w;

  // get the standardizing constant
  double log_std_const = log_std_constant_onesided(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], K, J);

  return log_lik - log_std_const;
}

double cpp_wmnorm_2s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int K, const int J)
{
  // obtain product of the weights (on log scale)
  double log_w = 0;
  for(int k = 0; k < K; k++){
    log_w += log_weight_twosided(&x[k], &crit_x[k * (J - 1)], &omega[0], J);
  }

  // the log weighted normal likelihood
  double log_lik = cpp_mnorm_lpdf(&x[0], &mu[0], &sigma[0], K) + log_w;

  // get the standardizing constant
  double log_std_const = log_std_constant_twosided(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], K, J);

  return log_lik - log_std_const;
}
