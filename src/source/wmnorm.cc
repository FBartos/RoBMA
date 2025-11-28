#include "wmnorm.h"
#include "mnorm.h"
#include "tools.h"

#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

double cpp_wnorm_1s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int J)
{

  // obtain the weight (on log scale)
  double log_w = log_weight_onesided(&x[0], &crit_x[0], &omega[0], J);;

  // the log weighted normal likelihood
  double log_lik = dnorm(*x, *mu, *sigma, true) + log_w;

  // get the standardizing constant
  double log_std_const = log_std_constant_onesided(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], J);

  return log_lik - log_std_const;
}

double cpp_wnorm_2s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int J)
{

  // obtain the weight (on log scale)
  double log_w = log_weight_twosided(&x[0], &crit_x[0], &omega[0], J);;

  // the log weighted normal likelihood
  double log_lik = dnorm(*x, *mu, *sigma, true) + log_w;

  // get the standardizing constant
  double log_std_const = log_std_constant_twosided(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], J);

  return log_lik - log_std_const;
}

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
  double log_std_const = log_std_m_constant_onesided(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], K, J);

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
  double log_std_const = log_std_m_constant_twosided(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], K, J);

  return log_lik - log_std_const;
}


// Modified Log-PDF for Hierarchical Selection Model
double cpp_wnorm_hierarchical_lpdf(double const *x, 
                                   double const *theta,   // Study-specific True Effect (theta_i)
                                   double const *se,      // Study Standard Error (sigma_i)
                                   double const *pop_mu,  // Population Mean (mu)
                                   double const *pop_tau, // Population Heterogeneity (tau)
                                   double const *crit_x, 
                                   double const *omega, 
                                   const int J)
{
  // 1. Numerator: Conditional Likelihood
  // log f(y | theta_i) + log w(y)
  // This term fits the data 'x' to the specific latent effect 'theta'
  double log_w = log_weight_onesided(x, crit_x, omega, J);
  double log_lik = dnorm(*x, *theta, *se, true) + log_w;

  // 2. Denominator: Marginal Normalization
  // We normalize by the probability that a study exists in the population.
  // This depends on pop_mu and the total marginal variance (tau^2 + se^2)
  
  double marg_sd_val = std::sqrt(std::pow(*pop_tau, 2) + std::pow(*se, 2));
  
  // Reuse your existing helper, but pass the POPULATION parameters
  // Note: We pass the address of our local variable 'marg_sd_val'
  double log_std_const = log_std_constant_onesided(x, pop_mu, &marg_sd_val, crit_x, omega, J);

  // 3. Return Corrected Log-Likelihood
  return log_lik - log_std_const;
}