#include "wmnorm.h"

#include "../functions/mnorm.h"
#include "../functions/tools.h"

#include <iostream>
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::vector;
using std::string;
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


namespace jags {
  namespace RoBMA {

    // one-sided multivariate normal
    wmnorm_1s_lpdf::wmnorm_1s_lpdf() :ArrayFunction("wmnorm_1s_lpdf", 5)
    {}
    void wmnorm_1s_lpdf::evaluate(double *value, vector<double const *> const &args, vector<vector<unsigned int> > const &dims) const
    {
      // reassign the addresses to pointers
      const double *x      = args[0];
      const double *mu     = args[1];
      const double *sigma  = args[2];
      const double *crit_x = args[3];
      const double *omega  = args[4];

      // information about the dimensions
      const int K = dims[0][0]; // of the outcome
      const int J = dims[4][0]; // of the weights


      *value = cpp_wmnorm_1s_lpdf(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], K, J);
    }

    bool wmnorm_1s_lpdf::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

		bool wmnorm_1s_lpdf::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

    vector<unsigned int> wmnorm_1s_lpdf::dim(vector<vector<unsigned int> > const &dims,	vector<double const *> const &values) const
    {
	    return vector<unsigned int>(1,1);
    }


    // two-sided multivariate normal
    wmnorm_2s_lpdf::wmnorm_2s_lpdf() :ArrayFunction("wmnorm_2s_lpdf", 5)
    {}
    void wmnorm_2s_lpdf::evaluate(double *value, vector<double const *> const &args, vector<vector<unsigned int> > const &dims) const
    {
      // reassign the addresses to pointers
      const double *x      = args[0];
      const double *mu     = args[1];
      const double *sigma  = args[2];
      const double *crit_x = args[3];
      const double *omega  = args[4];

      // information about the dimensions
      const int K = dims[0][0]; // of the outcome
      const int J = dims[4][0]; // of the weights


      *value = cpp_wmnorm_2s_lpdf(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], K, J);
    }

    bool wmnorm_2s_lpdf::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

		bool wmnorm_2s_lpdf::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

    vector<unsigned int> wmnorm_2s_lpdf::dim(vector<vector<unsigned int> > const &dims,	vector<double const *> const &values) const
    {
	    return vector<unsigned int>(1,1);
    }

  }
}
