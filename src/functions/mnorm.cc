#include "mnorm.h"
#include "../source/mnorm.h"

#include <iostream>
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::vector;
using std::cout;
using std::endl;


namespace jags {
  namespace RoBMA {

    //// normal input type ----
    // one-sided multivariate normal
    mnorm_lpdf::mnorm_lpdf() :ArrayFunction("mnorm_lpdf", 3)
    {}
    void mnorm_lpdf::evaluate(double *value, vector<double const *> const &args, vector<vector<unsigned int> > const &dims) const
    {
      // reassign the addresses to pointers
      const double *x      = args[0];
      const double *mu     = args[1];
      const double *sigma  = args[2];

      // information about the dimensions
      const int K = dims[0][0]; // of the outcome


      *value = cpp_mnorm_lpdf(&x[0], &mu[0], &sigma[0], K);;
    }

    bool mnorm_lpdf::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

		bool mnorm_lpdf::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

    vector<unsigned int> mnorm_lpdf::dim(vector<vector<unsigned int> > const &dims,	vector<double const *> const &values) const
    {
	    return vector<unsigned int>(1,1);
    }


    //// vector input type ----
    mnorm_v_lpdf::mnorm_v_lpdf() :ArrayFunction("mnorm_v_lpdf", 6)
    {}
    void mnorm_v_lpdf::evaluate(double *value, vector<double const *> const &args, vector<vector<unsigned int> > const &dims) const
    {
      // reassign the addresses to pointers
      const double *x      = args[0];
      const double *mu     = args[1];
      const double *se2    = args[2];
      const double *tau2   = args[3];
      const double *rho    = args[4];
      const double *indx   = args[5];

      // information about the dimensions
      const int I = dims[5][0]; // of the clusters

     // precompute the covariance
      double cov = *tau2 * *rho;

      double log_lik = 0;

      for(int i = 0; i < I; i++){
        int temp_K;
        if(i == 0){
          temp_K = *(indx + i);
        }else{
          temp_K = *(indx + i) - *(indx + i - 1);
        }
        int indx_start = *(indx + i) - temp_K;

        // construct the covariance matricies and mean vectors for each set of observations
        double * temp_x;
        double * temp_mu;
        double * temp_sigma;

        temp_x      = new double [temp_K];
        temp_mu     = new double [temp_K];
        temp_sigma  = new double [temp_K * temp_K];

        for(int k = 0; k < temp_K; k++){

          // observations and means
          *(temp_x  + k) = *(x  + indx_start + k);
          *(temp_mu + k) = *(mu + indx_start + k);

          // sigma
          for(int j = 0; j < temp_K; j++){
            if(k == j){
              *(temp_sigma + k * temp_K + j) = *(se2 + indx_start + k) + *tau2;
            }else{
              *(temp_sigma + k * temp_K + j) = cov;
            }
          }
        }


        // the log weighted normal likelihood
        log_lik += cpp_mnorm_lpdf(&temp_x[0], &temp_mu[0], &temp_sigma[0], temp_K);

        // clean the memory
        delete[] temp_x;
        delete[] temp_mu;
        delete[] temp_sigma;
      }

      *value = log_lik;
    }

    bool mnorm_v_lpdf::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

		bool mnorm_v_lpdf::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

    vector<unsigned int> mnorm_v_lpdf::dim(vector<vector<unsigned int> > const &dims,	vector<double const *> const &values) const
    {
	    return vector<unsigned int>(1,1);
    }

  }
}
