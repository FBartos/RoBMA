#include "wmnorm.h"
#include "../source/wmnorm.h"

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


    //// vector input type ----
    // one-sided multivariate normal
    wmnorm_1s_v_lpdf::wmnorm_1s_v_lpdf() :ArrayFunction("wmnorm_1s_v_lpdf", 8)
    {}
    void wmnorm_1s_v_lpdf::evaluate(double *value, vector<double const *> const &args, vector<vector<unsigned int> > const &dims) const
    {
      // reassign the addresses to pointers
      const double *x      = args[0];
      const double *mu     = args[1];
      const double *se2    = args[2];
      const double *tau2   = args[3];
      const double *rho    = args[4];
      const double *crit_x = args[5];
      const double *omega  = args[6];
      const double *indx   = args[7];

      // information about the dimensions
      const int J = dims[6][0]; // of the weights
      const int I = dims[7][0]; // of the clusters

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
        double * temp_crit_x;

        temp_x      = new double [temp_K];
        temp_mu     = new double [temp_K];
        temp_sigma  = new double [temp_K * temp_K];
        temp_crit_x = new double [temp_K * (J - 1)];

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

          // crit_x
          for(int j = 0; j < J - 1 ; j++){
            *(temp_crit_x + k * (J - 1) + j) = *(crit_x + (indx_start + k) * (J - 1) + j);
          }
        }


        // the log weighted normal likelihood
        log_lik += cpp_wmnorm_1s_lpdf(&temp_x[0], &temp_mu[0], &temp_sigma[0], &temp_crit_x[0], &omega[0], temp_K, J);

        // clean the memory
        delete[] temp_x;
        delete[] temp_mu;
        delete[] temp_sigma;
        delete[] temp_crit_x;
      }

      *value = log_lik;
    }

    bool wmnorm_1s_v_lpdf::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

		bool wmnorm_1s_v_lpdf::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

    vector<unsigned int> wmnorm_1s_v_lpdf::dim(vector<vector<unsigned int> > const &dims,	vector<double const *> const &values) const
    {
	    return vector<unsigned int>(1,1);
    }


    // two-sided multivariate normal
    wmnorm_2s_v_lpdf::wmnorm_2s_v_lpdf() :ArrayFunction("wmnorm_2s_v_lpdf", 8)
    {}
    void wmnorm_2s_v_lpdf::evaluate(double *value, vector<double const *> const &args, vector<vector<unsigned int> > const &dims) const
    {
      // reassign the addresses to pointers
      const double *x      = args[0];
      const double *mu     = args[1];
      const double *se2    = args[2];
      const double *tau2   = args[3];
      const double *rho    = args[4];
      const double *crit_x = args[5];
      const double *omega  = args[6];
      const double *indx   = args[7];

      // information about the dimensions
      const int J = dims[6][0]; // of the weights
      const int I = dims[7][0]; // of the clusters

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
        double * temp_crit_x;

        temp_x      = new double [temp_K];
        temp_mu     = new double [temp_K];
        temp_sigma  = new double [temp_K * temp_K];
        temp_crit_x = new double [temp_K * (J - 1)];

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

          // crit_x
          for(int j = 0; j < J - 1 ; j++){
            *(temp_crit_x + k * (J - 1) + j) = *(crit_x + (indx_start + k) * (J - 1) + j);
          }
        }


        // the log weighted normal likelihood
        log_lik += cpp_wmnorm_1s_lpdf(&temp_x[0], &temp_mu[0], &temp_sigma[0], &temp_crit_x[0], &omega[0], temp_K, J);

        // clean the memory
        delete[] temp_x;
        delete[] temp_mu;
        delete[] temp_sigma;
        delete[] temp_crit_x;
      }

      *value = log_lik;
    }

    bool wmnorm_2s_v_lpdf::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

		bool wmnorm_2s_v_lpdf::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
    {
      return true;
    }

    vector<unsigned int> wmnorm_2s_v_lpdf::dim(vector<vector<unsigned int> > const &dims,	vector<double const *> const &values) const
    {
	    return vector<unsigned int>(1,1);
    }
  }
}
