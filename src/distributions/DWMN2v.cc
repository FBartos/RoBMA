#include "DWMN2v.h"

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cfloat>
#include <cmath>
#include <vector>
#include <array>


#include <JRmath.h>
#include "../functions/mnorm.h"
#include "../functions/wmnorm.h"
#include "../functions/tools.h"

#include <iostream>


using std::vector;
using std::log;
using std::exp;
using std::sqrt;
using std::fabs;
using std::cout;
using std::endl;


namespace jags {
  namespace RoBMA {

    vector<unsigned int> DWMN2v::dim(vector<vector<unsigned int> > const &dims) const
    {
      return vector<unsigned int>(1,dims[0][0]);
    }


    bool DWMN2v::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      bool se2_OK    = true; // check that standard errors squared and mu dimension matches
      bool tau2_OK   = true; // check that tau squared is a single double
      bool rho2_OK   = true; // check that rho squared is a single double
      bool omega_OK  = true; // check that omega and crit_x dimension matches
      bool crit_x_OK = true; // check that crit_x and mu dimension matches

      se2_OK    = dims[0][0] == dims[1][0];
      tau2_OK   = dims[2][0] == 1;
      rho2_OK   = dims[3][0] == 1;
      if(dims[5][0] == 2){
        crit_x_OK = dims[0][0] == dims[4][0];
        omega_OK  = true;
      }else{
        crit_x_OK = dims[0][0] == dims[4][1];
        omega_OK  = dims[5][0] == (dims[5][0] - 1);
      }

      return se2_OK && tau2_OK && rho2_OK && omega_OK && crit_x_OK;
    }


    bool DWMN2v::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
    {
      const double *tau2   = par[2];
      const double *rho2   = par[3];
      const double *omega  = par[5];

      const int J = dims[5][0];

      bool tau2_OK  = true;  // check that tau2 is a positive double
      bool rho2_OK  = true;  // check that rho2 is between 0 and 1
      bool omega_OK = true;  // check that omega is between 0 and 1

      tau2_OK = *tau2 >= 0;
      rho2_OK = *rho2 >= 0 && *rho2 <= 1;

      for(int i = 0; i < J; i++){
        omega_OK = omega_OK && omega[i] >= 0 && omega[i] <= 1;
      }

      return tau2_OK && rho2_OK && omega_OK;
    }


    DWMN2v::DWMN2v():ArrayDist("dwmnorm_2s_v", 7) {}

    double DWMN2v::logDensity(double const *x, unsigned int length, PDFType type, vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims, double const *lower, double const *upper) const
    {
      // reassign the addresses to pointers
      const double *mu     = par[0];
      const double *se2    = par[1];
      const double *tau2   = par[2];
      const double *rho2   = par[3];
      const double *crit_x = par[4];
      const double *omega  = par[5];
      const double *indx   = par[6];

      // information about the dimensions
      const int J = dims[5][0]; // of the weights
      const int I = dims[6][0]; // of the clusters

      // precompute the covariance
      double cov = *tau2 * *rho2;

      double log_lik = 0;

      for(int i = 0; i < I; i++){
        cout << i << endl;
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

        // parameter printing for debuging
        /*
        cout << "i = " << i << endl;
        cout << "K = " << temp_K << endl;
        cout << "J = " << J << endl;
        cout << "indx_start = " << indx_start << endl;

        //cout << "log_lik = " << log_lik << endl;
        cout << "mu = ";
        for(int k = 0; k < temp_K; k++){
          cout << *(temp_mu + k) << "\t";
        }
        cout << endl;

        cout << "omega = ";
        for(int k = 0; k < J; k++){
          cout << *(omega + k) << "\t";
        }
        cout << endl;

        for(int k = 0; k < temp_K; k++){
          cout << "sigma[" << k << ",] = ";
          for(int j = 0; j < temp_K; j++){
            cout << *(temp_sigma + k * temp_K + j) << "\t";
          }
          cout << endl;
        }

        for(int k = 0; k < temp_K; k++){
          cout << "crit_x[" << k << ",] = ";
          for(int j = 0; j < J-1; j++){
            cout << *(temp_crit_x + k * (J-1) + j) << "\t";
          }
          cout << endl;
        }

        cout << endl;
        cout << "computing loglik" << endl;

        cout << "log_lik = " << log_lik << endl;
        cout << endl;
        cout << endl;
        */ 

        // the log weighted normal likelihood
        log_lik += cpp_wmnorm_2s_lpdf(&temp_x[0], &temp_mu[0], &temp_sigma[0], &temp_crit_x[0], &omega[0], temp_K, J);

        // clean the memory
        delete[] temp_x;
        delete[] temp_mu;
        delete[] temp_sigma;
        delete[] temp_crit_x;
      }

      return log_lik;
    }

    void DWMN2v::randomSample(double *x, unsigned int length, vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims,
              double const *lower, double const *upper,
              RNG *rng) const
    {
      // not implemented
    }

    void DWMN2v::support(double *lower, double *upper, unsigned int length,
              vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims) const
    {
      // no idea whether this is correct
      for (unsigned int i = 0; i < length; ++i) {
        lower[i] = JAGS_NEGINF;
        upper[i] = JAGS_POSINF;
      }
    }

    void DWMN2v::typicalValue(double *x, unsigned int length,
              vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims,
              double const *lower, double const *upper) const
    {
      // not implemented
    }

    bool DWMN2v::isSupportFixed(vector<bool> const &fixmask) const
    {
      return true;
    }


  }
}
