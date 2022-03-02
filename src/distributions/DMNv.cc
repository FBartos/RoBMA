#include "DMNv.h"

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cfloat>
#include <cmath>
#include <vector>
#include <array>


#include <JRmath.h>
#include "../source/mnorm.h"
#include "../source/tools.h"

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

    vector<unsigned int> DMNv::dim(vector<vector<unsigned int> > const &dims) const
    {
      return vector<unsigned int>(1,dims[0][0]);
    }


    bool DMNv::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      bool se2_OK   = true; // check that standard errors squared and mu dimension matches
      bool tau2_OK  = true; // check that tau squared is a single double
      bool rho_OK   = true; // check that rho squared is a single double

      se2_OK   = dims[0][0] == dims[1][0];
      tau2_OK  = dims[2][0] == 1;
      rho_OK   = dims[3][0] == 1;

      return se2_OK && tau2_OK && rho_OK;
    }


    bool DMNv::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
    {
      const double *tau2  = par[2];
      const double *rho   = par[3];

      bool tau2_OK = true;  // check that tau2 is a positive double
      bool rho_OK  = true;  // check that rho is between 0 and 1

      tau2_OK = *tau2 >= 0;
      rho_OK = *rho >= 0 && *rho <= 1;

      return tau2_OK && rho_OK;
    }


    DMNv::DMNv():ArrayDist("dmnorm_v", 5) {}

    double DMNv::logDensity(double const *x, unsigned int length, PDFType type, vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims, double const *lower, double const *upper) const
    {
      // reassign the addresses to pointers
      const double *mu     = par[0];
      const double *se2    = par[1];
      const double *tau2   = par[2];
      const double *rho    = par[3];
      const double *indx   = par[4];

      // information about the dimensions
      const int I = dims[4][0]; // of the clusters

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
        double * temp_x       = extract_x_v(&x[0], indx_start, temp_K);
        double * temp_mu      = extract_mu_v(&mu[0], indx_start, temp_K);
        double * temp_sigma   = extract_sigma_v(&se2[0], &tau2[0], cov, indx_start, temp_K);

        // the log weighted normal likelihood
        log_lik += cpp_mnorm_lpdf(&temp_x[0], &temp_mu[0], &temp_sigma[0], temp_K);

        // clean the memory
        delete[] temp_x;
        delete[] temp_mu;
        delete[] temp_sigma;
      }

      return log_lik;
    }

    void DMNv::randomSample(double *x, unsigned int length, vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims,
              double const *lower, double const *upper,
              RNG *rng) const
    {
      // not implemented
    }

    void DMNv::support(double *lower, double *upper, unsigned int length,
              vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims) const
    {
      // no idea whether this is correct
      for (unsigned int i = 0; i < length; ++i) {
        lower[i] = JAGS_NEGINF;
        upper[i] = JAGS_POSINF;
      }
    }

    void DMNv::typicalValue(double *x, unsigned int length,
              vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims,
              double const *lower, double const *upper) const
    {
      // not implemented
    }

    bool DMNv::isSupportFixed(vector<bool> const &fixmask) const
    {
      return true;
    }


  }
}
