#include "DWMN2v.h"

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cfloat>
#include <cmath>
#include <vector>
#include <array>

#include <JRmath.h>
#include "../source/mnorm.h"
#include "../source/wmnorm.h"
#include "../source/tools.h"

namespace jags {
  namespace RoBMA {

    std::vector<unsigned int> DWMN2v::dim(std::vector<std::vector<unsigned int> > const &dims) const
    {
      return std::vector<unsigned int>(1,dims[0][0]);
    }


    bool DWMN2v::checkParameterDim (std::vector<std::vector<unsigned int> > const &dims) const
    {
      bool se2_OK    = true; // check that standard errors squared and mu dimension matches
      bool tau2_OK   = true; // check that tau squared is a single double
      bool rho_OK    = true; // check that rho squared is a single double
      bool omega_OK  = true; // check that omega and crit_x dimension matches
      bool crit_x_OK = true; // check that crit_x and mu dimension matches

      se2_OK    = dims[0][0] == dims[1][0];
      tau2_OK   = dims[2][0] == 1;
      rho_OK    = dims[3][0] == 1;
      if(dims[5][0] == 2){
        crit_x_OK = dims[0][0] == dims[4][0];
        omega_OK  = dims[5][0] == 2;
      }else{
        crit_x_OK = dims[0][0] == dims[4][1];
        omega_OK  = dims[5][0] == (dims[4][0] + 1);
      }

      return se2_OK && tau2_OK && rho_OK  && omega_OK && crit_x_OK;
    }


    bool DWMN2v::checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned int> > const &dims) const
    {
      const double *tau2   = par[2];
      const double *rho    = par[3];
      const double *omega  = par[5];

      const int J = dims[5][0];

      bool tau2_OK  = true;  // check that tau2 is a positive double
      bool rho_OK   = true;  // check that rho  is between 0 and 1
      bool omega_OK = true;  // check that omega is between 0 and 1

      tau2_OK = *tau2 >= 0;
      rho_OK  = *rho  >= 0 && *rho <= 1;

      for(int i = 0; i < J; i++){
        omega_OK = omega_OK && omega[i] >= 0 && omega[i] <= 1;
      }

      return tau2_OK && rho_OK  && omega_OK;
    }


    DWMN2v::DWMN2v():ArrayDist("dwmnorm_2s_v", 7) {}

    double DWMN2v::logDensity(double const *x, unsigned int length, PDFType type, std::vector<double const *> const &par,
              std::vector<std::vector<unsigned int> > const &dims, double const *lower, double const *upper) const
    {
      // reassign the addresses to pointers
      const double *mu     = par[0];
      const double *se2    = par[1];
      const double *tau2   = par[2];
      const double *rho    = par[3];
      const double *crit_x = par[4];
      const double *omega  = par[5];
      const double *indx   = par[6];

      // information about the dimensions
      const int J = dims[5][0]; // of the weights
      const int I = dims[6][0]; // of the clusters

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
        double * temp_crit_x  = extract_crit_x_v(&crit_x[0], indx_start, temp_K, J);

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

    void DWMN2v::randomSample(double *x, unsigned int length, std::vector<double const *> const &par,
              std::vector<std::vector<unsigned int> > const &dims,
              double const *lower, double const *upper,
              RNG *rng) const
    {
      // not implemented
    }

    void DWMN2v::support(double *lower, double *upper, unsigned int length,
              std::vector<double const *> const &par,
              std::vector<std::vector<unsigned int> > const &dims) const
    {
      // no idea whether this is correct
      for (unsigned int i = 0; i < length; ++i) {
        lower[i] = JAGS_NEGINF;
        upper[i] = JAGS_POSINF;
      }
    }

    void DWMN2v::typicalValue(double *x, unsigned int length,
              std::vector<double const *> const &par,
              std::vector<std::vector<unsigned int> > const &dims,
              double const *lower, double const *upper) const
    {
      // not implemented
    }

    bool DWMN2v::isSupportFixed(std::vector<bool> const &fixmask) const
    {
      return true;
    }


  }
}
