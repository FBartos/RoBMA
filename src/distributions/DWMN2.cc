#include "DWMN2.h"

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

    vector<unsigned int> DWMN2::dim(vector<vector<unsigned int> > const &dims) const
    {
      return vector<unsigned int>(1,dims[0][0]);
    }


    bool DWMN2::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      bool sigma_OK  = true; // check that sigma and mu dimension matches
      bool omega_OK  = true; // check that omega and crit_x dimension matches
      bool crit_x_OK = true; // check that crit_x and mu dimension matches

      sigma_OK  = dims[0][0] == dims[1][0] && dims[1][0] == dims[1][1];
      if(dims[3][0] == 2){
        crit_x_OK = dims[0][0] == dims[2][0];
        omega_OK  = true;
      }else{
        crit_x_OK = dims[0][0] == dims[2][1];
        omega_OK  = dims[2][0] == (dims[3][0] - 1);
      }

      return sigma_OK && omega_OK && crit_x_OK;
    }


    bool DWMN2::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
    {
      const double *sigma  = par[1];
      const double *omega  = par[3];

      const int K = dims[0][0];
      const int J = dims[3][0];

      bool sigma_OK = true;  // check that sigma is symmetric and positive (not that it is semidefinite)
      bool omega_OK = true;  // check that omega is between 0 and 1

      for(int i = 0; i < K; i++){
        for(int j = 0; j < K && j <= i; j++){
          sigma_OK = sigma_OK && sigma[K * i + j] == sigma[K * j + i] && sigma[K * i + j] >= 0;
        }
      }

      for(int i = 0; i < J; i++){
        omega_OK = omega_OK && omega[i] >= 0 && omega[i] <= 1;
      }

      return sigma_OK && omega_OK;
    }


    DWMN2::DWMN2():ArrayDist("dwmnorm_2s", 4) {}

    double DWMN2::logDensity(double const *x, unsigned int length, PDFType type, vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims, double const *lower, double const *upper) const
    {
      // reassign the addresses to pointers
      const double *mu     = par[0];
      const double *sigma  = par[1];
      const double *crit_x = par[2];
      const double *omega  = par[3];

      // information about the dimensions
      const int K = dims[0][0]; // of the outcome
      const int J = dims[3][0]; // of the weights

      // the log weighted normal likelihood
      double log_lik = cpp_wmnorm_2s_lpdf(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], K, J);

      return log_lik;
    }

    void DWMN2::randomSample(double *x, unsigned int length, vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims,
              double const *lower, double const *upper,
              RNG *rng) const
    {
      // not implemented
    }

    void DWMN2::support(double *lower, double *upper, unsigned int length,
              vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims) const
    {
      // no idea whether this is correct
      for (unsigned int i = 0; i < length; ++i) {
        lower[i] = JAGS_NEGINF;
        upper[i] = JAGS_POSINF;
      }
    }

    void DWMN2::typicalValue(double *x, unsigned int length,
              vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims,
              double const *lower, double const *upper) const
    {
      // not implemented
    }

    bool DWMN2::isSupportFixed(vector<bool> const &fixmask) const
    {
      return true;
    }


  }
}
