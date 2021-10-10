#include "DMN.h"

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

    vector<unsigned int> DMN::dim(vector<vector<unsigned int> > const &dims) const
    {
      return vector<unsigned int>(1,dims[0][0]);
    }

    bool DMN::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      bool sigma_OK  = true; // check that sigma and mu dimension matches

      sigma_OK  = dims[0][0] == dims[1][0] && dims[1][0] == dims[1][1];

      return sigma_OK;
    }

    bool DMN::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
    {
      const double *sigma  = par[1];

      const int K = dims[0][0];

      bool sigma_OK = true;  // check that sigma is symmetric and positive (not that it is semidefinite)

      for(int i = 0; i < K; i++){
        for(int j = 0; j < K && j <= i; j++){
          sigma_OK = sigma_OK && sigma[K * i + j] == sigma[K * j + i] && sigma[K * i + j] >= 0;
        }
      }

      return sigma_OK;
    }


    DMN::DMN():ArrayDist("dmnorm", 2) {}

    double DMN::logDensity(double const *x, unsigned int length, PDFType type, vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims, double const *lower, double const *upper) const
    {
      // reassign the addresses to pointers
      const double *mu     = par[0];
      const double *sigma  = par[1];

      // information about the dimensions
      const int K = dims[0][0]; // of the outcome

      double log_lik = cpp_mnorm_lpdf(&x[0], &mu[0], &sigma[0], K);

      return log_lik;
    }

    void DMN::randomSample(double *x, unsigned int length, vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims,
              double const *lower, double const *upper,
              RNG *rng) const
    {
      // not implemented
    }

    void DMN::support(double *lower, double *upper, unsigned int length,
              vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims) const
    {
      // no idea whether this is correct
      for (unsigned int i = 0; i < length; ++i) {
        lower[i] = JAGS_NEGINF;
        upper[i] = JAGS_POSINF;
      }
    }

    void DMN::typicalValue(double *x, unsigned int length,
              vector<double const *> const &par,
              vector<vector<unsigned int> > const &dims,
              double const *lower, double const *upper) const
    {
      // not implemented
    }

    bool DMN::isSupportFixed(vector<bool> const &fixmask) const
    {
      return true;
    }


  }
}
