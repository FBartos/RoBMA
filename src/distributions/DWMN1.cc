#include "DWMN1.h"

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cfloat>
#include <cmath>
#include <vector>
#include <array>

//#include <lapack.h>
//#include <R_ext/Lapack.h>


#include <JRmath.h>
#include "../functions/mnorm.h"
#include "../functions/tools.h"

#include <iostream>


using std::vector;
using std::log;
using std::exp;
using std::sqrt;
using std::fabs;
using std::cout;
using std::endl;




// define parameters
// mu  = par[0]
// var = 1/par[1]
#define crit_x(par) (par[2])
#define omega(par) (par[3])
// and their dimensions
#define n_crit_x(len) (len[2])
#define n_omega(len) (len[3])

#define SCALE(par) (par[0])
#define DF(par)    (*par[1])
#define NROW(dims) (dims[0][0])

namespace jags {
namespace RoBMA {

// identify whether the matrix with crit_y values collapsed to a vector
bool steps_collapsed (vector<vector<unsigned int> > const &dims)
{
  return dims[3].size() == 1;
}


vector<unsigned int> DWMN1::dim(vector<vector<unsigned int> > const &dims) const
{
  return vector<unsigned int>(1,dims[0][0]);
}


bool DWMN1::checkParameterDim (vector<vector<unsigned int> > const &dims) const
{

  bool sigma_OK  = true; // check that sigma and mu dimension matches
  bool omega_OK  = true; // check that omega and crit_x dimension matches
  bool crit_x_OK = true; // check that crit_x and mu dimension matches

  sigma_OK  = dims[0][0] == dims[1][0] && dims[1][0] == dims[1][1];
  crit_x_OK = dims[0][0] == dims[3][0];
  if(steps_collapsed(dims)){
    omega_OK  = dims[2][0] == 2;
  }else{
    omega_OK  = dims[2][0] == (dims[3][1] + 1);
  }

  return sigma_OK && omega_OK && crit_x_OK;
}


bool DWMN1::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
{

  const double *sigma  = par[1];
  const double *omega  = par[2];

  const int K = dims[0][0];
  const int J = dims[2][0];

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


DWMN1::DWMN1():ArrayDist("dwmnorm_1s", 4) {}

double DWMN1::logDensity(double const *x, unsigned int length, PDFType type, vector<double const *> const &par,
			    vector<vector<unsigned int> > const &dims, double const *lower, double const *upper) const
{
  // reassign the addresses to pointers
  const double *mu     = par[0];
  const double *sigma  = par[1];
  const double *omega  = par[2];
  const double *crit_x = par[3];

  // information about the dimensions
  const int K = dims[0][0]; // of the outcome
  const int J = dims[2][0]; // of the weights

  // obtain product of the weights (on log scale)
  double log_w = 0;
  for(int i = 0; i < K; i++){
    log_w += log_weight_onesided(&x[i], &crit_x[i * (J - 1)], &omega[0], J);
  }

  // the log weighted normal likelihood
  double log_lik = dmnorm(&x[0], &mu[0], &sigma[0], K) + log_w;

  // get the standardizing constant
  double log_std_const = log_std_constant_onesided(&x[0], &mu[0], &sigma[0], &crit_x[0], &omega[0], K, J);

  return log_lik - log_std_const;
}

void DWMN1::randomSample(double *x, unsigned int length, vector<double const *> const &par,
			    vector<vector<unsigned int> > const &dims,
			    double const *lower, double const *upper,
			    RNG *rng) const
{
  // not implemented
}

void DWMN1::support(double *lower, double *upper, unsigned int length,
			    vector<double const *> const &par,
			    vector<vector<unsigned int> > const &dims) const
{
  // no idea whether this is correct
  for (unsigned int i = 0; i < length; ++i) {
	  lower[i] = JAGS_NEGINF;
	  upper[i] = JAGS_POSINF;
  }
}

void DWMN1::typicalValue(double *x, unsigned int length,
			    vector<double const *> const &par,
			    vector<vector<unsigned int> > const &dims,
			    double const *lower, double const *upper) const
{
  // not implemented
}

bool DWMN1::isSupportFixed(vector<bool> const &fixmask) const
{
  return true;
}


}}
