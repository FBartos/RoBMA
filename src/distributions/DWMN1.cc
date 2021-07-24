#include "DWMN1.h"

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cfloat>
#include <cmath>
#include <vector>

#include <algorithm>
#include <JRmath.h>
#include "../functions/pmwnorm.h"

using std::vector;
using std::log;
using std::exp;
using std::fabs;
using std::copy;

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



vector<unsigned int> DWMN1::dim(vector<vector<unsigned int> > const &dims) const
{
  if (isScalar(dims[0])) {
    return vector<unsigned int>(1,1);
  }
  else {
    return vector<unsigned int> (2, dims[0][0]);
  }
}

bool DWMN1::checkParameterDim (vector<vector<unsigned int> > const &dims) const
{
  return true;
}

bool DWMN1::checkParameterValue(vector<double const *> const &par, vector<vector<unsigned int> > const &dims) const
{
  return true;
}

DWMN1::DWMN1():ArrayDist("dwmnorm_1s", 4) {}

double DWMN1::logDensity(double const *x, unsigned int length, PDFType type, vector<double const *> const &par,
			    vector<vector<unsigned int> > const &dims,
			    double const *lower, double const *upper) const
{

  double mu     = *par[0];
  double sigma  = *par[1];
  double omega  = *par[2];
  double crit_x = *par[3];

  double p;

  vector<double> t_lower;
  vector<double> t_upper;

  t_lower.push_back(0);
  t_lower.push_back(0);
  t_upper.push_back(1);
  t_upper.push_back(1);
  p = pmnorm(t_lower, t_upper,
             mu, sigma);

/*
  double loglik = (k - p - 1) * logdet(x, p) / 2;
  for (unsigned int i = 0; i < p; ++i) {
    loglik -= (k+1) * log(df * x[i*p + i] + 1/(A[i]*A[i])) / 2;

  }
*/
  return 0;
/*
  // imported
  int n = 10;




  double correlationMatrix = 1.0;
  int* infin = new int[n];
  double* delta = new double[n];

  for (int i = 0; i < n; ++i) {
    infin[i] = 0; // (-inf, bound]
 //   lower[i] = 0.0;
    delta[i] = 0.0;
  }


  int n =      // number of observations
  int nu = 0;  // degrees of freedom, 0 = normal

  // mvtnorm settings
  double releps = 0;      // default in mvtnorm: 0
  int maxpts = 25000;     // default in mvtnorm: 25000
  double abseps = 1e-3;   // default in mvtnorm: 0.001, absolute precision
  int rnd = 1;            // Get/PutRNGstate
  int nu = 0; // coresponds to multivariate normal

  // return values
  double error;
  double value;
  int inform;

  mvtnorm_C_mvtdst(&n,
                   &nu,
                   lower,
                   upper,
                   infin, correlationMatrix, delta,
                   &maxpts, &abseps, &releps,
                   &error, &value, &inform, &rnd);
                   delete[] (upper);
                   delete[] (lower);
                   delete[] (infin);
                   delete[] (delta);

                   return value;
                   */
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
