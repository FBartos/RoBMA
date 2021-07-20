#include "DWMN1.h"

//#include <mvtnormAPI.h>
/*using mvtnorm::C_mvtdst;
// return values
double error;
double value;
int inform;

mvtnorm_C_mvtdst(&n, &nu, lower, upper_,
infin, correlationMatrix, delta,
&maxpts, &abseps, &releps,
&error, &value, &inform, &rnd);
delete[] (lower);
delete[] (infin);
delete[] (delta);

return value;


#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>

using std::vector;
using std::log;
using std::exp;
using std::fabs;
*/

// from wishard
//#include <config.h>

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
//#include <module/ModuleError.h>

#include "DWMN1.h"

#include <cfloat>
#include <cmath>
#include <vector>
#include <algorithm>

#include <JRmath.h>
#include "mvtnorm_workaround.h"



using std::vector;
using std::log;
using std::copy;

#define SCALE(par) (par[0])
#define DF(par)    (*par[1])
#define NROW(dims)  (dims[0][0])
/*
#define F77_DSYEV F77_FUNC(dsyev,DSYEV)

extern "C" {
void F77_DSYEV (const char* jobz, const char* uplo,
const int* n, double* a, const int* lda,
double* w, 
double* work, const int* lwork, int* info);
}
*/
namespace jags {
namespace RoBMA {

//FIXME: This is copy-pasted from the bugs module
static double logdet(double const *a, int n)
{
// Log determinant of n x n symmetric positive matrix a */

  int N = n*n;
  vector<double> acopy(N);
  copy(a, a+N, acopy.begin());

  vector<double> w(n);
  int lwork = -1;
  double worktest = 0;
  int info = 0;
  //F77_DSYEV("N","U", &n, &acopy[0], &n, &w[0], &worktest, &lwork, &info);
  //	if (info != 0) {
  //	    throwRuntimeError("unable to calculate workspace size for dsyev");
  //	}
  lwork = static_cast<int>(worktest);
  double *work = new double[lwork];
  //F77_DSYEV("N","U", &n, &acopy[0], &n, &w[0], work, &lwork, &info);
  delete [] work;
  if (info != 0) {
  //throwRuntimeError("unable to calculate eigenvalues in dsyev");
  }

  if (w[0] <= 0) {
  //throwRuntimeError("Non positive definite matrix in call to logdet");
  }

  double logdet = 0;
  for (int i = 0; i < n; i++) {
    logdet += log(w[i]);
  }

  return logdet;
}

static double log_multigamma(double n, unsigned int p)
{
  double y =  (p * (p-1) * log(M_PI))/4;
  for (unsigned int j = 0; j < p; ++j) {
    y += lgammafn((n-j)/2);
  }
  return y;
}

DWMN1::DWMN1():ArrayDist("dwmnorm_1s", 2) 
{}

double
DWMN1::logDensity(double const *x, unsigned int length,
PDFType type,
vector<double const *> const &par,
vector<vector<unsigned int> > const &dims,
double const *lower, double const *upper) const
{
  double const *A = SCALE(par);
  unsigned int p = NROW(dims);
  double df = DF(par);
  double k = p + df - 1;

  double loglik = (k - p - 1) * logdet(x, p) / 2;
  for (unsigned int i = 0; i < p; ++i) {
    loglik -= (k+1) * log(df * x[i*p + i] + 1/(A[i]*A[i])) / 2;

  }

  if (type != PDF_PRIOR) {
    //Normalize density
    loglik += p * k * log(df) / 2;
    for (unsigned int i = 0; i < p; ++i) {
    loglik -= log(A[i]);
    }
    loglik += p * lgammafn((k+1)/2) -  p * lgammafn(0.5) -
    - log_multigamma(k, p);
  }

return loglik;
}

void DWMN1::sampleWishart(double *x, unsigned int length,
double const *R, unsigned int nrow,
double k, RNG *rng)
{
/* 
Generate random Wishart variable with diagonal prior matrix
Equivalent to DWish(diag(R), k)

The algorithm is adapted from bugs::DWish::randomSample but
is considerably simpler.
*/

//if (length != nrow * nrow) {
//throwLogicError("invalid length in DWMN1::sampleWishart");
//}

/* Generate square root of Wishart random variable:
- diagonal elements are square root of Chi square
- upper off-diagonal elements are normal
- lower off-diagonal elements are zero
*/
  vector<vector<double> > Z(nrow, vector<double>(nrow));

  for (unsigned int i = 0; i < nrow; ++i) {
    for (unsigned int l = 0; l < i; ++l) {
      Z[i][l] = rnorm(0, 1, rng);
    }
    Z[i][i] = sqrt(rchisq(k - i, rng));    
  }

/* Take inverse square root of R */
  vector<double> scale(nrow);
  for (unsigned int i = 0; i < nrow; ++i) {
    scale[i] = 1/sqrt(R[i]);
  }

/* Now put cross-product into x */
  for (unsigned int i = 0; i < nrow; i++) {
    for (unsigned int j = 0; j <= i; j++) {
      double xx = 0;
      for (unsigned int l = 0; l <= j; l++) {
        xx += Z[i][l] * Z[j][l];
      }
      x[nrow * j + i] = x[nrow * i + j] = scale[j] * scale[i] * xx;
    }
  }
}

void DWMN1::randomSample(double *x, unsigned int length,
vector<double const *> const &par,
vector<vector<unsigned int> > const &dims,
double const *lower, double const *upper,
RNG *rng) const
{
  unsigned int nrow = NROW(dims);
  double df = DF(par);
  double const *scale = SCALE(par);

  vector<double> R(nrow);
  double k = nrow + df - 1;
  for (unsigned int i = 0; i < nrow; ++i) {
  //NB Rmath implementation of gamma is parameterized in terms of
  //shape and scale, unlike BUGS which uses shape and rate
    R[i] = 2 * df * rgamma(0.5, scale[i]*scale[i], rng);
  }
//sampleWishart(x, length, &R[0], nrow, k, rng);
}

bool
DWMN1::checkParameterDim (vector<vector<unsigned int> > const &dims)
const
{
  return (isVector(dims[0]) || isScalar(dims[0])) && isScalar(dims[1]);
}

vector<unsigned int> 
DWMN1::dim(vector<vector<unsigned int> > const &dims) const
{
  if (isScalar(dims[0])) {
    return vector<unsigned int>(1,1);
  }
  else {
    return vector<unsigned int> (2, dims[0][0]);
  }
}

bool 
DWMN1::checkParameterValue(vector<double const *> const &par,
vector<vector<unsigned int> > const &dims) const
{
  if (DF(par) < 1) return false;
  double const *scale = SCALE(par);
  unsigned int n = NROW(dims);
  for (unsigned int i = 0; i < n; ++i) {
    if (scale[i] <= 0) return false;
  }
  return true;
}


void DWMN1::support(double *lower, double *upper, unsigned int length,
vector<double const *> const &par,
vector<vector<unsigned int> > const &dims) const
{
  for (unsigned int i = 0; i < length; ++i) {
    if (i % NROW(dims) == i / NROW(dims)) {
    //Diagonal elements
      lower[i] =  0;
    }
    else {
      lower[i] =  JAGS_NEGINF;
    }
  upper[i] = JAGS_POSINF; 
  }
}

void DWMN1::typicalValue(double *x, unsigned int length,
  vector<double const *> const &par,
  vector<vector<unsigned int> > const &dims,
  double const *lower, double const *upper) const
{
/* Returns the mean as a typical value. */

  for (unsigned int i = 0; i < length; ++i) {
    x[i] = 0;
  }
  for (unsigned int i = 0; i < NROW(dims); ++i) {
    unsigned int k = i * NROW(dims) + i;
    x[k] = DF(par)/(SCALE(par)[i] * SCALE(par)[i]);
  }
}

bool DWMN1::isSupportFixed(vector<bool> const &fixmask) const
{
  return true;
}

unsigned int DWMN1::df(vector<vector<unsigned int> > const &dims) const
{   
  return dims[0][0] * (dims[0][0] + 1) / 2;
}

}}
