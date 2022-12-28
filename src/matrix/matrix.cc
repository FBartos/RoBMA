#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "matrix.h"

#include <module/ModuleError.h>

/* lapack prototypes */
extern "C"
{
    void dsyev_(const char* jobz, const char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );

    void dgesv_(const int* n, const int* nrhs, double* a,
		            const int* lda,	int* ipiv, double* b, const int* ldb,
		            int* info);

    void dpotrf_(const char *uplo, const int *n, double *a,
		             const int *lda, const int *info);

    void dpotri_(const char *uplo, const int *n, double *a,
		             const int *lda, const int *info);
}

/* adapted from the BUGS module from JAGS*/

namespace jags {
namespace RoBMA {

double logdet(double const *a, int n)
{
  // Log determinant of n x n symmetric positive matrix a */

  int N = n*n;
  double *acopy = new double[N];
  for (int i = 0; i < N; i++) {
    acopy[i] = a[i];
  }

  double *w = new double[n];
  int lwork = -1;
  double worktest = 0;
  int info = 0;
  dsyev_("N","U", &n, acopy, &n, w, &worktest, &lwork, &info);
  if (info != 0) {
    delete [] acopy;
    delete [] w;
    throwRuntimeError("unable to calculate workspace size for dsyev");
  }
  lwork = static_cast<int>(worktest);
  double *work = new double[lwork];
  dsyev_("N","U", &n, acopy, &n, w, work, &lwork, &info);
  delete [] acopy;
  delete [] work;
  if (info != 0) {
    delete [] w;
    throwRuntimeError("unable to calculate eigenvalues in dsyev");
  }

  if (w[0] <= 0) {
      throwRuntimeError("Non positive definite matrix in call to logdet");
  }

  double logdet = 0;
  for (int i = 0; i < n; i++) {
    logdet += std::log(w[i]);
  }
  delete [] w;

  return logdet;
}

bool check_symmetric_ispd(double const *a, int n)
{
  /* Checks that an n x n symmetric matrix is positive definite.
      The code is essentially the same as logdet, but we return
      false if the smallest eigenvalue is less than zero.
  */

  int N = n*n;
  std::vector<double> acopy(N);
  std::copy(a, a+N, acopy.begin());

  //Workspace query to get optimal workspace
  std::vector<double> w(n);
  int lwork = -1;
  double worktest = 0;
  int info = 0;
  dsyev_("N","U", &n, &acopy[0], &n, &w[0], &worktest, &lwork, &info);
  if (info != 0) {
    throwRuntimeError("unable to calculate workspace size for dsyev");
  }
  lwork = static_cast<int>(worktest);
  std::vector<double> work(lwork);

  //Calculate eigenvalues
  dsyev_("N","U", &n, &acopy[0], &n, &w[0], &work[0], &lwork, &info);
  if (info != 0) {
    throwRuntimeError("unable to calculate eigenvalues in dsyev");
  }

  return w[0] > 0;
}


bool inverse_spd (double *X, double const *A, int n)
{
  /* invert n x n symmetric positive definite matrix A. Put result in X*/

  int N = n*n;
  double *Acopy = new double[N];
  for (int i = 0; i < N; i++) {
    Acopy[i] = A[i];
  }

  int info = 0;
  dpotrf_("L", &n, Acopy, &n, &info);
  if (info < 0) {
    throwLogicError("Illegal argument in inverse_spd");
  }
  else if (info > 0) {
    delete [] Acopy;
    throwRuntimeError("Cannot invert matrix: not positive definite");
  }
  dpotri_("L", &n, Acopy, &n, &info);

  for (int i = 0; i < n; ++i) {
    X[i*n + i] = Acopy[i*n + i];
    for (int j = 0; j < i; ++j) {
      X[i*n + j] = X[j*n + i] = Acopy[j*n + i];
    }
  }
  delete [] Acopy;

  if (info != 0) {
    throwRuntimeError("Unable to invert symmetric positive definite matrix");
  }
  return true;
}


bool inverse (double *X, double const *A, int n)
{
  /* invert n x n matrix A. Put result in X */

  int N = n*n;
  double *Acopy = new double[N];
  for (int i = 0; i < N; i++) {
    Acopy[i] = A[i];
    X[i] = 0;
  }
  for (int i = 0; i < n; i++) {
    X[i*n + i] = 1;
  }

  int info = 0;
  int *ipiv = new int[n];
  dgesv_(&n, &n, Acopy, &n, ipiv, X, &n, &info);

  delete [] ipiv;
  delete [] Acopy;

  if (info != 0) {
    return false;
  }
  return true;
}

bool check_symmetry(double const *x, unsigned int n, double tol)
{
  for (unsigned int i = 1; i < n; ++i) {
    double const *xp = x + i;
    double const *yp = x + n*i;
    for (unsigned int j = 0; j < i; ++j) {
      if (std::fabs(*xp - *yp) > tol) return false;
      xp += n;
      yp++;
    }
  }
  return true;
}

}}
