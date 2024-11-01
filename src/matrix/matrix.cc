#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "matrix.h"
#include <iostream>   // For std::cout
#include <iomanip>    // For std::setw, std::setprecision
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

bool inverse_spd(double *X, double const *A, int n)
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

bool inverse(double *X, double const *A, int n)
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

bool compute_eigenvectors(double *vectors, double const *A, int n)
{
  /* Computes eigenvectors of an n x n symmetric matrix A. The eigenvectors are stored in 'vectors'. */

  int N = n * n;
  std::vector<double> Acopy(N);
  std::copy(A, A + N, Acopy.begin());

  // Workspace query to get optimal workspace size
  int lwork = -1;
  double work_query;
  int info = 0;
  std::vector<double> w(n); // Eigenvalues (can be used if needed)

  dsyev_("V", "U", &n, &Acopy[0], &n, &w[0], &work_query, &lwork, &info);
  if (info != 0) {
    throwRuntimeError("Unable to calculate workspace size for dsyev");
  }
  lwork = static_cast<int>(work_query);
  std::vector<double> work(lwork);

  // Compute eigenvalues and eigenvectors
  dsyev_("V", "U", &n, &Acopy[0], &n, &w[0], &work[0], &lwork, &info);
  if (info != 0) {
    throwRuntimeError("Unable to compute eigenvectors in dsyev");
  }

  // The eigenvectors are stored in columns of Acopy
  std::copy(Acopy.begin(), Acopy.end(), vectors);

  return true;
}

bool cholesky_decomposition(double *U, double const *A, int n)
{
  // Copy A into U because dpotrf_ overwrites the input matrix
  int N = n * n;
  for (int i = 0; i < N; ++i) {
    U[i] = A[i];
  }

  int info = 0;
  dpotrf_("U", &n, U, &n, &info);
  if (info != 0) {
    return false; // Decomposition failed
  }

  // The upper triangle of U now contains the Cholesky factor
  // The lower triangle is not needed
  return true; // Decomposition succeeded
}

void print_matrix(const double *matrix, int n, const std::string &name)
{
  std::cout << "Matrix " << name << " (" << n << "x" << n << "):\n";
  std::cout << std::fixed << std::setprecision(4);

  for (int i = 0; i < n; ++i) { // Rows
    for (int j = 0; j < n; ++j) { // Columns
      std::cout << std::setw(12) << matrix[i + j * n] << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;

  // Reset format flags if necessary
  std::cout.unsetf(std::ios_base::fixed);
  std::cout.precision(6); // Reset to default precision
}

void print_vector(const double *vector, int n, const std::string &name)
{
  std::cout << "Vector " << name << " (" << n << " elements):\n";
  std::cout << std::fixed << std::setprecision(4);

  for (int i = 0; i < n; ++i) {
    std::cout << std::setw(12) << vector[i] << "\n";
  }
  std::cout << std::endl;

  // Reset format flags if necessary
  std::cout.unsetf(std::ios_base::fixed);
  std::cout.precision(6); // Reset to default precision
}

bool check_upper_triangular(const double *mat, int K)
{
  for (int i = 1; i < K; ++i) {
    for (int j = 0; j < i; ++j) {
      if (std::fabs(mat[i + j * K]) > 1e-8) {
        return false;
      }
    }
  }
  return true;
}

void simulate_mnorm_chol(double *x, const double *mu, const double *chol, int K, RNG *rng) {
  // Generate standard normal variates
  double *z = new double[K];
  for (int i = 0; i < K; ++i) {
    z[i] = rng->normal();
  }

  // Compute x = mu + chol^T * z
  for (int i = 0; i < K; ++i) {
    x[i] = mu[i];
    for (int j = 0; j <= i; ++j) {
      x[i] += chol[j + i * K] * z[j];
    }
  }

  delete[] z;
}


}}
