#ifndef MATRIX_H_
#define MATRIX_H_

#include <string>
#include <rng/RNG.h>

namespace jags {
namespace RoBMA {

/**
 * Inverts a general square matrix using the LAPACK routine DGESV
 *
 * @param X Pointer to an array of length n squared, which will contain
 * the inverse on exit.
 *
 * @param A pointer to array containing the values of the matrix
 *
 * @param n number or rows or columns in the matrix
 */
bool inverse (double *X, double const *A, int n);

/**
 * Inverts a symmetrix positive definite matrix by Cholesky
 * decomposition using the LAPACK routines DPOTRF and DPOTRI.
 *
 * @param X Pointer to an array of length n squared, which will contain
 * the inverse on exit.
 *
 * @param A pointer to array containing the values of the matrix. Only
 * the lower triangle of the matrix (in column-major order) is used.
 *
 * @param n number or rows or columns in the matrix
 */
bool inverse_spd (double *X, double const *A, int n);

/**
 * Checks whether a symmetric matrix is positive definite
 *
 * @param A pointer to array containing the values of the matrix. Only
 * the lower triangle (in column-major order) is used.
 *
 * @param n number or rows or columns in the matrix
 */
bool check_symmetric_ispd(double const *a, int n);

/**
 * Log determinant of a symmetric positive definite matrix
 *
 * @param A pointer to array containing the values of the matrix. Only
 * the lower triangle (in column-major order) is used.
 *
 * @param n number or rows or columns in the matrix
 */
double logdet(double const *A, int n);

/**
 * Checks the symmetry of a square matrix
 *
 * @param A pointer to array containing the values of the matrix
 *
 * @param n number or rows or columns in the matrix
 *
 * @param tol tolerance for symmetry test
 */
bool check_symmetry(double const *X, unsigned int n, double tol=1e-7);

/**
 * Computes eigenvectors of a symmetric matrix using the LAPACK routine DSYEV
 *
 * @param vectors Pointer to an array of length n squared, which will contain
 * the eigenvectors on exit. The eigenvectors are stored column-wise.
 *
 * @param A Pointer to array containing the values of the matrix. Only
 * the upper triangle of the matrix (in column-major order) is used.
 *
 * @param n Number of rows or columns in the matrix
 *
 * @return True if the computation is successful, false otherwise
 */
bool compute_eigenvectors(double *vectors, double const *A, int n);

/**
 * Computes the Cholesky decomposition of a symmetric positive definite matrix.
 *
 * @param L Pointer to an array of length n squared, which will contain
 * the Cholesky factor on exit. Only the lower triangle is used.
 *
 * @param A Pointer to array containing the values of the matrix.
 * Only the lower triangle of the matrix (in column-major order) is used.
 *
 * @param n Number of rows or columns in the matrix
 *
 * @return True if the decomposition succeeded, false otherwise.
 */
bool cholesky_decomposition(double *U, double const *A, int n);

/**
 * Performs the multiplication of a triangular matrix with a vector.
 *
 * Computes x := A * x or x := A^T * x, where A is a triangular matrix.
 *
 * @param uplo 'U' if A is an upper triangular matrix, 'L' if lower triangular.
 * @param trans 'N' for no transpose, 'T' for transpose.
 * @param diag 'N' if A is non-unit triangular, 'U' if unit triangular.
 * @param n The order of the matrix A (number of rows or columns).
 * @param A Pointer to the matrix A, stored in column-major order.
 * @param x Pointer to the vector x, which will be overwritten with the result.
 */
void triangular_matrix_vector_multiply(const char uplo, const char trans, const char diag,
                                       int n, const double *A, double *x);

/**
 * Prints a matrix stored in column-major order.
 *
 * @param matrix Pointer to the matrix data.
 * @param n Size of the matrix (n x n).
 * @param name Name of the matrix to be printed.
 */
void print_matrix(const double *matrix, int n, const std::string &name);

/**
 * Prints a vector.
 *
 * @param vector Pointer to the vector data.
 * @param n Size of the vector.
 * @param name Name of the vector to be printed.
 */
void print_vector(const double *vector, int n, const std::string &name);

bool check_upper_triangular(const double *mat, int K);

void simulate_mnorm_chol(double *x, const double *mu, const double *chol, int K, RNG *rng);

}}

#endif /* MATRIX_H_ */
