#ifndef MATRIX_H_
#define MATRIX_H_

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

}}

#endif /* MATRIX_H_ */
