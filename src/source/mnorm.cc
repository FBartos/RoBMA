#include "mnorm.h"
#include <mvtnormAPI.h>
#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cfloat>
#include <cmath>
#include <vector>
#include <array>
#include <JRmath.h>

// wrapper around the mvtnorm package
double cpp_mnorm_cdf(double *lower, double *upper, int *infin, double *mu, double *sigma_stdev, double *sigma_corr, int K)
{
  // create dynamically allocated arrays for the standardized locations
  double * lower_std;
  double * upper_std;
  double * delta;

  lower_std = new double [K];
  upper_std = new double [K];
  delta     = new double [K];

  // standardized boundary points
  for(int k = 0; k < K; k++){
    *(lower_std + k) = ( *(lower + k) - *(mu + k)) / *(sigma_stdev + k);
    *(upper_std + k) = ( *(upper + k) - *(mu + k)) / *(sigma_stdev + k);
    *(delta + k)     = 0;
  }

  // mvtnorm settings
  double releps = 0;      // default in mvtnorm: 0
  int    maxpts = 25000;  // default in mvtnorm: 25000
  double abseps = 1e-3;   // default in mvtnorm: 0.001, absolute precision
  int    rnd    = 1;      // Get/PutRNGstate
  int    nu     = 0;      // degrees of freedom, 0 = normal

  // return values
  double error  = 0;
  double value  = 0;
  int    inform = 0;

  mvtnorm_C_mvtdst(&K,
                   &nu,
                   lower_std,
                   upper_std,
                   infin,
                   sigma_corr,
                   delta,
                   &maxpts, &abseps, &releps,
                   &error, &value, &inform, &rnd);

  // clean the memory
  delete[] lower_std;
  delete[] upper_std;
  delete[] delta;

  return value;
}

double cpp_mnorm_lpdf(double const *x, double const *mu, double const *sigma, const int K)
{

  double * sigma_chol;
  double * chol_inv;
  sigma_chol = new double [K*K];
  chol_inv   = new double [K*K];

  for(int i = 0; i < K*K; i++){
    *(sigma_chol + i) = 0;
    *(chol_inv + i)   = 0;
  }

  // lower triangle cholesky decomposition
  chol(&sigma[0], K, sigma_chol);

  // inverse
  inverse(&sigma_chol[0], K, chol_inv);

  // product of the diagonal ellements
  double diag_prod = 0;
  for(int i = 0; i < K; i++){
    diag_prod += std::log(*(chol_inv+K*i+i));
  }

  // standardized means? (based on arma::inplace_tri_mat_mult)
  double * z;
  z = new double [K];
  for(int i = 0; i < K; i++){
    *(z+i) = 0;
  }
  for(int i = 0; i < K; i++){
    *(z+i) += *(x+i) - *(mu+i);
  }
  for(int i = K - 1; i >= 0; i--){
    double temp = 0;
    for(int j = 0; j <= i; j++){
      temp += *(chol_inv+K*i+j) * *(z+j);
    }
    *(z+i) = temp;
  }

  // log lik
  double log_lik = 0;
  for(int i = 0; i < K; i++){
    log_lik += std::pow(*(z+i), 2);
  }
  log_lik = - 0.5 * log_lik + diag_prod - K * 0.9189385;

  // clean the memory
  delete[] sigma_chol;
  delete[] chol_inv;
  delete[] z;

  return log_lik;
}


// based on: https://www.geeksforgeeks.org/cholesky-decomposition-matrix-decomposition/
void chol(double const *matrix, const int n, double *lower)
{
	// Decomposing a matrix into Lower Triangular
	for (int i = 0; i < n; i++) {
		for (int j = 0; j <= i; j++) {
			double sum = 0;

			if (j == i) {
        // summation for diagnols
				for (int k = 0; k < j; k++)
					sum += std::pow(lower[n*j+k], 2);
				lower[n*j+j] = std::sqrt(*(matrix+n*j+j) -	sum);
			} else {
				// Evaluating L(i, j) using L(j, j)
				for (int k = 0; k < j; k++)
					sum += (lower[n*i+k] * lower[n*j+k]);
				lower[n*i+j] = (*(matrix+n*i+j) - sum) / lower[n*j+j];
			}
		}
	}
}

void cofactor(double const *matrix, double *temp, int p, int q, int n, int const K)
{
  int i = 0, j = 0;

  // Looping for each element of the matrix
  for (int row = 0; row < n; row++){
    for (int col = 0; col < n; col++){
      //  Copying into temporary matrix only those element
      //  which are not in given row and column
      if (row != p && col != q){
        temp[K*i+j++] = *(matrix+K*row+col);

        // Row is filled, so increase row index and
        // reset col index
        if (j == n - 1){
            j = 0;
            i++;
        }
      }
    }
  }
}

// Recursive function for finding determinant of matrix.
//   n is current dimension of A[][].
double determinant(double const *matrix, int n, int const K)
{
  double D = 0; // Initialize result

  //  Base case : if matrix contains single element
  if (n == 1)
    return *matrix;

  // To store cofactors
  double * temp;
  temp = new double [K*K];

  int sign = 1;  // To store sign multiplier

  // Iterate for each element of first row
  for (int f = 0; f < n; f++){
    // Getting Cofactor of A[0][f]
    cofactor(&matrix[0], temp, 0, f, n, K);
    D += sign * *(matrix+f) * determinant(temp, n - 1, K);

    // terms are to be added with alternate sign
    sign = -sign;
  }

  // clean the memory
  delete[] temp;

  return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(double const *matrix, double *adj, int const K)
{
    if (K == 1)
    {
        *adj = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    int sign = 1;

    double * temp;
    temp = new double [K*K];

    for (int i=0; i<K; i++)
    {
        for (int j=0; j<K; j++)
        {
            // Get cofactor of A[i][j]
            cofactor(&matrix[0], temp, i, j, K, K);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[K*j+i] = (sign)*(determinant(temp, K-1, K));
        }
    }

  // clean the memory
  delete[] temp;
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(double const *matrix, int const K, double *inverse)
{
  // Find determinant of A[][]
  double det = determinant(&matrix[0], K, K);
  if (det == 0){
    //cout << "Singular matrix, can't find its inverse" << endl;
    return false;
  }

  // Find adjoint
  double * adj;
  adj = new double [K*K];

  adjoint(&matrix[0], adj, K);

  // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
  for (int i=0; i<K; i++)
    for (int j=0; j<K; j++)
      inverse[K*i+j] = *(adj+K*i+j)/det;

  // clean the memory
  delete[] adj;

  return true;
}

