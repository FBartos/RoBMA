#include "mnorm.h"
#include <iostream>
#include <mvtnormAPI.h>


using std::cout;
using std::endl;

// wrapper around the mvtnorm package
double pmnorm(double const *lower, double const *upper, double const *mu, double const *sigma_stdev, double const *sigma_corr, const int K)
{
  /*
  // create dynamically allocated arrays for the standardized locations
  double * lower_std;
  double * upper_std;
  double * delta;
  double * infin;

  lower_std = new double [K];
  upper_std = new double [K];
  delta     = new double [K];
  infin     = new double [K];

  // standardized boundary points
  for(int i = 0; i < K; i++){
    lower_std[i] = (lower[i] - mean[i]) / sigma_stdev[i];
    upper_std[i] = (upper[i] - mean[i]) / sigma_stdev[i];
  }

  // return 0 on the same upper and lower bounds
  for(int i = 0; i < K; i++){
    if(fabs(lower_std[i] - upper_std[i]) < EPSILON)
      return 0;
  }




  // mvtnorm settings
  double releps = 0;      // default in mvtnorm: 0
  int maxpts    = 25000;  // default in mvtnorm: 25000
  double abseps = 1e-3;   // default in mvtnorm: 0.001, absolute precision
  int rnd       = 1;      // Get/PutRNGstate
  int nu        = 0;      // degrees of freedom, 0 = normal

  // return values
  double error;
  double value;
  int    inform;

  mvtnorm_C_mvtdst(K,
                   nu,
                   lower_std,
                   upper_std,
                   infin,
                   sigma_corr,
                   delta,
                   maxpts, abseps, releps,
                   error, value, inform, rnd);

  return value;

*/

  double y = 0;
  return y;
}

double dmnorm(double const *x, double const *mu, double const *sigma, const int K)
{

  double * sigma_chol;
  double * chol_inv;
  sigma_chol = new double [K*K];
  chol_inv   = new double [K*K];

  // lower triangle cholesky decomposition
  chol(&sigma[0], K, sigma_chol);
  inverse(&sigma_chol[0], K, chol_inv);

  for(int i = 0; i < K*K; i++){
    cout << "sigma_chol[" << i << "]" << sigma_chol[i] << endl;
  }
    for(int i = 0; i < K*K; i++){
    cout << "chol_inv[" << i << "]" << chol_inv[i] << endl;
  }
  cout << "ggeeeg" << endl;
  /*
  Matrix3d m = Matrix3d::Random();
  m = (m + Matrix3d::Constant(1.2)) * 50;
  cout << "m =" << endl << m << endl;
  Vector3d v(1,2,3);

  cout << "m * v =" << endl << m * v << endl;


  arma::vec dmvnrm_arma_fast(arma::mat const &x,
                           arma::rowvec const &mean,
                           arma::mat const &sigma,
                           bool const logd = false) {
    using arma::uword;
    uword const n = x.n_rows,
             xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())),
                constants = -(double)xdim/2.0 * log2pi,
              other_terms = rootisum + constants;

    arma::rowvec z;
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);
        out(i) = other_terms - 0.5 * arma::dot(z, z);
    }

    if (logd)
      return out;
    return exp(out);
    */
   return 0;
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
					sum += pow(lower[n*j+k], 2);
				lower[n*j+j] = sqrt(*(matrix+n*j+j) -	sum);
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
        [temp+K*i+j++] = *(matrix+K*row+col);

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
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(double const *matrix, double const K, double *inverse)
{
  // Find determinant of A[][]
  double det = determinant(&matrix[0], K, K);
  if (det == 0){
    cout << "Singular matrix, can't find its inverse";
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

  return true;
}
