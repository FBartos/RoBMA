#ifndef MNORM_H_
#define MNORM_H_

// wrapper around the mvtnorm
double cpp_mnorm_cdf(double *lower, double *upper, int *infin, double *mu, double *sigma_stdev, double *sigma_corr, int K);
// adapted from https://gallery.rcpp.org/articles/dmvnorm_arma/
double cpp_mnorm_lpdf(double const *x, double const *mu, double const *sigma, const int K);

// adapted from https://www.geeksforgeeks.org/
void chol(double const *matrix, const int n, double *lower);
void cofactor(double const *matrix, double *temp, int p, int q, int n, int const K);
double determinant(double const *matrix, int n, int const K);
void adjoint(double const *matrix, double *adj, int const K);
bool inverse(double const *matrix, int const K, double *inverse);


#endif /* MNORM_H_ */
