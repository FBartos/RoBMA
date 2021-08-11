#ifndef mnorm_H_
#define mnorm_H_

// wrapper around the mvtnorm
double pmnorm(double *lower, double *upper, int *infin, double *mu, double *sigma_stdev, double *sigma_corr, int K);
// adapted from https://gallery.rcpp.org/articles/dmvnorm_arma/
double dmnorm(double const *x, double const *mu, double const *sigma, const int K);

// adapted from https://www.geeksforgeeks.org/
void chol(double const *matrix, const int n, double *lower);
void cofactor(double const *matrix, double *temp, int p, int q, int n, int const K);
double determinant(double const *matrix, int n, int const K);
void adjoint(double const *matrix, double *adj, int const K);
bool inverse(double const *matrix, int const K, double *inverse);


#endif /* mnorm_H_ */
