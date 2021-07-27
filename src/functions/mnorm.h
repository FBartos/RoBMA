#ifndef mnorm_H_
#define mnorm_H_

// wrapper around the mvtnorm
double pmnorm(double const *lower, double const *upper, double const *mu, double const *sigma_stdev, double const *sigma_corr, const int K);
// adapted from https://gallery.rcpp.org/articles/dmvnorm_arma/
double dmnorm(double const *x, double const *mu, double const *sigma_stdev, double const *sigma_corr, const int K);

#endif /* mnorm_H_ */
