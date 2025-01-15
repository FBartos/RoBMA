#ifndef WMNORM_H_
#define WMNORM_H_

double cpp_wnorm_1s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int J);
double cpp_wnorm_2s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int J);
double cpp_wmnorm_1s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int K, const int J);
double cpp_wmnorm_2s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int K, const int J);

#endif /* WMNORM_H_ */
