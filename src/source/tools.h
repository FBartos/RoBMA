#ifndef tools_H_
#define tools_H_

#include <vector>

double log_weight_onesided(double const *x, double const *crit_x, double const *omega, const int J);
double log_weight_twosided(double const *x, double const *crit_x, double const *omega, const int J);

double log_std_constant_onesided(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int J);
double log_std_constant_twosided(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int J);
double log_std_m_constant_onesided(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int K, const int J);
double log_std_m_constant_twosided(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int K, const int J);

void increase_index(int *index, int i, int max_i);

double * extract_x_v(const double *x_v, int indx_start, int K);
double * extract_mu_v(const double *mu_v,  int indx_start, int K);
double * extract_sigma_v(const double *se2_v, const double *tau2, double cov, int indx_start, int K);
double * extract_crit_x_v(const double *crit_x_v, int indx_start, int K, const int J);

double ddirichlet(const std::vector<double>& x, const std::vector<double>& alpha);

#endif /* tools_H_ */
