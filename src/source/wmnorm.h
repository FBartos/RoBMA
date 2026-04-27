#ifndef WMNORM_H_
#define WMNORM_H_

double cpp_wnorm_1s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int J);
double cpp_wnorm_2s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int J);
double cpp_wmnorm_1s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int K, const int J);
double cpp_wmnorm_2s_lpdf(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega, const int K, const int J);

double cpp_wnorm_mix_lpdf(double const x, double const mu, double const sigma, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, const int omega_stride = 1);
double cpp_wnorm_mix_cdf(double const q, double const mu, double const sigma, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, const bool lower_tail, const bool log_p, const int omega_stride = 1);
double cpp_wnorm_mix_log_norm(double const mean, double const sd, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, const int omega_stride = 1);
double cpp_wnorm_mix_lpdf_precomputed(double const x, double const mu, double const sigma, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, const double log_std_const, const int omega_stride = 1);
double cpp_wnorm_mix_cdf_precomputed(double const q, double const mu, double const sigma, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, const double log_denom_total, const bool lower_tail, const bool log_p, const int omega_stride = 1);
void cpp_wnorm_mix_moments(double const mu, double const sigma, double const *crit_y_all, double const *omega_all, double const *crit_y_mapping, const int crit_y_mapping_max, double *moment_mean, double *moment_second, const int omega_stride = 1);

#endif /* WMNORM_H_ */
