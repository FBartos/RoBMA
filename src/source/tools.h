#ifndef tools_H_
#define tools_H_

double log_weight_onesided(double const *x, double const *crit_x, double const *omega, const int J);
double log_weight_twosided(double const *x, double const *crit_x, double const *omega, const int J);

double log_std_constant_onesided(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega ,const int K, const int J);
double log_std_constant_twosided(double const *x, double const *mu, double const *sigma, double const *crit_x, double const *omega ,const int K, const int J);

int* create_index(int J);
void increase_index(int *index, int i, int max_i);

#endif /* tools_H_ */
