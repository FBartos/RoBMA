#ifndef pmwnorm_H_
#define pmwnorm_H_

#include <vector>
using std::vector;

// wrapper around the mvtnorm
extern void mvtnorm_C_mvtdst(int *n, int *nu, double *lower, double *upper,
                      int *infin, double *corr, double *delta,
                      int *maxpts, double *abseps, double *releps,
                      double *error, double *value, int *inform, int *rnd);

double pmnorm(vector<double> const&lower, vector<double> const&upper,
                     vector<double> const&mu, vector<vector<double> > const&sigma);

#endif /* pmwnorm_H_ */
