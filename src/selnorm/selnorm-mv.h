#ifndef ROBMA_SELNORM_MV_H
#define ROBMA_SELNORM_MV_H

double cpp_selnorm_mnorm_step_lpdf(
  const double *x,
  const double *mean,
  const double *covariance_lower,
  int dimension,
  const double *selection_se,
  const double *omega,
  int n_bins,
  const double *z_lower,
  const double *z_upper,
  const int *obs_bin,
  int effect_sign,
  bool telescope_probabilities,
  int kernel_mode,
  const double *qmc,
  int points,
  int scrambles,
  double *relative_mcse
);

#endif
