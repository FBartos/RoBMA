#ifndef ROBMA_SELNORM_MV_H
#define ROBMA_SELNORM_MV_H

// Rao-Blackwellized scalar marginals of a finite-vector selected Gaussian.
// The error is the largest absolute density MCSE, relative to the curve peak.
struct SelNormZProjection {
  const double *z;
  int size;
  double *density;
  double relative_error;
  bool probability;
  const double *factor_loading;
  const double *residual_sd;
  int rank;
  const double *nodes;
  const double *log_weights;
  int order;
};

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
  double *relative_mcse,
  SelNormZProjection *projection = nullptr
);

double cpp_selnorm_cluster_step_lpdf(
  const double *x,
  const double *mean,
  const double *residual_sd,
  const double *loading,
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
  const double *quadrature_nodes,
  const double *quadrature_log_weights,
  const double *quadrature_orders,
  int quadrature_rule_count,
  double relative_tolerance,
  double *relative_change,
  const long double *quadrature_weights = nullptr
);

double cpp_selnorm_factor_step_lpdf(
  const double *x,
  const double *mean,
  const double *residual_sd,
  const double *loading,
  int dimension,
  int rank,
  const double *selection_se,
  const double *omega,
  int n_bins,
  const double *z_lower,
  const double *z_upper,
  const int *obs_bin,
  int effect_sign,
  bool telescope_probabilities,
  int kernel_mode,
  const double *quadrature_nodes,
  const double *quadrature_log_weights,
  const double *quadrature_orders,
  int quadrature_rule_count,
  const double *qmc,
  int initial_points,
  int max_points,
  int scrambles,
  double relative_tolerance,
  double *relative_mcse,
  double *relative_change
);

#endif
