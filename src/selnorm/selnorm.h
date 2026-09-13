#ifndef SELNORM_H_
#define SELNORM_H_

#include <cstddef>
#include "selnorm-fma.h"

enum SelKernelMode {
  SELKERNEL_NORMAL           = 0,
  SELKERNEL_STEP             = 1,
  SELKERNEL_PHACK_POWER      = 2,
  SELKERNEL_STEP_PHACK_POWER = 3
};

enum SelVectorRule {
  SELVECTOR_PRODUCT = 0,
  SELVECTOR_BEST_ONE_SIDED = 1,
  SELVECTOR_BEST_TWO_SIDED = 2
};

double cpp_selnorm_normal_interval_log_prob(
  double lower, double upper, double mean, double sd
);

double cpp_selnorm_normal_interval_quantile(
  double lower, double upper, double mean, double sd, double u
);

struct SelNormKernelData {
  int n_bins;
  int n_segments;
  int effect_sign;
  int q;
  const double *z_lower;
  const double *z_upper;
  const double *phack_z_source;
  const double *phack_z_dest;
  const double *segment_bounds;
  const int *segment_step_bin;
  const int *segment_phack_region;
  const double *segment_step_bin_real;
  const double *segment_phack_region_real;
  bool trusted_step_partition;
  bool telescope_probabilities;
  // Set only on a call-owned R post-fit integration copy.
  bool bounded_cdf = false;
};

double cpp_selnorm_bounded_cdf_absolute_error();

bool selnorm_is_descending_step_partition(const double *z_lower,
                                          const double *z_upper,
                                          int n_bins);

double cpp_selnorm_normal_lpdf(double x, double mean, double sd);

double cpp_selnorm_affine_normal_lpdf(double z, double sei, double mean,
                                      double sd);

double cpp_selnorm_affine_normal_lpdf_log_scale(
  double z, double sei, double mean, double sd, double log_scale
);

double cpp_selnorm_affine_normal_pdf_log_scale(
  double z, double sei, double mean, double sd, double log_scale
);

double cpp_selnorm_kernel_lpdf(
  double y,
  double mu_num,
  double sigma_num,
  double mu_norm,
  double sigma_norm,
  double sei,
  double weight,
  const double *omega,
  int obs_bin,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride = 1,
  bool validate_omega = true
);

double cpp_selnorm_kernel_log_norm(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride = 1,
  bool validate_omega = true
);

double cpp_selnorm_kernel_cdf(
  double q,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride = 1,
  bool lower_tail = true,
  bool validate_omega = true
);

bool cpp_selnorm_step_cdf_plan(
  double mean,
  double sd,
  double sei,
  const double *omega,
  const SelNormKernelData &data,
  double *lower_score,
  double *upper_score,
  double *log_weight,
  int *n_groups,
  int omega_stride = 1,
  bool validate_omega = true
);

double cpp_selnorm_step_log_norm(
  double mean,
  double sd,
  double sei,
  const double *omega,
  const SelNormKernelData &data,
  int omega_stride = 1,
  bool validate_omega = true
);

bool cpp_selnorm_step_cdf_telescope_plan(
  double mean,
  double sd,
  double sei,
  const double *omega,
  const SelNormKernelData &data,
  double *boundary_tail,
  double *omega_diff,
  double *omega_last,
  double *normalizer,
  int omega_stride = 1,
  bool validate_omega = true
);

// Multiply one row's normalizers over a quadrature grid. The caller validates
// omega once per likelihood state; false retains the scalar/log-scale fallback.
bool cpp_selnorm_step_normalizer_product(
  const double *mean, std::size_t count, double sd, double sei,
  const double *omega, const SelNormKernelData &data, long double *product
);

double cpp_selnorm_step_cdf_from_telescope_plan(
  double q,
  double mean,
  double sd,
  double sei,
  const SelNormKernelData &data,
  const double *boundary_tail,
  const double *omega_diff,
  double omega_last,
  double normalizer,
  bool lower_tail = true
);

double cpp_selnorm_step_cdf_from_plan(
  double q,
  double mean,
  double sd,
  double sei,
  const SelNormKernelData &data,
  const double *lower_score,
  const double *upper_score,
  const double *log_weight,
  int n_groups,
  bool lower_tail = true
);

void cpp_selnorm_kernel_log_tail_pair(
  double q,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *log_lower,
  double *log_upper,
  int omega_stride = 1,
  bool validate_omega = true
);

void cpp_selnorm_kernel_moments(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *moment_mean,
  double *moment_second,
  int omega_stride = 1,
  bool validate_omega = true,
  double *moment_variance = nullptr
);

double cpp_selnorm_kernel_threshold(
  double z_threshold,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *inverse_weight,
  int omega_stride = 1,
  bool validate_omega = true
);

double cpp_selnorm_kernel_rng(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double u_bin,
  double u_interval,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride = 1,
  bool validate_omega = true
);

double cpp_selnorm_kernel_rng_workspace(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double u_bin,
  double u_interval,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *mass,
  double *lower,
  double *upper,
  int omega_stride = 1,
  bool validate_omega = true
);

bool cpp_selnorm_step_rng_prepare_partition(
  double mean, double sd, double sei, const double *omega,
  const SelNormKernelData &data, double *mass, double *lower, double *upper,
  int *n_groups, double *normalizer, double *log_normalizer
);

double cpp_selnorm_step_rng_from_partition(
  double mean, double sd, int effect_sign, double u_bin, double u_interval,
  const double *mass, const double *lower, const double *upper,
  int n_groups, double normalizer
);

double cpp_selnorm_step_log_norm_rng_workspace(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double u_bin,
  double u_interval,
  const SelNormKernelData &data,
  double *mass,
  double *lower,
  double *upper,
  double *log_normalizer,
  int omega_stride = 1,
  bool validate_omega = true
);

#endif
