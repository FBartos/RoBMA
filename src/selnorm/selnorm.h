#ifndef SELNORM_H_
#define SELNORM_H_

enum SelKernelMode {
  SELKERNEL_NORMAL           = 0,
  SELKERNEL_STEP             = 1,
  SELKERNEL_PHACK_POWER      = 2,
  SELKERNEL_STEP_PHACK_POWER = 3
};

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
};

bool selnorm_is_descending_step_partition(const double *z_lower,
                                          const double *z_upper,
                                          int n_bins);

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

double cpp_selnorm_kernel_log_cdf(
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

double cpp_selnorm_kernel_cdf_with_log_norm(
  double q,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double log_normalizer,
  int omega_stride = 1,
  bool lower_tail = true,
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

void cpp_selnorm_kernel_summary(
  double q,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *cdf,
  double *moment_mean,
  double *moment_second,
  int omega_stride = 1,
  bool lower_tail = true,
  bool validate_omega = true
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

#endif
