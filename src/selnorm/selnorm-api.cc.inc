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
  int omega_stride,
  bool validate_omega)
{
  if (!(sigma_num > 0) || !(sigma_norm > 0) || !(sei > 0) || !(weight > 0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (!valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const int mode = kernel_mode;
  const int q = phack_kind > 0 ? phack_kind : data.q;
  const bool use_step = mode == SELKERNEL_STEP ||
    mode == SELKERNEL_STEP_PHACK_POWER;
  const bool use_phack = (mode == SELKERNEL_PHACK_POWER ||
    mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;

  const double signed_y       = data.effect_sign * y;
  const double signed_mu_num  = data.effect_sign * mu_num;
  double log_lik = dnorm(signed_y, signed_mu_num, sigma_num, true);

  if (mode == SELKERNEL_NORMAL) {
    return weight * log_lik;
  }

  if (use_phack && (!(alpha >= 0) || !(alpha < 1) || !valid_phack_kind(q))) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double signed_mu_norm = data.effect_sign * mu_norm;
  const double mean_z         = signed_mu_norm / sei;
  const double sd_z           = sigma_norm / sei;

  double log_beta = -INFINITY;
  if (use_phack) {
    log_beta = phack_log_beta(alpha, q, data, mean_z, sd_z);
    const double z = signed_y / sei;
    const double log_r = phack_log_kernel_weight(
      z, alpha, log_beta, q, data
    );
    if (log_r == -INFINITY) {
      return -INFINITY;
    }
    if (!std::isfinite(log_r)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    log_lik += log_r;
  }

  if (use_step) {
    if (obs_bin < 1 || obs_bin > data.n_bins) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    if (validate_omega) {
      if (!validate_omega_values(omega, data, omega_stride)) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    }

    const double observed_weight = omega_at(omega, obs_bin, omega_stride);
    if (!(observed_weight > 0)) {
      return -INFINITY;
    }
    if (observed_weight != 1) {
      log_lik += std::log(observed_weight);
    }
  }

  double log_normalizer = 0;
  if (mode == SELKERNEL_STEP ||
      (mode == SELKERNEL_STEP_PHACK_POWER && !use_phack)) {
    log_normalizer = step_log_normalizer(
      omega, data, mean_z, sd_z, omega_stride
    );
  } else if (mode == SELKERNEL_STEP_PHACK_POWER) {
    log_normalizer = combined_log_normalizer(
      omega, use_phack ? alpha : 0, log_beta, use_phack ? q : 0,
      data, mean_z, sd_z, omega_stride
    );
  }

  if (!std::isfinite(log_normalizer)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  log_lik -= log_normalizer;
  return weight == 1 ? log_lik : weight * log_lik;
}

double cpp_selnorm_kernel_log_norm(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride,
  bool validate_omega)
{
  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (kernel_mode == SELKERNEL_NORMAL ||
      kernel_mode == SELKERNEL_PHACK_POWER) {
    return 0;
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const int q = phack_kind > 0 ? phack_kind : data.q;
  const bool use_phack = kernel_mode == SELKERNEL_STEP_PHACK_POWER &&
    phack_kind > 0 && alpha > 0;

  if (use_phack && (!(alpha >= 0) || !(alpha < 1) || !valid_phack_kind(q))) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double mean_z = data.effect_sign * mean / sei;
  const double sd_z   = sd / sei;
  double log_normalizer = 0;

  if (kernel_mode == SELKERNEL_STEP ||
      (kernel_mode == SELKERNEL_STEP_PHACK_POWER && !use_phack)) {
    log_normalizer = step_log_normalizer(
      omega, data, mean_z, sd_z, omega_stride
    );
  } else if (kernel_mode == SELKERNEL_STEP_PHACK_POWER) {
    const double log_beta = use_phack ?
      phack_log_beta(alpha, q, data, mean_z, sd_z) : -INFINITY;
    log_normalizer = combined_log_normalizer(
      omega, use_phack ? alpha : 0, log_beta, use_phack ? q : 0,
      data, mean_z, sd_z, omega_stride
    );
  }

  if (!std::isfinite(log_normalizer)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return log_normalizer;
}

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
  int omega_stride,
  bool lower_tail,
  bool validate_omega)
{
  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    return pnorm(q, mean, sd, lower_tail, false);
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    return pnorm(q, mean, sd, lower_tail, false);
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double signed_q = data.effect_sign * q;
  const double q_z      = signed_q / sei;
  const double mean_z   = data.effect_sign * mean / sei;
  const double sd_z     = sd / sei;
  const bool signed_lower_tail = data.effect_sign == 1 ?
    lower_tail : !lower_tail;

  if (data.trusted_step_partition && data.telescope_probabilities) {
    bool ok = false;
    const double out = step_cdf_trusted_telescoping(
      omega, data, q_z, mean_z, sd_z, omega_stride, signed_lower_tail, &ok
    );
    if (ok) {
      return out;
    }
  }

  return step_cdf_logspace(
    omega, data, q_z, mean_z, sd_z, omega_stride, signed_lower_tail,
    data.trusted_step_partition
  );
}

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
  int omega_stride,
  bool lower_tail,
  bool validate_omega)
{
  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    return pnorm(q, mean, sd, lower_tail, true);
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    return pnorm(q, mean, sd, lower_tail, true);
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double signed_q = data.effect_sign * q;
  const double q_z      = signed_q / sei;
  const double mean_z   = data.effect_sign * mean / sei;
  const double sd_z     = sd / sei;
  const bool signed_lower_tail = data.effect_sign == 1 ?
    lower_tail : !lower_tail;

  return step_log_cdf_logspace(
    omega, data, q_z, mean_z, sd_z, omega_stride, signed_lower_tail,
    data.trusted_step_partition
  );
}

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
  int omega_stride,
  bool lower_tail,
  bool validate_omega)
{
  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    return pnorm(q, mean, sd, lower_tail, false);
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    return pnorm(q, mean, sd, lower_tail, false);
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double signed_q = data.effect_sign * q;
  const double q_z      = signed_q / sei;
  const double mean_z   = data.effect_sign * mean / sei;
  const double sd_z     = sd / sei;
  const bool signed_lower_tail = data.effect_sign == 1 ?
    lower_tail : !lower_tail;

  if (data.trusted_step_partition && data.telescope_probabilities) {
    bool ok = false;
    const double out = step_cdf_trusted_telescoping_with_log_norm(
      omega, data, q_z, mean_z, sd_z, omega_stride, signed_lower_tail,
      log_normalizer, &ok
    );
    if (ok) {
      return out;
    }
  }

  return step_cdf_logspace_with_log_norm(
    omega, data, q_z, mean_z, sd_z, omega_stride, signed_lower_tail,
    data.trusted_step_partition, log_normalizer
  );
}

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
  int omega_stride,
  bool validate_omega,
  double *moment_variance)
{
  *moment_mean   = std::numeric_limits<double>::quiet_NaN();
  *moment_second = std::numeric_limits<double>::quiet_NaN();
  if (moment_variance != NULL) {
    *moment_variance = std::numeric_limits<double>::quiet_NaN();
  }

  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return;
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    *moment_mean   = mean;
    *moment_second = sd * sd + mean * mean;
    if (moment_variance != NULL) {
      *moment_variance = sd * sd;
    }
    return;
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return;
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    *moment_mean   = mean;
    *moment_second = sd * sd + mean * mean;
    if (moment_variance != NULL) {
      *moment_variance = sd * sd;
    }
    return;
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return;
  }

  const double signed_mean = data.effect_sign * mean;
  double denom;
  double signed_first;
  double signed_second;
  double centered_first;
  double centered_second;

  if (!(moment_variance == NULL && data.trusted_step_partition &&
        data.telescope_probabilities &&
        step_moments_trusted_telescoping(
          omega, data, signed_mean, sd, sei, omega_stride,
          &denom, &signed_first, &signed_second
        ))) {
    step_moments_partition(
      omega, data, signed_mean, sd, sei, omega_stride,
      data.trusted_step_partition, &denom, &signed_first, &signed_second,
      moment_variance == NULL ? NULL : &centered_first,
      moment_variance == NULL ? NULL : &centered_second
    );
  }

  if (!(denom > 0)) {
    return;
  }

  if (moment_variance != NULL) {
    const double centered_mean = centered_first / denom;
    *moment_mean     = data.effect_sign * (signed_mean + centered_mean);
    *moment_variance = centered_second / denom - centered_mean * centered_mean;
  } else {
    *moment_mean = data.effect_sign * signed_first / denom;
  }
  *moment_second = signed_second / denom;
}

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
  int omega_stride,
  bool lower_tail,
  bool validate_omega)
{
  *cdf           = std::numeric_limits<double>::quiet_NaN();
  *moment_mean   = std::numeric_limits<double>::quiet_NaN();
  *moment_second = std::numeric_limits<double>::quiet_NaN();

  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return;
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    *cdf           = pnorm(q, mean, sd, lower_tail, false);
    *moment_mean   = mean;
    *moment_second = sd * sd + mean * mean;
    return;
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return;
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    *cdf           = pnorm(q, mean, sd, lower_tail, false);
    *moment_mean   = mean;
    *moment_second = sd * sd + mean * mean;
    return;
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return;
  }

  const double signed_q          = data.effect_sign * q;
  const double q_z               = signed_q / sei;
  const double signed_mean       = data.effect_sign * mean;
  const double mean_z            = signed_mean / sei;
  const double sd_z              = sd / sei;
  const bool signed_lower_tail   = data.effect_sign == 1 ?
    lower_tail : !lower_tail;

  if (data.trusted_step_partition && data.telescope_probabilities &&
      step_summary_trusted_telescoping(
        omega, data, q_z, mean_z, sd_z, signed_mean, sd, sei,
        omega_stride, signed_lower_tail, cdf, moment_mean, moment_second
      )) {
    return;
  }

  cpp_selnorm_kernel_moments(
    mean, sd, sei, omega, alpha, phack_kind, kernel_mode, data,
    moment_mean, moment_second, omega_stride, false
  );
  *cdf = cpp_selnorm_kernel_cdf(
    q, mean, sd, sei, omega, alpha, phack_kind, kernel_mode, data,
    omega_stride, lower_tail, false
  );
}

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
  int omega_stride,
  bool validate_omega)
{
  *inverse_weight = std::numeric_limits<double>::quiet_NaN();

  if (!(z_threshold >= 0) || !(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double q_upper = z_threshold * sei;
  const double q_lower = -q_upper;

  if (kernel_mode == SELKERNEL_NORMAL) {
    *inverse_weight = 1;
    return pnorm(q_upper, mean, sd, false, false) +
      pnorm(q_lower, mean, sd, true, false);
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    *inverse_weight = 1;
    return pnorm(q_upper, mean, sd, false, false) +
      pnorm(q_lower, mean, sd, true, false);
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double mean_z = data.effect_sign * mean / sei;
  const double sd_z   = sd / sei;

  if (data.trusted_step_partition && data.telescope_probabilities) {
    double upper_prob;
    double lower_prob;
    const double q_upper_z = data.effect_sign * z_threshold;
    const double q_lower_z = -q_upper_z;
    const bool upper_signed_lower_tail = data.effect_sign == 1 ? false : true;
    const bool lower_signed_lower_tail = data.effect_sign == 1 ? true : false;

    if (step_cdf_pair_trusted_telescoping(
          omega, data,
          q_upper_z, upper_signed_lower_tail,
          q_lower_z, lower_signed_lower_tail,
          mean_z, sd_z, omega_stride,
          &upper_prob, &lower_prob, inverse_weight
        )) {
      return upper_prob + lower_prob;
    }
  }

  double upper_prob;
  double lower_prob;
  if (step_cdf_pair_logspace(
        omega, data,
        data.effect_sign * z_threshold,
        data.effect_sign == 1 ? false : true,
        -data.effect_sign * z_threshold,
        data.effect_sign == 1 ? true : false,
        mean_z, sd_z, omega_stride, data.trusted_step_partition,
        &upper_prob, &lower_prob, inverse_weight
      )) {
    return upper_prob + lower_prob;
  }

  return std::numeric_limits<double>::quiet_NaN();
}

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
  int omega_stride,
  bool validate_omega)
{
  if (!(sd > 0) || !(sei > 0) ||
      !(u_bin >= 0) || !(u_bin <= 1) ||
      !(u_interval >= 0) || !(u_interval <= 1) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    return qnorm(u_interval, mean, sd, true, false);
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    return qnorm(u_interval, mean, sd, true, false);
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double mean_z = data.effect_sign * mean / sei;
  const double sd_z   = sd / sei;

  if (data.trusted_step_partition && data.telescope_probabilities) {
    bool ok = false;
    const double z = step_rng_trusted_partition(
      omega, data, mean_z, sd_z, u_bin, u_interval, omega_stride, &ok
    );
    if (ok) {
      return data.effect_sign * z * sei;
    }
  }

  const double z = step_rng_logspace(
    omega, data, mean_z, sd_z, u_bin, u_interval, omega_stride,
    data.trusted_step_partition
  );

  return data.effect_sign * z * sei;
}
