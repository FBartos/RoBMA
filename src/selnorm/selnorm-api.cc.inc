#include "selnorm-cdf-table-data.h"

double cpp_selnorm_bounded_cdf_absolute_error()
{
  return selnorm_cdf_table::absolute_error;
}

double cpp_selnorm_normal_lpdf(double x, double mean, double sd)
{
  if (!(sd > 0) || !std::isfinite(sd) || !std::isfinite(mean) ||
      std::isnan(x)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double score = stable_normal_score(x, mean, sd);
  return dnorm(score, 0, 1, true) - std::log(sd);
}

double cpp_selnorm_affine_normal_lpdf(double z, double sei, double mean,
                                      double sd)
{
  return cpp_selnorm_affine_normal_lpdf_log_scale(
    z, sei, mean, sd, std::log(sei) - std::log(sd)
  );
}

double cpp_selnorm_affine_normal_lpdf_log_scale(
  double z, double sei, double mean, double sd, double log_scale)
{
  if (!(sei > 0) || !std::isfinite(sei) ||
      !(sd > 0) || !std::isfinite(sd) || !std::isfinite(mean) ||
      std::isnan(z)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double score = stable_affine_score(z, sei, mean, sd);
  return dnorm(score, 0, 1, true) + log_scale;
}

double cpp_selnorm_affine_normal_pdf_log_scale(
  double z, double sei, double mean, double sd, double log_scale)
{
  if (!(sei > 0) || !std::isfinite(sei) ||
      !(sd > 0) || !std::isfinite(sd) || !std::isfinite(mean) ||
      std::isnan(z)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double score = stable_affine_score(z, sei, mean, sd);
  const double density = Rf_dnorm4(score, 0, 1, false);
  const double scale = sei / sd;
  if (std::isnormal(density) && std::isnormal(scale)) {
    return density * scale;
  }
  // Retain the log-scale evaluation when either factor would lose tail mass
  // to underflow, or the Jacobian would overflow before multiplication.
  return std::exp(dnorm(score, 0, 1, true) + log_scale);
}

double cpp_selnorm_normal_interval_log_prob(
  double lower, double upper, double mean, double sd)
{
  return interval_log_prob(lower, upper, mean, sd);
}

double cpp_selnorm_normal_interval_quantile(
  double lower, double upper, double mean, double sd, double u)
{
  return truncated_normal_quantile(lower, upper, mean, sd, u);
}

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
  if (!(sigma_num > 0) || !std::isfinite(sigma_num) ||
      !(sigma_norm > 0) || !std::isfinite(sigma_norm) ||
      !(sei > 0) || !std::isfinite(sei) ||
      !std::isfinite(mu_num) || !std::isfinite(mu_norm) ||
      !(weight > 0) || !std::isfinite(weight)) {
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
  double log_lik = cpp_selnorm_normal_lpdf(
    signed_y, signed_mu_num, sigma_num
  );

  if (mode == SELKERNEL_NORMAL) {
    return weight * log_lik;
  }

  if (use_phack && (!(alpha >= 0) || !(alpha < 1) || !valid_phack_kind(q))) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double signed_mu_norm = data.effect_sign * mu_norm;
  double mean_z = 0;
  double sd_z   = 1;
  double log_beta = -INFINITY;
  if (use_phack) {
    mean_z  = signed_mu_norm / sei;
    sd_z    = sigma_norm / sei;
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
    log_normalizer = cpp_selnorm_step_log_norm(
      mu_norm, sigma_norm, sei, omega, data, omega_stride, false
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
  if (!(sd > 0) || !std::isfinite(sd) ||
      !(sei > 0) || !std::isfinite(sei) || !std::isfinite(mean) ||
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

  double mean_z = 0;
  double sd_z   = 1;
  if (use_phack) {
    mean_z = data.effect_sign * mean / sei;
    sd_z   = sd / sei;
  }
  double log_normalizer = 0;

  if (kernel_mode == SELKERNEL_STEP ||
      (kernel_mode == SELKERNEL_STEP_PHACK_POWER && !use_phack)) {
    log_normalizer = cpp_selnorm_step_log_norm(
      mean, sd, sei, omega, data, omega_stride, false
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


double cpp_selnorm_step_log_norm(
  double mean,
  double sd,
  double sei,
  const double *omega,
  const SelNormKernelData &data,
  int omega_stride,
  bool validate_omega)
{
  if (!(sd > 0) || !std::isfinite(sd) ||
      !(sei > 0) || !std::isfinite(sei) || !std::isfinite(mean) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      (validate_omega && !validate_omega_values(omega, data, omega_stride))) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (data.trusted_step_partition && data.telescope_probabilities) {
    double omega_last = 0;
    double normalizer = 0;
    if (step_normalizer_trusted_telescoping(
          omega, data, data.effect_sign * mean, sd, sei, omega_stride,
          NULL, NULL, &omega_last, &normalizer
        )) {
      return std::log(normalizer);
    }

    int n_groups = 0;
    bool ok = false;
    const double fallback = step_normalizer_trusted_partition(
      omega, data, data.effect_sign * mean, sd, sei, omega_stride,
      NULL, NULL, NULL, &n_groups, &ok
    );
    if (ok) return std::log(fallback);
  }

  return step_log_normalizer_score_space(
    omega, data, data.effect_sign * mean, sd, sei, omega_stride
  );
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
  double log_lower;
  double log_upper;
  cpp_selnorm_kernel_log_tail_pair(
    q, mean, sd, sei, omega, alpha, phack_kind, kernel_mode, data,
    &log_lower, &log_upper, omega_stride, validate_omega
  );
  const double log_probability = lower_tail ? log_lower : log_upper;
  if (std::isnan(log_probability)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return probability_from_log_ratio(log_probability, 0);
}

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
  int omega_stride,
  bool validate_omega)
{
  *n_groups = 0;
  if (!(sd > 0) || !std::isfinite(sd) ||
      !(sei > 0) || !std::isfinite(sei) ||
      !std::isfinite(mean) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      (validate_omega &&
       !validate_omega_values(omega, data, omega_stride))) {
    return false;
  }

  const double signed_mean = data.effect_sign * mean;
  for (int b = 0; b < data.n_bins; ) {
    const double w = omega_at(omega, b + 1, omega_stride);
    int last = b;

    if (data.trusted_step_partition) {
      while (last + 1 < data.n_bins &&
             omega_at(omega, last + 2, omega_stride) == w) {
        ++last;
      }
    }

    if (w > 0) {
      const double lower_z = data.trusted_step_partition ?
        data.z_lower[last] : data.z_lower[b];
      const double upper_z = data.z_upper[b];
      const double lower = stable_affine_score(
        lower_z, sei, signed_mean, sd
      );
      const double upper = stable_affine_score(
        upper_z, sei, signed_mean, sd
      );
      if (std::isnan(lower) || std::isnan(upper) || !(lower < upper)) {
        return false;
      }

      const int group = *n_groups;
      lower_score[group] = lower;
      upper_score[group] = upper;
      log_weight[group]  = std::log(w);
      ++(*n_groups);
    } else if (w < 0 || std::isnan(w)) {
      return false;
    }

    b = data.trusted_step_partition ? last + 1 : b + 1;
  }

  return *n_groups > 0;
}


// Telescoped step normalizers for a whole quadrature batch. For one row the
// cutoff scores are affine in the quadrature location, so the row's scale and
// scaled cutoffs hoist out of the loop and every tail is evaluated through the
// batched kernel. Returns false without a usable result whenever the row's
// geometry or any normalizer falls outside the direct form, and the caller's
// unchanged log-scale evaluator takes over.
bool step_normalizer_batch_telescoping(
  const double *mean, std::size_t count, double sd, double sei,
  const double *omega, const SelNormKernelData &data, long double *product)
{
  if (!data.trusted_step_partition || !data.telescope_probabilities ||
      data.n_bins < 1) {
    return false;
  }
  const double last = omega[data.n_bins - 1];
  if (!std::isfinite(last)) return false;
  const double inverse_sd = data.effect_sign / sd;
  if (!std::isfinite(inverse_sd)) return false;

  // Scaled cutoffs and weight differences are fixed over the whole batch.
  int boundaries = 0;
  double offset[SELNORM_TAIL_MAX_BINS];
  double difference[SELNORM_TAIL_MAX_BINS];
  if (data.n_bins - 1 > SELNORM_TAIL_MAX_BINS) return false;
  for (int bin = 0; bin < data.n_bins - 1; ++bin) {
    const double step = omega[bin] - omega[bin + 1];
    if (!std::isfinite(step)) return false;
    if (step == 0.0) continue;
    const double scaled = data.z_lower[bin] * sei / sd;
    if (!std::isfinite(scaled)) return false;
    offset[boundaries] = scaled;
    difference[boundaries] = step;
    ++boundaries;
  }

  const std::size_t block = 256;
  double tail[block];
  double normalizer[block];
  double magnitude[block];
  for (std::size_t start = 0; start < count; start += block) {
    const std::size_t size = std::min(block, count - start);
    for (std::size_t index = 0; index < size; ++index) {
      if (!std::isfinite(mean[start + index])) return false;
      normalizer[index] = last;
      magnitude[index] = std::abs(last);
    }
    for (int bin = 0; bin < boundaries; ++bin) {
      selnorm_tail::upper_tail_affine(mean + start, size, offset[bin],
                                      inverse_sd, tail);
      for (std::size_t index = 0; index < size; ++index) {
        const double contribution = difference[bin] * tail[index];
        normalizer[index] += contribution;
        magnitude[index] += std::abs(contribution);
      }
    }
    for (std::size_t index = 0; index < size; ++index) {
      const double value = normalizer[index];
      if (!std::isfinite(value) ||
          !(value > std::numeric_limits<double>::min()) ||
          !(value > 1e-8 * magnitude[index])) {
        return false;
      }
      long double &target = product[start + index];
      target *= static_cast<long double>(value);
      if (!(target > 0.0L) || !std::isfinite(target)) return false;
    }
  }
  return true;
}


bool cpp_selnorm_step_normalizer_product(
  const double *mean, std::size_t count, double sd, double sei,
  const double *omega, const SelNormKernelData &data, long double *product)
{
  if (!(sd > 0) || !std::isfinite(sd) ||
      !(sei > 0) || !std::isfinite(sei) || !data.trusted_step_partition ||
      !data.telescope_probabilities || data.n_bins < 1 ||
      !(data.effect_sign == 1 || data.effect_sign == -1)) {
    return false;
  }
  if (data.bounded_cdf && data.n_bins == 2 &&
      std::isfinite(data.z_lower[0]) && std::isfinite(omega[1])) {
    // These weights and their difference are fixed over every quadrature node.
    const bool bounded = data.bounded_cdf &&
      omega[0] >= std::numeric_limits<double>::min() &&
      omega[1] >= std::numeric_limits<double>::min() &&
      std::isfinite(omega[0]) && std::isfinite(omega[1]);
    const double last = bounded ? std::min(omega[0], omega[1]) : omega[1];
    const double difference = bounded ? std::max(omega[0], omega[1]) - last : omega[0] - last;
    const double direction = bounded && omega[0] < omega[1] ? -1.0 : 1.0;
    const double absolute_last = std::abs(last);
    for (std::size_t i = 0; i < count; ++i) {
      if (!std::isfinite(mean[i])) return false;
      const double score = stable_affine_score(
        data.z_lower[0], sei, data.effect_sign * mean[i], sd
      );
      double probability;
      if (!bounded || !selnorm_cdf_table::evaluate(-direction * score, &probability)) {
        probability = 0.5 * std::erfc(direction * score * M_SQRT1_2);
      }
      const double contribution = difference * probability;
      const double normalizer = last + contribution;
      if (!std::isfinite(normalizer) ||
          !(normalizer > std::numeric_limits<double>::min()) ||
          !(normalizer > 1e-8 * (absolute_last + std::abs(contribution)))) {
        return false;
      }
      product[i] *= static_cast<long double>(normalizer);
      if (!(product[i] > 0.0L) || !std::isfinite(product[i])) return false;
    }
    return true;
  }
  if (step_normalizer_batch_telescoping(mean, count, sd, sei, omega, data,
                                       product)) {
    return true;
  }

  for (std::size_t i = 0; i < count; ++i) {
    if (!std::isfinite(mean[i])) return false;
    double omega_last = 0.0;
    double normalizer = 0.0;
    if (!step_normalizer_trusted_telescoping(
          omega, data, data.effect_sign * mean[i], sd, sei, 1,
          nullptr, nullptr, &omega_last, &normalizer)) {
      return false;
    }
    product[i] *= static_cast<long double>(normalizer);
    if (!(product[i] > 0.0L) || !std::isfinite(product[i])) return false;
  }
  return true;
}

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
  int omega_stride,
  bool validate_omega)
{
  if (!(sd > 0) || !std::isfinite(sd) ||
      !(sei > 0) || !std::isfinite(sei) ||
      !std::isfinite(mean) || !data.trusted_step_partition ||
      !data.telescope_probabilities || data.n_bins < 1 ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      (validate_omega &&
       !validate_omega_values(omega, data, omega_stride))) {
    return false;
  }

  return step_normalizer_trusted_telescoping(
    omega, data, data.effect_sign * mean, sd, sei, omega_stride,
    boundary_tail, omega_diff, omega_last, normalizer
  );
}

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
  bool lower_tail)
{
  if (!(sd > 0) || !std::isfinite(sd) ||
      !(sei > 0) || !std::isfinite(sei) ||
      !std::isfinite(mean) || std::isnan(q) ||
      !data.trusted_step_partition || !data.telescope_probabilities ||
      data.n_bins < 1 || !(normalizer > 0) ||
      !std::isfinite(normalizer) ||
      !(data.effect_sign == 1 || data.effect_sign == -1)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double signed_q    = data.effect_sign * q;
  const double signed_mean = data.effect_sign * mean;
  const double q_score     = stable_normal_score(
    signed_q, signed_mean, sd
  );
  const double q_tail = tail_at(q_score, 0, 1);
  const double q_cdf  = cdf_at(q_score, 0, 1);
  if (!std::isfinite(q_tail) || !std::isfinite(q_cdf)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const bool signed_lower_tail = data.effect_sign == 1 ?
    lower_tail : !lower_tail;
  double requested = omega_last *
    (signed_lower_tail ? q_cdf : q_tail);

  for (int b = 0; b < data.n_bins - 1; ++b) {
    const double boundary = data.z_lower[b] * sei;
    double contribution;
    if (signed_lower_tail) {
      contribution = 0;
      if (signed_q > boundary) {
        const double interval = boundary_tail[b] - q_tail;
        if (!(interval > 0)) {
          return std::numeric_limits<double>::quiet_NaN();
        }
        contribution = omega_diff[b] * interval;
      }
    } else {
      contribution = omega_diff[b] *
        (signed_q >= boundary ? q_tail : boundary_tail[b]);
    }
    requested += contribution;
  }

  if (!std::isfinite(requested) || requested < 0 || requested > normalizer) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double probability = requested / normalizer;
  return std::isfinite(probability) && probability >= 0 && probability <= 1 ?
    probability : std::numeric_limits<double>::quiet_NaN();
}

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
  bool lower_tail)
{
  if (!(sd > 0) || !std::isfinite(sd) ||
      !(sei > 0) || !std::isfinite(sei) ||
      !std::isfinite(mean) || std::isnan(q) || n_groups < 1 ||
      !(data.effect_sign == 1 || data.effect_sign == -1)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double signed_q    = data.effect_sign * q;
  const double signed_mean = data.effect_sign * mean;
  const double q_score     = stable_normal_score(
    signed_q, signed_mean, sd
  );
  double signed_log_lower;
  double signed_log_upper;
  if (!step_log_tail_pair_from_plan(
        q_score, lower_score, upper_score, log_weight, n_groups,
        &signed_log_lower, &signed_log_upper
      )) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double log_lower = data.effect_sign == 1 ?
    signed_log_lower : signed_log_upper;
  const double log_upper = data.effect_sign == 1 ?
    signed_log_upper : signed_log_lower;
  return probability_from_log_ratio(
    lower_tail ? log_lower : log_upper, 0
  );
}

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
  int omega_stride,
  bool validate_omega)
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  *log_lower = nan;
  *log_upper = nan;
  if (!(sd > 0) || !std::isfinite(sd) ||
      !(sei > 0) || !std::isfinite(sei) ||
      !std::isfinite(mean) || std::isnan(q) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return;
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    const double score = stable_normal_score(q, mean, sd);
    *log_lower = pnorm(score, 0, 1, true, true);
    *log_upper = pnorm(score, 0, 1, false, true);
    return;
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return;
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    const double score = stable_normal_score(q, mean, sd);
    *log_lower = pnorm(score, 0, 1, true, true);
    *log_upper = pnorm(score, 0, 1, false, true);
    return;
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return;
  }

  const double signed_q    = data.effect_sign * q;
  const double signed_mean = data.effect_sign * mean;
  const double q_score     = stable_normal_score(signed_q, signed_mean, sd);
  double signed_log_lower;
  double signed_log_upper;
  if (!step_log_tail_pair_logspace(
    omega, data, q_score, signed_mean, sd, sei, omega_stride,
    data.trusted_step_partition, &signed_log_lower, &signed_log_upper
  )) {
    return;
  }

  if (data.effect_sign == 1) {
    *log_lower = signed_log_lower;
    *log_upper = signed_log_upper;
  } else {
    *log_lower = signed_log_upper;
    *log_upper = signed_log_lower;
  }
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

  if (!std::isfinite(mean) || !std::isfinite(sd) ||
      !std::isfinite(sei) || !(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return;
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    const double second = sd * sd + mean * mean;
    if (moment_variance == NULL && !std::isfinite(second)) {
      return;
    }
    *moment_mean   = mean;
    if (std::isfinite(second)) {
      *moment_second = second;
    }
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
    const double second = sd * sd + mean * mean;
    if (moment_variance == NULL && !std::isfinite(second)) {
      return;
    }
    *moment_mean   = mean;
    if (std::isfinite(second)) {
      *moment_second = second;
    }
    if (moment_variance != NULL) {
      *moment_variance = sd * sd;
    }
    return;
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return;
  }

  const double signed_mean = data.effect_sign * mean;
  if (!std::isfinite(step_log_normalizer_score_space(
        omega, data, signed_mean, sd, sei, omega_stride
      ))) {
    return;
  }

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

  if (!std::isfinite(denom) || !(denom > 0) ||
      !std::isfinite(signed_first) ||
      (moment_variance == NULL &&
       (!std::isfinite(signed_second) || signed_second < 0))) {
    return;
  }

  if (moment_variance != NULL) {
    if (!std::isfinite(centered_first) ||
        !std::isfinite(centered_second) || centered_second < 0) {
      return;
    }
    const double centered_mean = centered_first / denom;
    const double result_mean =
      data.effect_sign * (signed_mean + centered_mean);
    const double result_variance =
      centered_second / denom - centered_mean * centered_mean;
    if (!std::isfinite(result_mean) || !std::isfinite(result_variance) ||
        result_variance < 0) {
      return;
    }
    *moment_mean     = result_mean;
    *moment_variance = result_variance;
    if (std::isfinite(signed_second) && signed_second >= 0) {
      const double result_second = signed_second / denom;
      if (std::isfinite(result_second) && result_second >= 0) {
        *moment_second = result_second;
      }
    }
  } else {
    const double result_mean = data.effect_sign * signed_first / denom;
    const double result_second = signed_second / denom;
    if (!std::isfinite(result_mean) || !std::isfinite(result_second) ||
        result_second < 0) {
      return;
    }
    *moment_mean   = result_mean;
    *moment_second = result_second;
  }
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

  const double upper_score = stable_affine_score(
    z_threshold, sei, mean, sd
  );
  const double lower_score = stable_affine_score(
    -z_threshold, sei, mean, sd
  );

  if (kernel_mode == SELKERNEL_NORMAL) {
    *inverse_weight = 1;
    return pnorm(upper_score, 0, 1, false, false) +
      pnorm(lower_score, 0, 1, true, false);
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    *inverse_weight = 1;
    return pnorm(upper_score, 0, 1, false, false) +
      pnorm(lower_score, 0, 1, true, false);
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Every continuous selected-normal response has |Z| > 0 almost surely.
  // Evaluate this boundary directly instead of summing rounded tail masses.
  if (z_threshold == 0.0) {
    const double log_norm = cpp_selnorm_kernel_log_norm(
      mean, sd, sei, omega, alpha, phack_kind, kernel_mode, data,
      omega_stride, false);
    if (!std::isfinite(log_norm)) return std::numeric_limits<double>::quiet_NaN();
    *inverse_weight = std::exp(-log_norm);
    return 1.0;
  }

  const double signed_mean = data.effect_sign * mean;

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
          signed_mean, sd, sei, omega_stride,
          &upper_prob, &lower_prob, inverse_weight
        )) {
      const double probability = upper_prob + lower_prob;
      return std::isfinite(probability) &&
        probability >= 0 && probability <= 1 ?
        probability : std::numeric_limits<double>::quiet_NaN();
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
        signed_mean, sd, sei, omega_stride, data.trusted_step_partition,
        &upper_prob, &lower_prob, inverse_weight
      )) {
    const double probability = upper_prob + lower_prob;
    return std::isfinite(probability) &&
      probability >= 0 && probability <= 1 ?
      probability : std::numeric_limits<double>::quiet_NaN();
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
  std::vector<double> mass(static_cast<std::size_t>(data.n_bins));
  std::vector<double> lower(static_cast<std::size_t>(data.n_bins));
  std::vector<double> upper(static_cast<std::size_t>(data.n_bins));

  return cpp_selnorm_kernel_rng_workspace(
    mean, sd, sei, omega, u_bin, u_interval, alpha, phack_kind,
    kernel_mode, data, mass.data(), lower.data(), upper.data(),
    omega_stride, validate_omega
  );
}

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

  const double signed_mean = data.effect_sign * mean;

  if (data.trusted_step_partition && data.telescope_probabilities) {
    bool ok = false;
    const double score = step_rng_trusted_partition(
      omega, data, signed_mean, sd, sei, u_bin, u_interval,
      omega_stride, mass, lower, upper, NULL, &ok
    );
    if (ok) {
      return selnorm_fma(data.effect_sign * score, sd, mean);
    }
  }

  const double score = step_rng_logspace(
    omega, data, signed_mean, sd, sei, u_bin, u_interval, omega_stride,
    data.trusted_step_partition, mass, lower, upper
  );

  return selnorm_fma(data.effect_sign * score, sd, mean);
}

bool cpp_selnorm_step_rng_prepare_partition(
  double mean, double sd, double sei, const double *omega,
  const SelNormKernelData &data, double *mass, double *lower, double *upper,
  int *n_groups, double *normalizer, double *log_normalizer)
{
  if (!(sd > 0) || !(sei > 0) || !std::isfinite(mean) ||
      !std::isfinite(sd) || !std::isfinite(sei) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !data.trusted_step_partition || !data.telescope_probabilities) {
    return false;
  }
  bool ok = false;
  *normalizer = step_normalizer_trusted_partition(
    omega, data, data.effect_sign * mean, sd, sei, 1,
    mass, lower, upper, n_groups, &ok
  );
  if (!ok) return false;
  *log_normalizer = std::log(*normalizer);
  return std::isfinite(*log_normalizer);
}

double cpp_selnorm_step_rng_from_partition(
  double mean, double sd, int effect_sign, double u_bin, double u_interval,
  const double *mass, const double *lower, const double *upper,
  int n_groups, double normalizer)
{
  if (!(u_bin >= 0) || !(u_bin <= 1) ||
      !(u_interval >= 0) || !(u_interval <= 1)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double score = step_rng_prepared_partition(
    u_bin, u_interval, mass, lower, upper, n_groups, normalizer
  );
  return selnorm_fma(effect_sign * score, sd, mean);
}

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
  int omega_stride,
  bool validate_omega)
{
  if (!(sd > 0) || !(sei > 0) || !std::isfinite(mean) ||
      !std::isfinite(sd) || !std::isfinite(sei) ||
      !(u_bin >= 0) || !(u_bin <= 1) ||
      !(u_interval >= 0) || !(u_interval <= 1) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !data.trusted_step_partition || !data.telescope_probabilities ||
      (validate_omega &&
       !validate_omega_values(omega, data, omega_stride))) {
    *log_normalizer = std::numeric_limits<double>::quiet_NaN();
    return std::numeric_limits<double>::quiet_NaN();
  }

  *log_normalizer = std::numeric_limits<double>::quiet_NaN();
  const double signed_mean = data.effect_sign * mean;
  bool ok = false;
  const double score = step_rng_trusted_partition(
    omega, data, signed_mean, sd, sei, u_bin, u_interval,
    omega_stride, mass, lower, upper, log_normalizer, &ok
  );
  if (ok && std::isfinite(*log_normalizer)) {
    return selnorm_fma(data.effect_sign * score, sd, mean);
  }

  *log_normalizer = step_log_normalizer_score_space(
    omega, data, signed_mean, sd, sei, omega_stride
  );
  const double fallback_score = step_rng_logspace(
    omega, data, signed_mean, sd, sei, u_bin, u_interval, omega_stride,
    true, mass, lower, upper
  );
  return selnorm_fma(data.effect_sign * fallback_score, sd, mean);
}
