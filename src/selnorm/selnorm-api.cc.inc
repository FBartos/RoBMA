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


bool cpp_selnorm_step_has_transitions(
  const double *omega, const SelNormKernelData &data, int omega_stride)
{
  bool changed = false;
  for (int bin = 0; bin + 1 < data.n_bins; ++bin) {
    changed = changed || (std::isfinite(data.z_lower[bin]) &&
      omega[bin * omega_stride] != omega[(bin + 1) * omega_stride]);
  }
  return changed;
}

bool cpp_selnorm_step_rule_resolves_transitions(
  const double *nodes, int offset, int order,
  const double *mean, const double *residual_sd, const double *loading,
  int dimension, int rank, const double *selection_se,
  const double *omega, const SelNormKernelData &data, int omega_stride,
  double width_scale, const double *axis_width_scales)
{
  const bool changed = cpp_selnorm_step_has_transitions(omega, data, omega_stride);
  if (!changed || rank == 0) return true;
  if (nodes == nullptr || order < 2 || offset < 0 ||
      !(width_scale > 0.0) || !std::isfinite(width_scale)) return false;
  const double *z = nodes + offset;
  // Axes are coupled only through a row that loads on two or more of them; the
  // central node spacing is needed for no other case.
  bool any_coupled = false;
  for (int row = 0; row < dimension && !any_coupled; ++row) {
    int loaded_axes = 0;
    for (int axis = 0; axis < rank; ++axis) {
      loaded_axes += loading[row + dimension * axis] != 0.0;
    }
    any_coupled = loaded_axes >= 2;
  }
  double coupled_gap = 0.0;
  for (int point = 0; point < order; ++point) {
    if (!std::isfinite(z[point]) || (point > 0 && !(z[point] > z[point - 1]))) return false;
    if (any_coupled && point > 0 && z[point - 1] <= 3.0 && z[point] >= -3.0) {
      coupled_gap = std::max(coupled_gap, z[point] - z[point - 1]);
    }
  }
  for (int axis = 0; axis < rank; ++axis) {
    const double axis_width_scale = axis_width_scales == nullptr ?
      width_scale : axis_width_scales[axis];
    if (!(axis_width_scale > 0.0) || !std::isfinite(axis_width_scale)) return false;
    bool coupled = false;
    for (int row = 0; any_coupled && row < dimension; ++row) {
      if (loading[row + dimension * axis] == 0.0) continue;
      for (int other = 0; other < rank; ++other) {
        coupled = coupled || (other != axis && loading[row + dimension * other] != 0.0);
      }
    }
    for (int row = 0; row < dimension; ++row) {
      const double coefficient = loading[row + dimension * axis];
      if (coefficient == 0.0) continue;
      if (!std::isfinite(coefficient) || !(residual_sd[row] > 0.0) ||
          !std::isfinite(residual_sd[row]) || !std::isfinite(mean[row])) return false;
      if (coupled) {
        if (!(coupled_gap > 0.0) ||
            coupled_gap * std::fabs(coefficient) > axis_width_scale * residual_sd[row]) return false;
        continue;
      }
      for (int bin = 0; bin + 1 < data.n_bins; ++bin) {
        if (!std::isfinite(data.z_lower[bin]) ||
            omega[bin * omega_stride] == omega[(bin + 1) * omega_stride]) continue;
        const double center = (coefficient < 0.0 ? -1.0 : 1.0) * stable_affine_score(
          data.effect_sign * data.z_lower[bin], selection_se[row], mean[row], std::fabs(coefficient));
        if (std::isnan(center)) return false;
        // Exterior transitions remain owned by the existing tail gates.
        if (center <= z[0] || center >= z[order - 1]) continue;
        // The nodes are strictly increasing and the center lies strictly
        // inside them: the first node not below it has index 1 .. order - 1.
        const int upper = static_cast<int>(std::lower_bound(z, z + order, center) - z);
        double gap = z[upper] - z[upper - 1];
        if (z[upper] == center && upper + 1 < order) {
          gap = std::max(gap, z[upper + 1] - z[upper]);
        }
        if (gap * std::fabs(coefficient) > axis_width_scale * residual_sd[row]) return false;
      }
    }
  }
  return true;
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
  if (data.n_bins < 1 || omega_stride < 1 || !(sd > 0) || !std::isfinite(sd) ||
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


// Mean and product arrays are contiguous (stride one) across this batch.
// SELNORM_TAIL_MAX_BINS bounds changed boundaries, rather than the bin count.
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
    // Keep the original direct/log routing at the strict normal-minimum floor.
    const bool bounded_tail = bounded && last > std::numeric_limits<double>::min();
    for (std::size_t i = 0; i < count; ++i) {
      if (!std::isfinite(mean[i])) return false;
      const double score = stable_affine_score(
        data.z_lower[0], sei, data.effect_sign * mean[i], sd
      );
      double probability;
      if (bounded_tail && std::isfinite(score) && std::abs(score) >= 8.0) {
        // Q(8) <= phi(8)/8 < 2^-35 < the existing absolute CDF error.
        // This disjoint tail region keeps that same uniform error bound;
        // it need not reproduce the ordinary erfc result bit for bit.
        probability = direction * score > 0.0 ? 0.0 : 1.0;
      } else if (!bounded || !selnorm_cdf_table::evaluate(-direction * score, &probability)) {
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
  std::vector<long double> attempted_product(product, product + count);
  if (step_normalizer_batch_telescoping(mean, count, sd, sei, omega, data,
                                       attempted_product.data())) {
    std::copy(attempted_product.begin(), attempted_product.end(), product);
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

void SelNormDoubleProducts::build(
  const double *residual_sd, const double *selection_se, int rows,
  const double *omega, const SelNormKernelData &data)
{
  usable_ = false;
  if (rows < 1 || !data.trusted_step_partition || !data.telescope_probabilities ||
      data.n_bins < 1 || data.bounded_cdf ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      data.n_bins - 1 > SELNORM_TAIL_MAX_BINS) {
    return;
  }
  double smallest = std::numeric_limits<double>::infinity();
  for (int bin = 0; bin < data.n_bins; ++bin) {
    if (!std::isfinite(omega[bin]) || omega[bin] < 0.0) return;
    if (omega[bin] > 0.0) smallest = std::min(smallest, omega[bin]);
  }
  if (!std::isfinite(smallest) ||
      static_cast<double>(rows) * std::log2(smallest) < -960.0) {
    return;
  }

  // The changed boundaries and weight differences of
  // step_normalizer_batch_telescoping, and per row the inverse residual SD and
  // the scaled cutoffs.
  const double last = omega[data.n_bins - 1];
  if (!std::isfinite(last)) return;
  last_ = last;
  boundaries_ = 0;
  int changed[SELNORM_TAIL_MAX_BINS];
  for (int bin = 0; bin < data.n_bins - 1; ++bin) {
    const double step = omega[bin] - omega[bin + 1];
    if (!std::isfinite(step)) return;
    if (step == 0.0) continue;
    difference_[boundaries_] = step;
    changed[boundaries_] = bin;
    ++boundaries_;
  }
  row_stride_ = 1 + boundaries_;
  rows_.assign(static_cast<std::size_t>(rows) * row_stride_, 0.0);
  for (int row = 0; row < rows; ++row) {
    const double sd = residual_sd[row];
    const double sei = selection_se[row];
    if (!(sd > 0) || !std::isfinite(sd) || !(sei > 0) || !std::isfinite(sei)) return;
    const double inverse_sd = data.effect_sign / sd;
    if (!std::isfinite(inverse_sd)) return;
    double *constants = rows_.data() + static_cast<std::size_t>(row) * row_stride_;
    constants[0] = inverse_sd;
    for (int k = 0; k < boundaries_; ++k) {
      const double scaled = data.z_lower[changed[k]] * sei / sd;
      if (!std::isfinite(scaled)) return;
      constants[1 + k] = scaled;
    }
  }
  usable_ = true;
}

bool SelNormDoubleProducts::multiply(
  int row, const double *mean, std::size_t count, double *product,
  double *normalizer) const
{
  const double *constants = rows_.data() + static_cast<std::size_t>(row) * row_stride_;
  const double inverse_sd = constants[0];
  const double last = last_;
  const double absolute_last = std::abs(last);
  const std::size_t block = 256;
  double tail[block];
  double magnitude[block];
  // The normalizers of the changed boundaries, block by block into the
  // caller's workspace, each block checked (those of the batch routine and a
  // finite location) as soon as it is complete.
  for (std::size_t start = 0; start < count; start += block) {
    const std::size_t size = std::min(block, count - start);
    const double *location = mean + start;
    double *value = normalizer + start;
    if (boundaries_ == 1) {
      if (!selnorm_tail::single_boundary_normalizers(location, size, constants[1],
            inverse_sd, last, difference_[0], value, tail, magnitude)) return false;
      continue;
    }
    std::fill(value, value + size, last);
    std::fill(magnitude, magnitude + size, absolute_last);
    for (int bin = 0; bin < boundaries_; ++bin) {
      selnorm_tail::upper_tail_affine(location, size, constants[1 + bin], inverse_sd, tail);
      selnorm_tail::accumulate_boundary(tail, size, difference_[bin], value, magnitude);
    }
    if (!selnorm_tail::validate_normalizers(location, value, magnitude, size)) return false;
  }
  return selnorm_tail::multiply_products(normalizer, count, product);
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
  double q_cdf;
  double q_tail;
  cdf_tail_at(q_score, &q_cdf, &q_tail);
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

namespace {

// One standard normal score with both tails and its density.
struct SharedNormalPoint {
  double score;
  double cdf;
  double tail;
  double phi;
  double score_phi;
};

SharedNormalPoint shared_normal_point(double score)
{
  SharedNormalPoint point;
  point.score = score;
  cdf_tail_at(score, &point.cdf, &point.tail);
  point.phi = std_phi(score);
  point.score_phi = score * point.phi;
  if (!std::isfinite(point.score_phi)) {
    point.score_phi = 0;
  }
  return point;
}

// Standard normal probability of [lower, upper] from the endpoints' shared
// tails. An interval inside either side of zero subtracts its smaller tails;
// one that touches or contains zero uses erf, whose terms do not cancel
// (as interval_log_prob() does). `operand` receives the magnitude of the
// subtracted values, which bounds the rounding error of the result. NaN marks
// a difference that rounding left negative.
double shared_interval_probability(const SharedNormalPoint &lower,
                                   const SharedNormalPoint &upper,
                                   double *operand)
{
  *operand = 0;
  if (!(lower.score < upper.score)) {
    return 0;
  }
  double probability;
  if (lower.score > 0) {
    probability = lower.tail - upper.tail;
    *operand    = lower.tail;
  } else if (upper.score < 0) {
    probability = upper.cdf - lower.cdf;
    *operand    = upper.cdf;
  } else {
    const double upper_erf = std::isinf(upper.score) ? 1.0 :
      std::erf(upper.score * M_SQRT1_2);
    const double lower_erf = std::isinf(lower.score) ? -1.0 :
      std::erf(lower.score * M_SQRT1_2);
    probability = 0.5 * (upper_erf - lower_erf);
    *operand    = probability;
  }
  return probability >= 0 && std::isfinite(probability) ?
    probability : std::numeric_limits<double>::quiet_NaN();
}

// Normalized paired log tails from the two tail masses, without forming a
// ratio that could overflow: each tail is taken relative to the larger mass.
void shared_log_tail_pair(double lower_mass, double upper_mass,
                          double *log_lower, double *log_upper)
{
  if (lower_mass == 0) {
    *log_lower = -INFINITY;
    *log_upper = 0;
  } else if (upper_mass == 0) {
    *log_lower = 0;
    *log_upper = -INFINITY;
  } else if (upper_mass <= lower_mass) {
    const double ratio = upper_mass / lower_mass;
    *log_lower = -std::log1p(ratio);
    *log_upper = std::log(upper_mass) - std::log(lower_mass) - std::log1p(ratio);
  } else {
    const double ratio = lower_mass / upper_mass;
    *log_lower = std::log(lower_mass) - std::log(upper_mass) - std::log1p(ratio);
    *log_upper = -std::log1p(ratio);
  }
}

// Log of a bound on the absolute error that underflow leaves in the terms one
// boundary at standardized score `score` contributes to a bin: its smaller
// tail (used by the bin and by a query slice), density, and score times
// density. Below |score| = 37 these are in the normal range. Beyond it each
// is at most |score| phi(score) in truth, and R's pnorm_both() and dnorm()
// return subnormal values only below |score| = 39 and exact zeros beyond.
double log_boundary_underflow_error(double score)
{
  const double magnitude = std::abs(score);
  if (!(magnitude > 37) || std::isinf(magnitude)) {
    return -INFINITY;
  }
  const double log_sqrt_2pi = 0.918938533204672741780329736406;
  const double log_true = std::log(magnitude) - 0.5 * magnitude * magnitude -
    log_sqrt_2pi;
  const double log_computed = magnitude < 39 ?
    std::log(std::numeric_limits<double>::min()) : -INFINITY;
  return std::log(4.0) + add_logspace(log_true, log_computed);
}

// The shared step evaluation in the signed direction. Returns false when the
// tails need the log-space kernel: a nonpositive or unavailable mass, a tail
// whose mass is close to the underflow range of probability space, or a
// weighted bin term lost to underflow that could reach a tail mass's last
// digit. `moments_unavailable` is set when such a loss could reach the
// normalizer's last digit: the probability-space moment kernels drop the same
// terms, so these moments are not representable. `moments_resolved` is false
// when the centered moments lose more than four digits to cancellation (a
// narrow interval carrying the mass, relative to the operands it is formed
// from); the partition-moment kernel then supplies them.
bool step_predictive_summary_shared(const double *omega,
                                    const SelNormKernelData &data,
                                    double q_score, double signed_mean,
                                    double sd, double sei, int omega_stride,
                                    double *log_lower, double *log_upper,
                                    double *centered_mean,
                                    double *centered_variance,
                                    bool *moments_resolved,
                                    bool *moments_unavailable)
{
  // Masses below this stay in probability space only far from underflow.
  const double minimum_mass = 1e-280;
  // Largest admitted ratio of the moments' operand scale to the variance.
  const double maximum_cancellation = 1e4;
  *moments_resolved    = false;
  *moments_unavailable = false;
  const SharedNormalPoint query = shared_normal_point(q_score);
  const double log_query_error = log_boundary_underflow_error(q_score);
  SharedNormalPoint upper = shared_normal_point(INFINITY);
  double log_underflow_error = -INFINITY;
  double lower_mass   = 0;
  double upper_mass   = 0;
  double first        = 0;
  double second       = 0;
  double first_scale  = 0;
  double second_scale = 0;
  bool lower_support  = false;
  bool upper_support  = false;

  for (int b = 0; b < data.n_bins; ) {
    const double w = omega_at(omega, b + 1, omega_stride);
    int last = b;
    while (last + 1 < data.n_bins &&
           omega_at(omega, last + 2, omega_stride) == w) {
      ++last;
    }
    if (!(w >= 0) || !std::isfinite(w)) {
      return false;
    }

    const SharedNormalPoint lower = last + 1 == data.n_bins ?
      shared_normal_point(-INFINITY) :
      shared_normal_point(stable_affine_score(
        data.z_lower[last], sei, signed_mean, sd
      ));
    if (std::isnan(lower.score) || !(lower.score < upper.score)) {
      return false;
    }

    if (w > 0) {
      double operand;
      double unused;
      const double probability = shared_interval_probability(
        lower, upper, &operand
      );
      if (std::isnan(probability)) {
        return false;
      }
      double log_bin_error = add_logspace(
        log_boundary_underflow_error(lower.score),
        log_boundary_underflow_error(upper.score)
      );
      if (lower.score < query.score && query.score < upper.score) {
        log_bin_error = add_logspace(log_bin_error, log_query_error);
      }
      if (log_bin_error != -INFINITY) {
        log_underflow_error = add_logspace(
          log_underflow_error, std::log(w) + log_bin_error
        );
      }
      if (query.score <= lower.score) {
        upper_mass   += w * probability;
        upper_support = true;
      } else if (query.score >= upper.score) {
        lower_mass   += w * probability;
        lower_support = true;
      } else {
        const double below = shared_interval_probability(lower, query, &unused);
        const double above = shared_interval_probability(query, upper, &unused);
        if (std::isnan(below) || std::isnan(above)) {
          return false;
        }
        lower_mass   += w * below;
        upper_mass   += w * above;
        lower_support = true;
        upper_support = true;
      }
      first        += w * (lower.phi - upper.phi);
      second       += w * (probability + lower.score_phi - upper.score_phi);
      first_scale  += w * (lower.phi + upper.phi);
      second_scale += w * (operand + std::abs(lower.score_phi) +
                           std::abs(upper.score_phi));
    }

    upper = lower;
    b = last + 1;
  }

  // A bin term lost to underflow can still carry a tail when its weight is
  // large enough (weights spanning beyond double precision); the tails are
  // then evaluated in log space. Weights bounded by one never reach this
  // with masses above `minimum_mass`.
  const double log_epsilon = std::log(std::numeric_limits<double>::epsilon());
  const auto underflow_material = [&](double mass) {
    return log_underflow_error != -INFINITY &&
      !(log_underflow_error <= log_epsilon + std::log(mass));
  };
  *moments_unavailable = underflow_material(lower_mass + upper_mass);
  if ((lower_support && !(lower_mass >= minimum_mass)) ||
      (upper_support && !(upper_mass >= minimum_mass)) ||
      !std::isfinite(lower_mass) || !std::isfinite(upper_mass) ||
      !std::isfinite(first) || !std::isfinite(second) ||
      !std::isfinite(first_scale) || !std::isfinite(second_scale) ||
      (lower_support && underflow_material(lower_mass)) ||
      (upper_support && underflow_material(upper_mass))) {
    return false;
  }

  const double normalizer = lower_mass + upper_mass;
  if (!(normalizer >= minimum_mass) || !std::isfinite(normalizer)) {
    return false;
  }
  shared_log_tail_pair(lower_mass, upper_mass, log_lower, log_upper);
  if (std::isnan(*log_lower) || std::isnan(*log_upper)) {
    return false;
  }

  const double mean_score     = first / normalizer;
  const double variance_score = second / normalizer - mean_score * mean_score;
  const double mean_scale     = first_scale / normalizer;
  const double variance_scale = second_scale / normalizer +
    2 * std::abs(mean_score) * mean_scale + mean_score * mean_score;
  *moments_resolved = variance_score > 0 &&
    variance_scale <= maximum_cancellation * variance_score &&
    mean_scale * mean_scale <=
      maximum_cancellation * maximum_cancellation * variance_score;
  *centered_mean     = sd * mean_score;
  *centered_variance = sd * sd * variance_score;
  return true;
}

}  // namespace

void cpp_selnorm_kernel_predictive_summary(
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
  double *moment_mean,
  double *moment_variance,
  int omega_stride,
  bool *moments_unavailable)
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  bool unrepresentable_moments = false;
  if (moments_unavailable != nullptr) *moments_unavailable = false;
  if (kernel_mode == SELKERNEL_STEP && data.trusted_step_partition &&
      sd > 0 && std::isfinite(sd) && sei > 0 && std::isfinite(sei) &&
      std::isfinite(mean) && std::isfinite(q) &&
      (data.effect_sign == 1 || data.effect_sign == -1) &&
      valid_phack_kind(phack_kind) && valid_phack_kind(data.q)) {
    const double signed_mean = data.effect_sign * mean;
    const double q_score = stable_normal_score(
      data.effect_sign * q, signed_mean, sd
    );
    double signed_log_lower;
    double signed_log_upper;
    double centered_mean;
    double centered_variance;
    bool moments_resolved;
    // A cell whose shared moments are ill-conditioned has a narrow interval
    // carrying its mass, which conditions its tails as well; it uses the
    // separate kernels below for all four values.
    if (std::isfinite(q_score) &&
        step_predictive_summary_shared(
          omega, data, q_score, signed_mean, sd, sei, omega_stride,
          &signed_log_lower, &signed_log_upper,
          &centered_mean, &centered_variance, &moments_resolved,
          &unrepresentable_moments
        ) && moments_resolved) {
      const double result_mean = data.effect_sign * (signed_mean + centered_mean);
      if (std::isfinite(result_mean) && std::isfinite(centered_variance) &&
          centered_variance >= 0) {
        *log_lower       = data.effect_sign == 1 ? signed_log_lower : signed_log_upper;
        *log_upper       = data.effect_sign == 1 ? signed_log_upper : signed_log_lower;
        *moment_mean     = result_mean;
        *moment_variance = centered_variance;
        return;
      }
    }
  }

  cpp_selnorm_kernel_log_tail_pair(
    q, mean, sd, sei, omega, alpha, phack_kind, kernel_mode, data,
    log_lower, log_upper, omega_stride, false
  );
  if (unrepresentable_moments) {
    // The partition-moment kernel would drop the same underflowed terms.
    // These moments lie beyond double precision; refuse them.
    *moment_mean     = nan;
    *moment_variance = nan;
    if (moments_unavailable != nullptr) *moments_unavailable = true;
    return;
  }
  double moment_second = nan;
  cpp_selnorm_kernel_moments(
    mean, sd, sei, omega, alpha, phack_kind, kernel_mode, data,
    moment_mean, &moment_second, omega_stride, false, moment_variance
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

  // A near-unit selected tail probability is determined by its central
  // complement. Retain the exact normalizer used by the successful tail pair.
  const auto finish_probability = [&](double upper_prob, double lower_prob,
                                      double log_normalizer) {
    const double tails = upper_prob + lower_prob;
    if (!std::isfinite(tails) || tails < 0.0 || !std::isfinite(log_normalizer)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (tails <= .5) return tails;
    const double central_lower = stable_affine_score(-z_threshold, sei, signed_mean, sd);
    const double central_upper = stable_affine_score(z_threshold, sei, signed_mean, sd);
    const double log_central = step_log_normalizer_score_space(
      omega, data, signed_mean, sd, sei, omega_stride, central_lower, central_upper);
    if (log_central == -INFINITY) return 1.0;
    if (!std::isfinite(log_central) || log_central > log_normalizer) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double probability = -std::expm1(log_central - log_normalizer);
    return std::isfinite(probability) && probability >= 0.0 && probability <= 1.0 ?
      probability : std::numeric_limits<double>::quiet_NaN();
  };

  if (data.trusted_step_partition && data.telescope_probabilities) {
    double upper_prob;
    double lower_prob;
    double log_normalizer;
    const double q_upper_z = data.effect_sign * z_threshold;
    const double q_lower_z = -q_upper_z;
    const bool upper_signed_lower_tail = data.effect_sign == 1 ? false : true;
    const bool lower_signed_lower_tail = data.effect_sign == 1 ? true : false;

    if (step_cdf_pair_trusted_telescoping(
          omega, data,
          q_upper_z, upper_signed_lower_tail,
          q_lower_z, lower_signed_lower_tail,
          signed_mean, sd, sei, omega_stride,
          &upper_prob, &lower_prob, inverse_weight, &log_normalizer
        )) {
      return finish_probability(upper_prob, lower_prob, log_normalizer);
    }
  }

  double upper_prob;
  double lower_prob;
  double log_normalizer;
  if (step_cdf_pair_logspace(
        omega, data,
        data.effect_sign * z_threshold,
        data.effect_sign == 1 ? false : true,
        -data.effect_sign * z_threshold,
        data.effect_sign == 1 ? true : false,
        signed_mean, sd, sei, omega_stride, data.trusted_step_partition,
        &upper_prob, &lower_prob, inverse_weight, &log_normalizer
      )) {
    return finish_probability(upper_prob, lower_prob, log_normalizer);
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
