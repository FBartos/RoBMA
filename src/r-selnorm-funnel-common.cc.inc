// Funnel summaries evaluate the posterior-row mixture distribution directly.
// Zero total SD creates an exact atom, active weight functions use the selected-
// normal CDF, and unweighted rows use the ordinary normal CDF. Quantiles retain
// generalized-inverse semantics when atoms or selection gaps are present.

int z_step_bin(double z, const SelNormKernelData &data)
{
  for (int b = 0; b < data.n_bins; ++b) {
    if (z >= data.z_lower[b] && z <= data.z_upper[b]) {
      return b + 1;
    }
  }

  return z > data.z_upper[0] ? 1 : data.n_bins;
}

bool row_has_active_phack(int mode, double alpha, int phack_kind)
{
  return (mode == SELKERNEL_PHACK_POWER ||
    mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
}

double validated_probability_or_na(double x)
{
  // Never repair an out-of-range result: surface invalid kernel output as NA.
  return std::isfinite(x) && x >= 0 && x <= 1 ? x : NA_REAL;
}

double funnel_interval_prob(double lower, double upper, double mean, double sd)
{
  if (!(sd > 0) || lower >= upper) {
    return 0;
  }

  if (std::isinf(lower) && lower < 0) {
    if (std::isinf(upper) && upper > 0) {
      return 1;
    }
    return pnorm(upper, mean, sd, true, false);
  }
  if (std::isinf(upper) && upper > 0) {
    return pnorm(lower, mean, sd, false, false);
  }

  if (lower > mean) {
    return pnorm(lower, mean, sd, false, false) -
      pnorm(upper, mean, sd, false, false);
  }

  return pnorm(upper, mean, sd, true, false) -
    pnorm(lower, mean, sd, true, false);
}

double funnel_empirical_quantile(const std::vector<double> &x, double p)
{
  const int n = static_cast<int>(x.size());
  std::vector<double> sorted(x);

  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(sorted[static_cast<size_t>(i)])) {
      return NA_REAL;
    }
  }
  std::sort(sorted.begin(), sorted.end());

  const int index = std::max(
    0,
    std::min(n - 1, static_cast<int>(std::ceil(static_cast<double>(n) * p)) - 1)
  );
  return sorted[static_cast<size_t>(index)];
}

struct FunnelSeContext {
  int S;
  double se;
  int direction;
  const double *mu;
  const double *tau;
  const double *PET;
  const double *PEESE;
  const int *is_weightfunction;
  const double *omega;
  const std::vector<double> *row_alpha;
  const std::vector<int> *row_phack;
  const std::vector<int> *row_mode;
  const std::vector<bool> *row_active_phack;
  const SelNormKernelData *data;
  std::vector<double> location;
  std::vector<double> total_sd;
  std::vector<char> zero_sd;
  std::vector<char> weighted;
  std::vector<char> step_plan_valid;
  std::vector<int> step_plan_groups;
  std::vector<double> step_lower_score;
  std::vector<double> step_upper_score;
  std::vector<double> step_log_weight;
};

double funnel_selected_cdf_zero_se(double q, int s, const FunnelSeContext &ctx)
{
  const SelNormKernelData &data = *ctx.data;
  const double *omega_s = ctx.omega + s;
  const double mean = data.effect_sign * ctx.mu[s];
  const double sd = ctx.total_sd[static_cast<size_t>(s)];
  const double q_signed = data.effect_sign * q;
  double requested_mass  = 0;
  double complement_mass = 0;

  for (int b = 0; b < data.n_bins; ++b) {
    const double z_lower = data.z_lower[b];
    const double z_upper = data.z_upper[b];
    double lower;
    double upper;

    if (std::isinf(z_lower) && z_lower < 0 &&
        std::isinf(z_upper) && z_upper > 0) {
      lower = -INFINITY;
      upper =  INFINITY;
    } else if (std::isinf(z_upper) && z_upper > 0) {
      lower = 0;
      upper = INFINITY;
    } else if (std::isinf(z_lower) && z_lower < 0) {
      lower = -INFINITY;
      upper = 0;
    } else {
      continue;
    }

    const double weight = omega_s[static_cast<size_t>(ctx.S) * b];
    if (data.effect_sign == 1) {
      requested_mass += weight * funnel_interval_prob(
        lower, std::min(upper, q_signed), mean, sd
      );
      complement_mass += weight * funnel_interval_prob(
        std::max(lower, q_signed), upper, mean, sd
      );
    } else {
      requested_mass += weight * funnel_interval_prob(
        std::max(lower, q_signed), upper, mean, sd
      );
      complement_mass += weight * funnel_interval_prob(
        lower, std::min(upper, q_signed), mean, sd
      );
    }
  }

  const double total_mass = requested_mass + complement_mass;
  if (!(total_mass > 0) || !std::isfinite(total_mass)) {
    return NA_REAL;
  }

  return validated_probability_or_na(requested_mass / total_mass);
}

void funnel_prepare_se_context(FunnelSeContext *ctx)
{
  const double se2 = ctx->se * ctx->se;

  ctx->location.assign(static_cast<size_t>(ctx->S), 0);
  ctx->total_sd.assign(static_cast<size_t>(ctx->S), 0);
  ctx->zero_sd.assign(static_cast<size_t>(ctx->S), 0);
  ctx->weighted.assign(static_cast<size_t>(ctx->S), 0);
  ctx->step_plan_valid.assign(static_cast<size_t>(ctx->S), 0);
  ctx->step_plan_groups.assign(static_cast<size_t>(ctx->S), 0);
  const int n_bins = ctx->data->n_bins;
  const size_t plan_size = static_cast<size_t>(ctx->S) *
    static_cast<size_t>(n_bins);
  ctx->step_lower_score.assign(plan_size, 0);
  ctx->step_upper_score.assign(plan_size, 0);
  ctx->step_log_weight.assign(plan_size, 0);

  for (int s = 0; s < ctx->S; ++s) {
    const size_t si = static_cast<size_t>(s);
    const double location = ctx->mu[s] +
      ctx->direction * ctx->PET[s] * ctx->se +
      ctx->direction * ctx->PEESE[s] * se2;
    const double total_sd = std::hypot(ctx->se, ctx->tau[s]);

    ctx->location[si] = location;
    ctx->total_sd[si] = total_sd;
    ctx->zero_sd[si]  = total_sd == 0;
    ctx->weighted[si] = ctx->is_weightfunction[s] != 0 && !ctx->zero_sd[si];

    const int mode = (*ctx->row_mode)[si];
    if (ctx->weighted[si] && ctx->se > 0 &&
        !(*ctx->row_active_phack)[si] &&
        (mode == SELKERNEL_STEP || mode == SELKERNEL_STEP_PHACK_POWER)) {
      const size_t offset = si * static_cast<size_t>(n_bins);
      int n_groups = 0;
      ctx->step_plan_valid[si] = cpp_selnorm_step_cdf_plan(
        location,
        total_sd,
        ctx->se,
        ctx->omega + s,
        *ctx->data,
        ctx->step_lower_score.data() + offset,
        ctx->step_upper_score.data() + offset,
        ctx->step_log_weight.data() + offset,
        &n_groups,
        ctx->S,
        false
      );
      ctx->step_plan_groups[si] = n_groups;
    }

  }
}

double funnel_model_averaged_cdf_native(double q, const FunnelSeContext &ctx)
{
  double cdf_sum = 0;

  for (int s = 0; s < ctx.S; ++s) {
    const size_t si = static_cast<size_t>(s);
    double cdf_s;

    if (ctx.zero_sd[si]) {
      cdf_s = q >= ctx.location[si] ? 1.0 : 0.0;
    } else if (ctx.weighted[si]) {
      if (ctx.se <= 0) {
        cdf_s = funnel_selected_cdf_zero_se(q, s, ctx);
      } else if ((*ctx.row_active_phack)[si]) {
        cdf_s = NA_REAL;
      } else if (ctx.step_plan_valid[si]) {
        const size_t offset = si * static_cast<size_t>(ctx.data->n_bins);
        cdf_s = cpp_selnorm_step_cdf_from_plan(
          q,
          ctx.location[si],
          ctx.total_sd[si],
          ctx.se,
          *ctx.data,
          ctx.step_lower_score.data() + offset,
          ctx.step_upper_score.data() + offset,
          ctx.step_log_weight.data() + offset,
          ctx.step_plan_groups[si],
          true
        );
      } else {
        cdf_s = cpp_selnorm_kernel_cdf(
          q,
          ctx.location[si],
          ctx.total_sd[si],
          ctx.se,
          ctx.omega + s,
          (*ctx.row_alpha)[si],
          (*ctx.row_phack)[si],
          (*ctx.row_mode)[si],
          *ctx.data,
          ctx.S,
          true,
          false
        );
      }
    } else {
      cdf_s = pnorm(q, ctx.location[si], ctx.total_sd[si], true, false);
    }

    cdf_sum += validated_probability_or_na(cdf_s);
  }

  return cdf_sum / static_cast<double>(ctx.S);
}

bool funnel_mixture_has_full_support(const FunnelSeContext &ctx)
{
  for (int s = 0; s < ctx.S; ++s) {
    const size_t si = static_cast<size_t>(s);
    if (!ctx.weighted[si]) {
      return true;
    }

    const int mode = (*ctx.row_mode)[si];
    if (mode == SELKERNEL_NORMAL ||
        (mode == SELKERNEL_PHACK_POWER && !(*ctx.row_active_phack)[si])) {
      return true;
    }
  }

  for (int b = 0; b < ctx.data->n_bins; ++b) {
    bool bin_has_support = false;
    for (int s = 0; s < ctx.S; ++s) {
      const size_t si = static_cast<size_t>(s);
      const int mode = (*ctx.row_mode)[si];
      if ((mode == SELKERNEL_STEP || mode == SELKERNEL_STEP_PHACK_POWER) &&
          ctx.omega[s + ctx.S * b] > 0) {
        bin_has_support = true;
        break;
      }
    }
    if (!bin_has_support) {
      return false;
    }
  }

  return true;
}

double funnel_model_averaged_quantile_native(double p,
                                             const FunnelSeContext &ctx,
                                             double previous,
                                             bool use_previous)
{
  bool all_zero_sd = true;
  bool all_positive_sd = true;
  double lower = INFINITY;
  double upper = -INFINITY;
  double step = 0;

  for (int s = 0; s < ctx.S; ++s) {
    const size_t si = static_cast<size_t>(s);
    all_zero_sd = all_zero_sd && ctx.zero_sd[si];
    all_positive_sd = all_positive_sd && !ctx.zero_sd[si];
    const double spread = ctx.total_sd[si];
    lower = std::min(lower, ctx.location[si] - 10 * spread);
    upper = std::max(upper, ctx.location[si] + 10 * spread);
    step  = std::max(step, spread);

  }

  if (all_zero_sd) {
    return funnel_empirical_quantile(ctx.location, p);
  }

  if (!std::isfinite(lower) || !std::isfinite(upper) || lower >= upper) {
    return NA_REAL;
  }

  const double global_lower = lower;
  const double global_upper = upper;
  const double global_step  = step;
  bool using_previous = use_previous && std::isfinite(previous);

  if (using_previous) {
    // The previous quantile only narrows the initial bracket. Failure retries
    // the global bracket, and the root/generalized inverse itself is unchanged.
    step = std::max(1e-3, 0.05 * global_step);
    lower = previous - step;
    upper = previous + step;
  }

  double lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
  double upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;

  if (!std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    if (!using_previous) {
      return NA_REAL;
    }
    using_previous = false;
    lower = global_lower;
    upper = global_upper;
    step = global_step;
    lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
    upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;
    if (!std::isfinite(lower_value) || !std::isfinite(upper_value)) {
      return NA_REAL;
    }
  }

  for (int i = 0; i < 25; ++i) {
    if (lower_value < 0 && upper_value >= 0) {
      break;
    }
    if (lower_value >= 0) {
      lower -= step;
      lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
    }
    if (upper_value < 0) {
      upper += step;
      upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;
    }
    step *= 2;
  }

  if (lower_value >= 0 || upper_value < 0 ||
      !std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    if (using_previous) {
      lower = global_lower;
      upper = global_upper;
      step = global_step;
      lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
      upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;

      for (int i = 0; i < 25; ++i) {
        if (lower_value < 0 && upper_value >= 0) {
          break;
        }
        if (lower_value >= 0) {
          lower -= step;
          lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
        }
        if (upper_value < 0) {
          upper += step;
          upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;
        }
        step *= 2;
      }
    }
  }

  if (lower_value >= 0 || upper_value < 0 ||
      !std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    return NA_REAL;
  }

  // Full-support continuous mixtures have a unique root and do not need exact
  // generalized-inverse bisection. Atoms and selection gaps retain the exact
  // path.
  if (all_positive_sd && funnel_mixture_has_full_support(ctx)) {
    const double tolerance = RoBMA::plot_root_tolerance(
      global_lower,
      global_upper,
      global_step
    );
    const double root = RoBMA::plot_brent_root(
      [&ctx, p](double q) {
        return funnel_model_averaged_cdf_native(q, ctx) - p;
      },
      lower,
      upper,
      lower_value,
      upper_value,
      tolerance
    );
    if (std::isfinite(root)) {
      return root;
    }
  }

  while (true) {
    const double mid = 0.5 * lower + 0.5 * upper;
    if (mid <= lower || mid >= upper) {
      break;
    }
    const double mid_value = funnel_model_averaged_cdf_native(mid, ctx) - p;

    if (!std::isfinite(mid_value)) {
      return NA_REAL;
    }

    if (mid_value >= 0) {
      upper = mid;
    } else {
      lower = mid;
    }

  }

  return upper;
}
