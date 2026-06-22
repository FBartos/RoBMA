namespace {

int matrix_nrow(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[0];
}

int matrix_ncol(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[1];
}

void check_real(SEXP x, const char *name)
{
  if (TYPEOF(x) != REALSXP) {
    Rf_error("'%s' must be numeric.", name);
  }
}

void check_integer(SEXP x, const char *name)
{
  if (TYPEOF(x) != INTSXP) {
    Rf_error("'%s' must be integer.", name);
  }
}

int scalar_or_row_int(const int *x, int n, int row, const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (row >= n) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
  return x[row];
}

double scalar_or_row_real(const double *x, int n, int row, const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (row >= n) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
  return x[row];
}

void validate_matrix_dims(SEXP x, const char *name, int S, int K)
{
  if (matrix_nrow(x, name) != S || matrix_ncol(x, name) != K) {
    Rf_error("'%s' dimensions must match 'mu_num'.", name);
  }
}

void validate_matrix_dims_named(SEXP x, const char *name, int S, int K,
                                const char *ref)
{
  if (matrix_nrow(x, name) != S || matrix_ncol(x, name) != K) {
    Rf_error("'%s' dimensions must match '%s'.", name, ref);
  }
}

void check_logical(SEXP x, const char *name)
{
  if (TYPEOF(x) != LGLSXP || Rf_length(x) != 1) {
    Rf_error("'%s' must be a single logical value.", name);
  }
}

double scalar_or_col_real(const double *x, int n, int col, const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (col >= n) {
    Rf_error("'%s' must have length 1 or one value per observation.", name);
  }
  return x[col];
}

void validate_omega_matrix(const double *omega, int S, int B)
{
  for (int b = 0; b < B; ++b) {
    for (int s = 0; s < S; ++s) {
      const double w = omega[s + S * b];
      if (!std::isfinite(w) || w < 0) {
        Rf_error("'omega' must contain finite non-negative values.");
      }
    }
  }
}

void set_kernel_data(SelNormKernelData *data, int B, int n_segments,
                     SEXP sign, SEXP q, SEXP z_lower, SEXP z_upper,
                     SEXP phack_z_source, SEXP phack_z_dest,
                     SEXP segment_bounds, SEXP segment_step_bin,
                     SEXP segment_phack_region,
                     bool telescope_probabilities)
{
  data->n_bins               = B;
  data->n_segments           = n_segments;
  data->effect_sign          = INTEGER(sign)[0];
  data->q                    = INTEGER(q)[0];
  data->z_lower              = REAL(z_lower);
  data->z_upper              = REAL(z_upper);
  data->phack_z_source       = REAL(phack_z_source);
  data->phack_z_dest         = REAL(phack_z_dest);
  data->segment_bounds       = REAL(segment_bounds);
  data->segment_step_bin     = INTEGER(segment_step_bin);
  data->segment_phack_region = INTEGER(segment_phack_region);
  data->segment_step_bin_real = 0;
  data->segment_phack_region_real = 0;
  data->trusted_step_partition = selnorm_is_descending_step_partition(
    data->z_lower, data->z_upper, B
  );
  data->telescope_probabilities = data->trusted_step_partition &&
    telescope_probabilities;
}

struct SelNormNativeArgs {
  int S;
  int B;
  int n_segments;
  int alpha_len;
  int phack_kind_len;
  int kernel_mode_len;
  const double *omega;
  const double *alpha;
  const int *phack_kind;
  const int *kernel_mode;
  SelNormKernelData data;

  const double *omega_row(int row) const
  {
    return omega + row;
  }

  double alpha_at(int row) const
  {
    return scalar_or_row_real(alpha, alpha_len, row, "alpha");
  }

  int phack_kind_at(int row) const
  {
    return scalar_or_row_int(phack_kind, phack_kind_len, row, "phack_kind");
  }

  int kernel_mode_at(int row) const
  {
    return scalar_or_row_int(kernel_mode, kernel_mode_len, row, "kernel_mode");
  }
};

void validate_selnorm_row_arg_length(SEXP x, int S, const char *name)
{
  if (Rf_length(x) != 1 && Rf_length(x) != S) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
}

SelNormNativeArgs validate_selnorm_native_args(
    SEXP omega, SEXP alpha, SEXP phack_kind, SEXP kernel_mode,
    SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
    SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
    SEXP segment_step_bin, SEXP segment_phack_region,
    SEXP telescope_probabilities, int S)
{
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int B = matrix_ncol(omega, "omega");
  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  validate_selnorm_row_arg_length(alpha, S, "alpha");
  validate_selnorm_row_arg_length(phack_kind, S, "phack_kind");
  validate_selnorm_row_arg_length(kernel_mode, S, "kernel_mode");

  SelNormNativeArgs args;
  args.S               = S;
  args.B               = B;
  args.n_segments      = n_segments;
  args.alpha_len       = Rf_length(alpha);
  args.phack_kind_len  = Rf_length(phack_kind);
  args.kernel_mode_len = Rf_length(kernel_mode);
  args.omega           = REAL(omega);
  args.alpha           = REAL(alpha);
  args.phack_kind      = INTEGER(phack_kind);
  args.kernel_mode     = INTEGER(kernel_mode);

  validate_omega_matrix(args.omega, S, B);
  set_kernel_data(
    &args.data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  return args;
}

int z_step_bin(double z, const SelNormKernelData &data)
{
  for (int b = 0; b < data.n_bins; ++b) {
    if (z >= data.z_lower[b] - 1e-12 && z <= data.z_upper[b] + 1e-12) {
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

double clamp_probability_value(double x)
{
  if (x < 0) {
    return 0;
  }
  if (x > 1) {
    return 1;
  }
  return x;
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

  double out = pnorm(upper, mean, sd, true, false) -
    pnorm(lower, mean, sd, true, false);
  if (out < 0 && out > -1e-14) {
    out = 0;
  }
  if (out > 1 && out < 1 + 1e-14) {
    out = 1;
  }
  return out;
}

double funnel_quantile_type8(const std::vector<double> &x, double p)
{
  const int n = static_cast<int>(x.size());
  std::vector<double> sorted(x);

  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(sorted[static_cast<size_t>(i)])) {
      return NA_REAL;
    }
  }
  std::sort(sorted.begin(), sorted.end());

  const double h = (static_cast<double>(n) + 1.0 / 3.0) * p + 1.0 / 3.0;
  if (h <= 1) {
    return sorted[0];
  }
  if (h >= n) {
    return sorted[static_cast<size_t>(n - 1)];
  }

  const int j = static_cast<int>(std::floor(h));
  const double g = h - static_cast<double>(j);
  return (1 - g) * sorted[static_cast<size_t>(j - 1)] +
    g * sorted[static_cast<size_t>(j)];
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
  std::vector<double> mean_z;
  std::vector<double> sd_z;
  std::vector<double> log_norm;
  std::vector<double> normalizer;
  std::vector<double> boundary_tail;
  std::vector<char> zero_sd;
  std::vector<char> weighted;
  std::vector<char> fast_step;
};

double funnel_selected_cdf_zero_se(double q, int s, const FunnelSeContext &ctx)
{
  const SelNormKernelData &data = *ctx.data;
  const double *omega_s = ctx.omega + s;
  const double mean = data.effect_sign * ctx.mu[s];
  const double sd = ctx.total_sd[static_cast<size_t>(s)];
  const double q_signed = data.effect_sign * q;
  double denom = 0;
  double mass = 0;

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

    const double bin_mass = funnel_interval_prob(lower, upper, mean, sd);
    const double weight = omega_s[static_cast<size_t>(ctx.S) * b];
    denom += weight * bin_mass;

    const double selected_mass = data.effect_sign == 1 ?
      funnel_interval_prob(lower, std::min(upper, q_signed), mean, sd) :
      funnel_interval_prob(std::max(lower, q_signed), upper, mean, sd);
    mass += weight * selected_mass;
  }

  if (!(denom > 0)) {
    return NA_REAL;
  }

  return clamp_probability_value(mass / denom);
}

void funnel_prepare_se_context(FunnelSeContext *ctx)
{
  const double eps_sd = std::sqrt(DBL_EPSILON);
  const SelNormKernelData &data = *ctx->data;
  const double se2 = ctx->se * ctx->se;
  const int n_boundaries = std::max(data.n_bins - 1, 0);

  ctx->location.assign(static_cast<size_t>(ctx->S), 0);
  ctx->total_sd.assign(static_cast<size_t>(ctx->S), 0);
  ctx->mean_z.assign(static_cast<size_t>(ctx->S), 0);
  ctx->sd_z.assign(static_cast<size_t>(ctx->S), 0);
  ctx->log_norm.assign(static_cast<size_t>(ctx->S), 0);
  ctx->normalizer.assign(static_cast<size_t>(ctx->S), 0);
  ctx->boundary_tail.assign(
    static_cast<size_t>(ctx->S) * static_cast<size_t>(n_boundaries), 0
  );
  ctx->zero_sd.assign(static_cast<size_t>(ctx->S), 0);
  ctx->weighted.assign(static_cast<size_t>(ctx->S), 0);
  ctx->fast_step.assign(static_cast<size_t>(ctx->S), 0);

  for (int s = 0; s < ctx->S; ++s) {
    const size_t si = static_cast<size_t>(s);
    const double location = ctx->mu[s] +
      ctx->direction * ctx->PET[s] * ctx->se +
      ctx->direction * ctx->PEESE[s] * se2;
    const double total_sd = std::sqrt(se2 + ctx->tau[s] * ctx->tau[s]);

    ctx->location[si] = location;
    ctx->total_sd[si] = total_sd;
    ctx->zero_sd[si]  = total_sd < eps_sd;
    ctx->weighted[si] = ctx->is_weightfunction[s] != 0 && !ctx->zero_sd[si];

    if (ctx->weighted[si] && ctx->se > 0) {
      if ((*ctx->row_active_phack)[si]) {
        ctx->log_norm[si] = std::numeric_limits<double>::quiet_NaN();
      } else {
        const int row_mode = (*ctx->row_mode)[si];
        const bool step_mode = row_mode == SELKERNEL_STEP ||
          row_mode == SELKERNEL_STEP_PHACK_POWER;
        const double mean_z = data.effect_sign * location / ctx->se;
        const double sd_z = total_sd / ctx->se;

        ctx->mean_z[si] = mean_z;
        ctx->sd_z[si] = sd_z;

        if (step_mode && data.trusted_step_partition &&
            data.telescope_probabilities) {
          const double *omega_s = ctx->omega + s;
          double normalizer = omega_s[static_cast<size_t>(ctx->S) *
            static_cast<size_t>(data.n_bins - 1)];
          double norm_abs = std::abs(normalizer);
          bool ok = std::isfinite(normalizer);

          for (int b = 0; ok && b < n_boundaries; ++b) {
            const double w_left = omega_s[static_cast<size_t>(ctx->S) *
              static_cast<size_t>(b)];
            const double w_right = omega_s[static_cast<size_t>(ctx->S) *
              static_cast<size_t>(b + 1)];
            if (w_left == w_right) {
              continue;
            }

            const double tail_boundary = pnorm(
              data.z_lower[b], mean_z, sd_z, false, false
            );
            if (!std::isfinite(tail_boundary)) {
              ok = false;
              break;
            }

            ctx->boundary_tail[
              si * static_cast<size_t>(n_boundaries) + static_cast<size_t>(b)
            ] = tail_boundary;

            const double contribution = (w_left - w_right) * tail_boundary;
            normalizer += contribution;
            norm_abs += std::abs(contribution);
          }

          if (ok && std::isfinite(normalizer) && normalizer > 0 &&
              normalizer > std::numeric_limits<double>::min() &&
              normalizer > 1e-8 * norm_abs) {
            ctx->normalizer[si] = normalizer;
            ctx->log_norm[si] = std::log(normalizer);
            ctx->fast_step[si] = 1;
            continue;
          }
        }

        ctx->log_norm[si] = cpp_selnorm_kernel_log_norm(
          location,
          total_sd,
          ctx->se,
          ctx->omega + s,
          (*ctx->row_alpha)[si],
          (*ctx->row_phack)[si],
          (*ctx->row_mode)[si],
          data,
          ctx->S,
          false
        );
      }
    }
  }
}

double funnel_step_selected_cdf_precomputed(double q, int s,
                                            const FunnelSeContext &ctx)
{
  const SelNormKernelData &data = *ctx.data;
  const size_t si = static_cast<size_t>(s);
  const int n_boundaries = std::max(data.n_bins - 1, 0);
  const double q_z = data.effect_sign * q / ctx.se;
  const double tail_q = pnorm(q_z, ctx.mean_z[si], ctx.sd_z[si], false, false);
  const double cdf_q = pnorm(q_z, ctx.mean_z[si], ctx.sd_z[si], true, false);

  if (!std::isfinite(tail_q) || !std::isfinite(cdf_q)) {
    return NA_REAL;
  }

  const double *omega_s = ctx.omega + s;
  const bool signed_lower_tail = data.effect_sign == 1;
  const double normalizer = ctx.normalizer[si];
  double requested = signed_lower_tail ?
    omega_s[static_cast<size_t>(ctx.S) * static_cast<size_t>(data.n_bins - 1)] *
      cdf_q :
    omega_s[static_cast<size_t>(ctx.S) * static_cast<size_t>(data.n_bins - 1)] *
      tail_q;
  double req_abs = std::abs(requested);

  for (int b = 0; b < n_boundaries; ++b) {
    const double w_left = omega_s[static_cast<size_t>(ctx.S) *
      static_cast<size_t>(b)];
    const double w_right = omega_s[static_cast<size_t>(ctx.S) *
      static_cast<size_t>(b + 1)];
    if (w_left == w_right) {
      continue;
    }

    const double diff = w_left - w_right;
    const double tail_boundary = ctx.boundary_tail[
      si * static_cast<size_t>(n_boundaries) + static_cast<size_t>(b)
    ];
    double req_contrib = 0;

    if (signed_lower_tail) {
      if (q_z > data.z_lower[b]) {
        double interval = tail_boundary - tail_q;
        const double scale = std::max(std::abs(tail_boundary), std::abs(tail_q));
        const double tol = 1e-14 * std::max(scale, 1.0);
        if (interval < 0 && interval > -tol) {
          interval = 0;
        } else if (interval < 0) {
          return NA_REAL;
        }
        req_contrib = diff * interval;
      }
    } else {
      req_contrib = diff * (q_z >= data.z_lower[b] ? tail_q : tail_boundary);
    }

    requested += req_contrib;
    req_abs += std::abs(req_contrib);
  }

  const double req_tol = 1e-10 *
    std::max(std::max(req_abs, normalizer), 1.0);
  if (requested < 0 && requested > -req_tol) {
    requested = 0;
  } else if (requested < 0) {
    return NA_REAL;
  }
  if (requested > normalizer && requested < normalizer + req_tol) {
    requested = normalizer;
  } else if (requested > normalizer) {
    return NA_REAL;
  }

  return clamp_probability_value(requested / normalizer);
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
      } else if (ctx.fast_step[si]) {
        cdf_s = funnel_step_selected_cdf_precomputed(q, s, ctx);
        if (!std::isfinite(cdf_s)) {
          cdf_s = cpp_selnorm_kernel_cdf_with_log_norm(
            q,
            ctx.location[si],
            ctx.total_sd[si],
            ctx.se,
            ctx.omega + s,
            (*ctx.row_alpha)[si],
            (*ctx.row_phack)[si],
            (*ctx.row_mode)[si],
            *ctx.data,
            ctx.log_norm[si],
            ctx.S,
            true,
            false
          );
        }
      } else {
        cdf_s = cpp_selnorm_kernel_cdf_with_log_norm(
          q,
          ctx.location[si],
          ctx.total_sd[si],
          ctx.se,
          ctx.omega + s,
          (*ctx.row_alpha)[si],
          (*ctx.row_phack)[si],
          (*ctx.row_mode)[si],
          *ctx.data,
          ctx.log_norm[si],
          ctx.S,
          true,
          false
        );
      }
    } else {
      cdf_s = pnorm(q, ctx.location[si], ctx.total_sd[si], true, false);
    }

    cdf_sum += clamp_probability_value(cdf_s);
  }

  return cdf_sum / static_cast<double>(ctx.S);
}

double funnel_grid_quantile_native(double p, double lower, double upper,
                                   const FunnelSeContext &ctx)
{
  const int n_grid = 1000;

  if (lower == upper) {
    return upper;
  }

  for (int i = 0; i < n_grid; ++i) {
    const double q = lower +
      (upper - lower) * static_cast<double>(i) / static_cast<double>(n_grid - 1);
    const double cdf = funnel_model_averaged_cdf_native(q, ctx);
    if (std::isfinite(cdf) && cdf >= p) {
      return q;
    }
  }

  return upper;
}

double funnel_model_averaged_quantile_native(double p,
                                             const FunnelSeContext &ctx,
                                             double previous,
                                             bool use_previous)
{
  const double eps_sd = std::sqrt(DBL_EPSILON);
  bool all_zero_sd = true;
  double lower = INFINITY;
  double upper = -INFINITY;
  double step = 0;
  double max_abs_location = 1;

  for (int s = 0; s < ctx.S; ++s) {
    const size_t si = static_cast<size_t>(s);
    all_zero_sd = all_zero_sd && ctx.zero_sd[si];
    const double spread = std::max(ctx.total_sd[si], eps_sd);
    lower = std::min(lower, ctx.location[si] - 10 * spread);
    upper = std::max(upper, ctx.location[si] + 10 * spread);
    step  = std::max(step, spread);

    if (std::isfinite(ctx.location[si])) {
      max_abs_location = std::max(max_abs_location, std::fabs(ctx.location[si]));
    }
  }

  if (all_zero_sd) {
    return funnel_quantile_type8(ctx.location, p);
  }

  if (!std::isfinite(lower) || !std::isfinite(upper)) {
    return NA_REAL;
  }
  if (lower >= upper) {
    lower -= 1;
    upper += 1;
  }
  if (!std::isfinite(step) || step <= 0) {
    step = max_abs_location;
  }

  const double global_lower = lower;
  const double global_upper = upper;
  const double global_step  = step;
  bool using_previous = use_previous && std::isfinite(previous);

  if (using_previous) {
    step = std::max(1e-3, 0.05 * global_step);
    lower = previous - step;
    upper = previous + step;
  }

  double lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
  double upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;

  if (!std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    if (!using_previous) {
      return funnel_grid_quantile_native(p, lower, upper, ctx);
    }
    using_previous = false;
    lower = global_lower;
    upper = global_upper;
    step = global_step;
    lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
    upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;
    if (!std::isfinite(lower_value) || !std::isfinite(upper_value)) {
      return funnel_grid_quantile_native(p, lower, upper, ctx);
    }
  }

  for (int i = 0; i < 25; ++i) {
    if (lower_value <= 0 && upper_value >= 0) {
      break;
    }
    if (lower_value > 0) {
      lower -= step;
      lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
    }
    if (upper_value < 0) {
      upper += step;
      upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;
    }
    step *= 2;
  }

  if (lower_value > 0 || upper_value < 0 ||
      !std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    if (using_previous) {
      lower = global_lower;
      upper = global_upper;
      step = global_step;
      lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
      upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;

      for (int i = 0; i < 25; ++i) {
        if (lower_value <= 0 && upper_value >= 0) {
          break;
        }
        if (lower_value > 0) {
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

  if (lower_value > 0 || upper_value < 0 ||
      !std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    return funnel_grid_quantile_native(p, lower, upper, ctx);
  }

  for (int i = 0; i < 100; ++i) {
    const double mid = lower + 0.5 * (upper - lower);
    const double mid_value = funnel_model_averaged_cdf_native(mid, ctx) - p;

    if (!std::isfinite(mid_value)) {
      return funnel_grid_quantile_native(p, lower, upper, ctx);
    }

    if (mid_value >= 0) {
      upper = mid;
      upper_value = mid_value;
    } else {
      lower = mid;
      lower_value = mid_value;
    }

    if (std::fabs(upper - lower) <= 1e-6 || mid_value == 0) {
      break;
    }
  }

  return lower + 0.5 * (upper - lower);
}

struct NormClusterArgs {
  int S;
  int K;
  int G;
  int n_gamma;
  bool selected;
  const double *yi;
  const double *sei;
  const double *mu;
  const double *tau_within;
  const double *tau_between;
  const int *cluster_index;
  const int *cluster_size;
  const double *weights;
  const double *gamma;
  const double *log_gamma;
  const double *omega;
  const double *alpha;
  const int *phack_kind;
  const int *kernel_mode;
  const int *obs_bin;
  int alpha_len;
  int phack_kind_len;
  int kernel_mode_len;
  SelNormKernelData selection_data;
};

double norm_cluster_logsumexp(const std::vector<double> &terms)
{
  double term_max = -INFINITY;

  for (std::vector<double>::const_iterator it = terms.begin(); it != terms.end(); ++it) {
    if (std::isnan(*it)) {
      return NA_REAL;
    }
    if (*it > term_max) {
      term_max = *it;
    }
  }

  if (std::isinf(term_max)) {
    return term_max;
  }

  double term_sum = 0;
  for (std::vector<double>::const_iterator it = terms.begin(); it != terms.end(); ++it) {
    term_sum += std::exp(*it - term_max);
  }

  return term_max + std::log(term_sum);
}

void validate_cluster_common_norm(SEXP yi, SEXP sei, SEXP mu_samples,
                                  SEXP tau_within, SEXP tau_between,
                                  SEXP cluster_index, SEXP cluster_size,
                                  SEXP weights, SEXP gamma_grid,
                                  SEXP log_gamma_weights,
                                  int *S, int *K, int *G, int *n_gamma)
{
  check_real(yi, "yi");
  check_real(sei, "sei");
  check_real(mu_samples, "mu_samples");
  check_real(tau_within, "tau_within");
  check_real(tau_between, "tau_between");
  check_integer(cluster_index, "cluster_index");
  check_integer(cluster_size, "cluster_size");
  check_real(gamma_grid, "gamma_grid");
  check_real(log_gamma_weights, "log_gamma_weights");

  if (weights != R_NilValue) {
    check_real(weights, "weights");
  }

  *S = matrix_nrow(mu_samples, "mu_samples");
  *K = matrix_ncol(mu_samples, "mu_samples");
  *G = Rf_length(cluster_size);
  *n_gamma = Rf_length(gamma_grid);

  validate_matrix_dims_named(tau_within, "tau_within", *S, *K, "mu_samples");
  validate_matrix_dims_named(tau_between, "tau_between", *S, *K, "mu_samples");

  if (Rf_length(yi) != *K || Rf_length(sei) != *K) {
    Rf_error("'yi' and 'sei' must have one value per observation.");
  }
  if (weights != R_NilValue && Rf_length(weights) != *K) {
    Rf_error("'weights' must have one value per observation.");
  }
  if (*n_gamma < 1 || Rf_length(log_gamma_weights) != *n_gamma) {
    Rf_error("'gamma_grid' and 'log_gamma_weights' must have matching positive lengths.");
  }

  const int *cluster_index_p = INTEGER(cluster_index);
  const int *cluster_size_p  = INTEGER(cluster_size);
  const int n_cluster_index  = Rf_length(cluster_index);
  int cluster_total = 0;

  for (int g = 0; g < *G; ++g) {
    if (cluster_size_p[g] <= 0) {
      Rf_error("'cluster_size' must contain positive cluster sizes.");
    }
    cluster_total += cluster_size_p[g];
  }
  if (cluster_total != n_cluster_index) {
    Rf_error("'cluster_size' must sum to the length of 'cluster_index'.");
  }
  for (int i = 0; i < n_cluster_index; ++i) {
    if (cluster_index_p[i] < 1 || cluster_index_p[i] > *K) {
      Rf_error("'cluster_index' values must be between 1 and the number of observations.");
    }
  }
}

double norm_cluster_cell_loglik(const NormClusterArgs &args, int s, int k,
                                double mu_node, double sd)
{
  const double weight = args.weights == 0 ? 1.0 : args.weights[k];

  if (args.selected) {
    return weight * cpp_selnorm_kernel_lpdf(
      args.yi[k],
      mu_node,
      sd,
      mu_node,
      sd,
      args.sei[k],
      1.0,
      args.omega + s,
      args.obs_bin[k],
      scalar_or_row_real(args.alpha, args.alpha_len, s, "alpha"),
      scalar_or_row_int(args.phack_kind, args.phack_kind_len, s, "phack_kind"),
      scalar_or_row_int(args.kernel_mode, args.kernel_mode_len, s, "kernel_mode"),
      args.selection_data,
      args.S,
      false
    );
  }

  return weight * Rf_dnorm4(args.yi[k], mu_node, sd, 1);
}

SEXP norm_cluster_accumulate(const NormClusterArgs &args)
{
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, args.S, args.G));
  double *out_p = REAL(out);

  std::vector<double> gamma_terms(static_cast<size_t>(args.n_gamma));
  std::vector<double> cluster_terms(
    static_cast<size_t>(args.S) * static_cast<size_t>(args.n_gamma)
  );

  int cluster_offset = 0;
  for (int g = 0; g < args.G; ++g) {
    const int cluster_n = args.cluster_size[g];
    std::fill(cluster_terms.begin(), cluster_terms.end(), 0.0);

    for (int m = 0; m < cluster_n; ++m) {
      const int k = args.cluster_index[cluster_offset + m] - 1;

      for (int j = 0; j < args.n_gamma; ++j) {
        const double gamma_j = args.gamma[j];

        for (int s = 0; s < args.S; ++s) {
          const int cell = s + args.S * k;
          const double mu_node = args.mu[cell] +
            gamma_j * args.tau_between[cell];
          const double sd = std::sqrt(
            args.tau_within[cell] * args.tau_within[cell] +
              args.sei[k] * args.sei[k]
          );

          cluster_terms[s + args.S * j] += norm_cluster_cell_loglik(
            args, s, k, mu_node, sd
          );
        }
      }
    }

    for (int s = 0; s < args.S; ++s) {
      for (int j = 0; j < args.n_gamma; ++j) {
        gamma_terms[j] = cluster_terms[s + args.S * j] + args.log_gamma[j];
      }
      out_p[s + args.S * g] = norm_cluster_logsumexp(gamma_terms);
    }

    cluster_offset += cluster_n;
  }

  UNPROTECT(1);
  return out;
}

SEXP norm_cluster_accumulate_row_sum(const NormClusterArgs &args)
{
  SEXP out = PROTECT(Rf_allocVector(REALSXP, args.S));
  double *out_p = REAL(out);
  std::fill(out_p, out_p + args.S, 0.0);

  std::vector<double> gamma_terms(static_cast<size_t>(args.n_gamma));
  std::vector<double> cluster_terms(
    static_cast<size_t>(args.S) * static_cast<size_t>(args.n_gamma)
  );

  int cluster_offset = 0;
  for (int g = 0; g < args.G; ++g) {
    const int cluster_n = args.cluster_size[g];
    std::fill(cluster_terms.begin(), cluster_terms.end(), 0.0);

    for (int m = 0; m < cluster_n; ++m) {
      const int k = args.cluster_index[cluster_offset + m] - 1;

      for (int j = 0; j < args.n_gamma; ++j) {
        const double gamma_j = args.gamma[j];

        for (int s = 0; s < args.S; ++s) {
          const int cell = s + args.S * k;
          const double mu_node = args.mu[cell] +
            gamma_j * args.tau_between[cell];
          const double sd = std::sqrt(
            args.tau_within[cell] * args.tau_within[cell] +
              args.sei[k] * args.sei[k]
          );

          cluster_terms[s + args.S * j] += norm_cluster_cell_loglik(
            args, s, k, mu_node, sd
          );
        }
      }
    }

    for (int s = 0; s < args.S; ++s) {
      for (int j = 0; j < args.n_gamma; ++j) {
        gamma_terms[j] = cluster_terms[s + args.S * j] + args.log_gamma[j];
      }
      out_p[s] += norm_cluster_logsumexp(gamma_terms);
    }

    cluster_offset += cluster_n;
  }

  UNPROTECT(1);
  return out;
}

}
