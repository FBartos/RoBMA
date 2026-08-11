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

double clamp_probability_value(double x)
{
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

  return clamp_probability_value(requested_mass / total_mass);
}

void funnel_prepare_se_context(FunnelSeContext *ctx)
{
  const double se2 = ctx->se * ctx->se;

  ctx->location.assign(static_cast<size_t>(ctx->S), 0);
  ctx->total_sd.assign(static_cast<size_t>(ctx->S), 0);
  ctx->zero_sd.assign(static_cast<size_t>(ctx->S), 0);
  ctx->weighted.assign(static_cast<size_t>(ctx->S), 0);

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

    cdf_sum += clamp_probability_value(cdf_s);
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

const int SELNORM_CLUSTER_RULE_7_OFFSET  = 0;
const int SELNORM_CLUSTER_RULE_15_OFFSET = 7;
const int SELNORM_CLUSTER_RULE_31_OFFSET = 22;
const int SELNORM_CLUSTER_RULE_TOTAL     = 53;
// The displaced-mode certification case has a 3.64e-4 successive-rule change
// at 31 nodes and a 5.61e-6 log-integral error against independent integration.
const double SELNORM_CLUSTER_QUADRATURE_TOLERANCE = 5e-4;

void validate_selnorm_cluster_quadrature_rules(int n_gamma)
{
  if (n_gamma != SELNORM_CLUSTER_RULE_TOTAL) {
    Rf_error(
      "Selected-normal cluster quadrature requires concatenated 7, 15, and 31 node rules."
    );
  }
}

struct NormClusterModeContext {
  const NormClusterArgs *args;
  int s;
  int cluster_offset;
  int cluster_n;
  const double *basis;
  double delta;
};

struct NormClusterModeBracket {
  double lower;
  double upper;
  double initial;
};

struct NormClusterIntegral {
  double value;
  double relative_change;
  int order;
};

double norm_cluster_logsumexp(const std::vector<double> &terms, int n = -1)
{
  const int n_terms = n < 0 ? static_cast<int>(terms.size()) : n;
  double term_max = -INFINITY;

  for (int j = 0; j < n_terms; ++j) {
    if (std::isnan(terms[j])) {
      return NA_REAL;
    }
    if (terms[j] > term_max) {
      term_max = terms[j];
    }
  }

  if (std::isinf(term_max)) {
    return term_max;
  }

  double term_sum = 0;
  for (int j = 0; j < n_terms; ++j) {
    term_sum += std::exp(terms[j] - term_max);
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

double norm_cluster_loglik_at_gamma(const NormClusterArgs &args, int s,
                                    int cluster_offset, int cluster_n,
                                    double gamma, const double *basis,
                                    double delta)
{
  double out = 0.0;

  for (int m = 0; m < cluster_n; ++m) {
    const int k = args.cluster_index[cluster_offset + m] - 1;
    const int cell = s + args.S * k;
    const double mu_node = args.mu[cell] +
      (basis == 0 ? 0.0 : basis[cell] * delta) +
      gamma * args.tau_between[cell];
    const double sd = std::sqrt(
      args.tau_within[cell] * args.tau_within[cell] +
        args.sei[k] * args.sei[k]
    );

    out += norm_cluster_cell_loglik(args, s, k, mu_node, sd);
  }

  return out;
}

double norm_cluster_log_posterior_gamma(const NormClusterArgs &args, int s,
                                        int cluster_offset, int cluster_n,
                                        double gamma, const double *basis,
                                        double delta)
{
  const double log_lik = norm_cluster_loglik_at_gamma(
    args, s, cluster_offset, cluster_n, gamma, basis, delta
  );
  if (!std::isfinite(log_lik)) {
    return log_lik;
  }

  return log_lik - 0.5 * gamma * gamma;
}

double norm_cluster_start_mode(const NormClusterArgs &args, int s,
                               int cluster_offset, int cluster_n,
                               const double *basis, double delta,
                               double *curvature)
{
  double precision = 1.0;
  double numerator = 0.0;

  for (int m = 0; m < cluster_n; ++m) {
    const int k = args.cluster_index[cluster_offset + m] - 1;
    const int cell = s + args.S * k;
    const double sd2 = args.tau_within[cell] * args.tau_within[cell] +
      args.sei[k] * args.sei[k];
    const double tau_between = args.tau_between[cell];
    const double mu_base = args.mu[cell] +
      (basis == 0 ? 0.0 : basis[cell] * delta);
    const double weight = args.weights == 0 ? 1.0 : args.weights[k];

    if (!std::isfinite(sd2) || sd2 <= 0 ||
        !std::isfinite(tau_between) || !std::isfinite(mu_base) ||
        !std::isfinite(weight) || weight < 0) {
      continue;
    }

    precision += weight * tau_between * tau_between / sd2;
    numerator += weight * tau_between * (args.yi[k] - mu_base) / sd2;
  }

  if (!std::isfinite(precision) || precision <= 0) {
    *curvature = 1.0;
    return 0.0;
  }

  *curvature = precision;
  const double mode = numerator / precision;
  return std::isfinite(mode) ? mode : 0.0;
}

double norm_cluster_mode_value(const NormClusterModeContext &context,
                               double gamma)
{
  return norm_cluster_log_posterior_gamma(
    *context.args, context.s, context.cluster_offset, context.cluster_n,
    gamma, context.basis, context.delta
  );
}

bool norm_cluster_mode_bracket(const NormClusterModeContext &context,
                               double start, double step,
                               NormClusterModeBracket *bracket)
{
  double center       = start;
  double center_value = norm_cluster_mode_value(context, center);
  double left         = center - step;
  double right        = center + step;
  const double left_value  = norm_cluster_mode_value(context, left);
  const double right_value = norm_cluster_mode_value(context, right);

  if (!std::isfinite(center_value) || std::isnan(left_value) ||
      std::isnan(right_value) || left_value == INFINITY ||
      right_value == INFINITY) {
    return false;
  }
  if (center_value >= left_value && center_value >= right_value) {
    bracket->lower   = left;
    bracket->upper   = right;
    bracket->initial = center;
    return true;
  }
  if (left_value > center_value && right_value > center_value) {
    return false;
  }

  const bool move_left = left_value > right_value;
  if (move_left) {
    right        = center;
    center       = left;
    center_value = left_value;
  } else {
    left         = center;
    center       = right;
    center_value = right_value;
  }

  for (int expansion = 0; expansion < 20; ++expansion) {
    step *= 2.0;
    const double candidate = move_left ? center - step : center + step;
    const double candidate_value = norm_cluster_mode_value(context, candidate);

    if (std::isnan(candidate_value) || candidate_value == INFINITY) {
      return false;
    }
    if (candidate_value <= center_value) {
      bracket->lower   = move_left ? candidate : left;
      bracket->upper   = move_left ? right : candidate;
      bracket->initial = center;
      return std::isfinite(bracket->lower) &&
        std::isfinite(bracket->upper) &&
        bracket->lower < bracket->initial &&
        bracket->initial < bracket->upper;
    }

    if (move_left) {
      right = center;
    } else {
      left = center;
    }
    center       = candidate;
    center_value = candidate_value;
  }

  return false;
}

bool norm_cluster_brent_mode(const NormClusterModeContext &context,
                             const NormClusterModeBracket &bracket,
                             double *mode, double *mode_value)
{
  const double golden = 0.3819660112501051;
  const double absolute_tolerance = 1e-8;
  const double sqrt_epsilon = std::sqrt(DBL_EPSILON);
  double lower = bracket.lower;
  double upper = bracket.upper;
  double x = bracket.initial;
  double w = x;
  double v = x;
  double fx = -norm_cluster_mode_value(context, x);
  double fw = fx;
  double fv = fx;
  double d = 0.0;
  double e = 0.0;

  if (!std::isfinite(fx)) {
    return false;
  }

  for (int iteration = 0; iteration < 100; ++iteration) {
    const double midpoint = 0.5 * (lower + upper);
    const double tol1 = sqrt_epsilon * std::fabs(x) + absolute_tolerance;
    const double tol2 = 2.0 * tol1;

    if (std::fabs(x - midpoint) <= tol2 - 0.5 * (upper - lower)) {
      *mode       = x;
      *mode_value = -fx;
      return std::isfinite(*mode_value);
    }

    bool use_parabola = false;
    if (std::fabs(e) > tol1 && std::isfinite(fw) && std::isfinite(fv)) {
      const double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) {
        p = -p;
      } else {
        q = -q;
      }
      const double previous_e = e;
      e = d;
      if (q > 0.0 && std::fabs(p) < std::fabs(0.5 * q * previous_e) &&
          p > q * (lower - x) && p < q * (upper - x)) {
        d = p / q;
        const double proposed = x + d;
        if (proposed - lower < tol2 || upper - proposed < tol2) {
          d = midpoint >= x ? tol1 : -tol1;
        }
        use_parabola = true;
      }
    }

    if (!use_parabola) {
      e = x >= midpoint ? lower - x : upper - x;
      d = golden * e;
    }

    const double u = std::fabs(d) >= tol1 ? x + d :
      x + (d > 0.0 ? tol1 : -tol1);
    const double log_fu = norm_cluster_mode_value(context, u);
    if (std::isnan(log_fu) || log_fu == INFINITY) {
      return false;
    }
    const double fu = -log_fu;

    if (fu <= fx) {
      if (u >= x) {
        lower = x;
      } else {
        upper = x;
      }
      v  = w;
      fv = fw;
      w  = x;
      fw = fx;
      x  = u;
      fx = fu;
    } else {
      if (u < x) {
        lower = u;
      } else {
        upper = u;
      }
      if (fu <= fw || w == x) {
        v  = w;
        fv = fw;
        w  = u;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v  = u;
        fv = fu;
      }
    }
  }

  return false;
}

double norm_cluster_rule_log_integral(const NormClusterArgs &args, int s,
                                      int cluster_offset, int cluster_n,
                                      const double *basis, double delta,
                                      double mode, double curvature,
                                      int rule_offset, int rule_order,
                                      std::vector<double> &terms)
{
  const double inv_sqrt_curvature = 1.0 / std::sqrt(curvature);
  const double log_jacobian = -0.5 * std::log(curvature);

  for (int j = 0; j < rule_order; ++j) {
    const int rule_index = rule_offset + j;
    const double z = args.gamma[rule_index];
    const double gamma = mode + z * inv_sqrt_curvature;
    const double log_lik = norm_cluster_loglik_at_gamma(
      args, s, cluster_offset, cluster_n, gamma, basis, delta
    );
    if (std::isnan(log_lik) || log_lik == INFINITY) {
      return NA_REAL;
    }
    terms[j] = args.log_gamma[rule_index] + log_lik -
      0.5 * gamma * gamma + 0.5 * z * z + log_jacobian;
  }

  return norm_cluster_logsumexp(terms, rule_order);
}

double norm_cluster_relative_change(double coarse, double fine)
{
  const double difference = coarse - fine;
  if (!std::isfinite(difference) || std::fabs(difference) > 50.0) {
    return INFINITY;
  }

  return std::fabs(std::expm1(difference));
}

NormClusterIntegral norm_cluster_adaptive_log_integral(
  const NormClusterArgs &args, int s, int cluster_id, int cluster_offset,
  int cluster_n, const double *basis, double delta,
  std::vector<double> &terms)
{
  double start_curvature = 1.0;
  const double start = norm_cluster_start_mode(
    args, s, cluster_offset, cluster_n, basis, delta, &start_curvature
  );
  if (!std::isfinite(start) || !std::isfinite(start_curvature) ||
      start_curvature <= 0.0) {
    Rf_error(
      "Selected-normal cluster quadrature could not initialize the mode for posterior row %d and cluster %d.",
      s + 1, cluster_id + 1
    );
  }

  NormClusterModeContext context;
  context.args           = &args;
  context.s              = s;
  context.cluster_offset = cluster_offset;
  context.cluster_n      = cluster_n;
  context.basis          = basis;
  context.delta          = delta;

  NormClusterModeBracket bracket;
  const double start_scale = 1.0 / std::sqrt(start_curvature);
  if (!std::isfinite(start_scale) || start_scale <= 0.0 ||
      !norm_cluster_mode_bracket(context, start, start_scale, &bracket)) {
    Rf_error(
      "Selected-normal cluster quadrature could not bracket a finite mode for posterior row %d and cluster %d.",
      s + 1, cluster_id + 1
    );
  }

  double mode = start;
  double mode_value = NA_REAL;
  if (!norm_cluster_brent_mode(context, bracket, &mode, &mode_value)) {
    Rf_error(
      "Selected-normal cluster quadrature mode refinement failed for posterior row %d and cluster %d.",
      s + 1, cluster_id + 1
    );
  }
  const double curvature_step = 1e-3 * std::max(1.0, std::fabs(mode));
  const double value_plus = norm_cluster_mode_value(context, mode + curvature_step);
  const double value_minus = norm_cluster_mode_value(context, mode - curvature_step);
  const double curvature = -(value_plus - 2.0 * mode_value + value_minus) /
    (curvature_step * curvature_step);
  if (!std::isfinite(curvature) || curvature <= 0.0) {
    Rf_error(
      "Selected-normal cluster quadrature found non-positive mode curvature for posterior row %d and cluster %d.",
      s + 1, cluster_id + 1
    );
  }

  const double integral_7 = norm_cluster_rule_log_integral(
    args, s, cluster_offset, cluster_n, basis, delta, mode, curvature,
    SELNORM_CLUSTER_RULE_7_OFFSET, 7, terms
  );
  const double integral_15 = norm_cluster_rule_log_integral(
    args, s, cluster_offset, cluster_n, basis, delta, mode, curvature,
    SELNORM_CLUSTER_RULE_15_OFFSET, 15, terms
  );
  if (!std::isfinite(integral_7) || !std::isfinite(integral_15)) {
    Rf_error(
      "Selected-normal cluster quadrature produced a non-finite integral for posterior row %d and cluster %d.",
      s + 1, cluster_id + 1
    );
  }

  NormClusterIntegral result;
  result.relative_change = norm_cluster_relative_change(integral_7, integral_15);
  if (result.relative_change <= SELNORM_CLUSTER_QUADRATURE_TOLERANCE) {
    result.value = integral_15;
    result.order = 15;
    return result;
  }

  const double integral_31 = norm_cluster_rule_log_integral(
    args, s, cluster_offset, cluster_n, basis, delta, mode, curvature,
    SELNORM_CLUSTER_RULE_31_OFFSET, 31, terms
  );
  if (!std::isfinite(integral_31)) {
    Rf_error(
      "Selected-normal cluster quadrature produced a non-finite integral for posterior row %d and cluster %d.",
      s + 1, cluster_id + 1
    );
  }

  result.relative_change = norm_cluster_relative_change(integral_15, integral_31);
  if (result.relative_change > SELNORM_CLUSTER_QUADRATURE_TOLERANCE) {
    Rf_error(
      "Selected-normal cluster quadrature failed to converge for posterior row %d and cluster %d (relative change %.6g after 31 nodes).",
      s + 1, cluster_id + 1, result.relative_change
    );
  }

  result.value = integral_31;
  result.order = 31;
  return result;
}

SEXP norm_cluster_accumulate(const NormClusterArgs &args)
{
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, args.S, args.G));
  double *out_p = REAL(out);

  std::vector<double> gamma_terms(
    static_cast<size_t>(args.selected ? 31 : args.n_gamma)
  );
  if (args.selected) {
    int cluster_offset = 0;
    for (int g = 0; g < args.G; ++g) {
      const int cluster_n = args.cluster_size[g];
      for (int s = 0; s < args.S; ++s) {
        const NormClusterIntegral integral = norm_cluster_adaptive_log_integral(
          args, s, g, cluster_offset, cluster_n, 0, 0.0, gamma_terms
        );
        out_p[s + args.S * g] = integral.value;
      }
      cluster_offset += cluster_n;
    }

    UNPROTECT(1);
    return out;
  }

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

  std::vector<double> gamma_terms(
    static_cast<size_t>(args.selected ? 31 : args.n_gamma)
  );
  if (args.selected) {
    int cluster_offset = 0;
    for (int g = 0; g < args.G; ++g) {
      const int cluster_n = args.cluster_size[g];
      for (int s = 0; s < args.S; ++s) {
        const NormClusterIntegral integral = norm_cluster_adaptive_log_integral(
          args, s, g, cluster_offset, cluster_n, 0, 0.0, gamma_terms
        );
        out_p[s] += integral.value;
      }
      cluster_offset += cluster_n;
    }

    UNPROTECT(1);
    return out;
  }

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
