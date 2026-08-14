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
