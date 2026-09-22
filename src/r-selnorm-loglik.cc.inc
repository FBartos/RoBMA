namespace {

// One candidate cell of the normal predictor grid, written exactly as
// RoBMA_norm_loglik_row_sum() writes it so that the two routes agree value for
// value: the same total standard deviation and the same weighted dnorm.
inline double norm_predictor_cell(double yi, double mu, double tau, double sei,
                                  double weight)
{
  const double sd = std::sqrt(tau * tau + sei * sei);
  return weight * Rf_dnorm4(yi, mu, sd, 1);
}

void check_predictor_matrix(SEXP x, const char *name, int S, int K)
{
  check_real(x, name);
  validate_matrix_dims(x, name, S, K);
}

}

// Weighted normal log-likelihood sums over a grid of candidate values for one
// posterior batch, without materializing the (S * G) x K candidate matrices
// .iwmde_predictor_candidates() builds. Column s of the result holds the G
// candidates of posterior row s, in the column-major order the R grid uses.
//
// The candidate location is mu[s, k] + mu_basis[s, k] * (values[g] -
// current[s]); its scale is values[g] when the grid replaces tau itself,
// otherwise tau[s, k], and either is carried through exp(log(tau) +
// log_tau_basis[s, k] * delta) when the batch has a log-tau basis. A log-tau
// basis of a logged formula intercept is additive in the logarithm of the
// focal coordinate instead, so 'log_tau_delta_log' replaces that delta with
// log(values[g]) - log(current[s]) and marks a non-positive pair invalid.
// 'valid' reproduces the candidate builder's validity mask, including the row
// sums it tests, so the caller's -Inf pattern is unchanged.
extern "C" SEXP RoBMA_norm_predictor_grid_loglik(
  SEXP yi, SEXP sei, SEXP weights, SEXP mu, SEXP mu_basis, SEXP tau,
  SEXP log_tau_basis, SEXP current, SEXP values, SEXP scale_update_tau,
  SEXP log_tau_delta_log)
{
  check_real(yi, "yi");
  check_real(sei, "sei");
  check_real(current, "current");
  check_real(values, "values");
  check_logical(scale_update_tau, "scale_update_tau");
  check_logical(log_tau_delta_log, "log_tau_delta_log");

  const int S = matrix_nrow(mu, "mu");
  const int K = matrix_ncol(mu, "mu");
  const int G = Rf_length(values);
  const bool scale_tau = LOGICAL(scale_update_tau)[0] == TRUE;
  const bool log_delta = LOGICAL(log_tau_delta_log)[0] == TRUE;
  if (log_delta && log_tau_basis == R_NilValue) {
    Rf_error("'log_tau_delta_log' needs a log-tau basis.");
  }
  if (log_delta && mu_basis != R_NilValue) {
    Rf_error("'log_tau_delta_log' cannot accompany a location basis.");
  }

  check_predictor_matrix(mu, "mu", S, K);
  if (mu_basis != R_NilValue) check_predictor_matrix(mu_basis, "mu_basis", S, K);
  if (log_tau_basis != R_NilValue) {
    check_predictor_matrix(log_tau_basis, "log_tau_basis", S, K);
  }
  if (scale_tau) {
    if (tau != R_NilValue) {
      Rf_error("'tau' must be absent when the grid replaces the scale itself.");
    }
  } else {
    check_predictor_matrix(tau, "tau", S, K);
  }
  if (Rf_length(yi) != K || Rf_length(sei) != K) {
    Rf_error("'yi' and 'sei' must have one value per observation.");
  }
  if (weights != R_NilValue) {
    check_real(weights, "weights");
    if (Rf_length(weights) != K) {
      Rf_error("'weights' must have one value per observation.");
    }
  }
  if (Rf_length(current) != S) {
    Rf_error("'current' must have one value per posterior sample.");
  }

  SEXP out_log_lik = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(G) * S));
  SEXP out_valid = PROTECT(Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(G) * S));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("log_lik"));
  SET_STRING_ELT(names, 1, Rf_mkChar("valid"));
  SET_VECTOR_ELT(out, 0, out_log_lik);
  SET_VECTOR_ELT(out, 1, out_valid);
  Rf_setAttrib(out, R_NamesSymbol, names);

  const double *yi_p        = REAL(yi);
  const double *sei_p       = REAL(sei);
  const double *weights_p   = weights == R_NilValue ? 0 : REAL(weights);
  const double *mu_p        = REAL(mu);
  const double *mu_basis_p  = mu_basis == R_NilValue ? 0 : REAL(mu_basis);
  const double *tau_p       = scale_tau ? 0 : REAL(tau);
  const double *log_tau_p   = log_tau_basis == R_NilValue ? 0 : REAL(log_tau_basis);
  const double *current_p   = REAL(current);
  const double *values_p    = REAL(values);
  double *log_lik_p         = REAL(out_log_lik);
  int *valid_p              = LOGICAL(out_valid);

  // The candidate matrices repeated one posterior row per grid value, so the
  // parts of the arithmetic that depend on only one of the two indices are
  // prepared once here. Every prepared value is the one the repeated row
  // carried: the logarithm of the same scale, and the row sums that
  // .iwmde_predictor_candidates() tested with rowSums(), accumulated over the
  // same observations in the same order and in the same long double.
  const bool has_log_tau = log_tau_p != 0;
  std::vector<double> log_base;
  if (has_log_tau) {
    if (scale_tau) {
      log_base.resize(static_cast<std::size_t>(G));
      for (int g = 0; g < G; ++g) log_base[g] = std::log(values_p[g]);
    } else {
      log_base.resize(static_cast<std::size_t>(S) * K);
      for (std::size_t i = 0; i < log_base.size(); ++i) log_base[i] = std::log(tau_p[i]);
    }
  }
  std::vector<double> mu_row;
  if (mu_basis_p == 0) {
    mu_row.resize(static_cast<std::size_t>(S));
    for (int s = 0; s < S; ++s) {
      long double sum = 0.0L;
      for (int k = 0; k < K; ++k) sum += static_cast<long double>(mu_p[s + S * k]);
      mu_row[static_cast<std::size_t>(s)] = static_cast<double>(sum);
    }
  }
  std::vector<double> tau_row;
  std::vector<char> tau_row_ok;
  if (!has_log_tau) {
    const int fixed = scale_tau ? G : S;
    tau_row.resize(static_cast<std::size_t>(fixed));
    tau_row_ok.resize(static_cast<std::size_t>(fixed));
    for (int i = 0; i < fixed; ++i) {
      long double sum = 0.0L;
      bool ok = true;
      for (int k = 0; k < K; ++k) {
        const double scale = scale_tau ? values_p[i] : tau_p[i + S * k];
        if (!std::isfinite(scale) || scale < 0.0) ok = false;
        sum += static_cast<long double>(scale);
      }
      tau_row[static_cast<std::size_t>(i)] = static_cast<double>(sum);
      tau_row_ok[static_cast<std::size_t>(i)] = ok ? 1 : 0;
    }
  }

  // Posterior rows are independent: row s reads call-common inputs and its own
  // slice, and writes only its own G output slots.
  const auto compute_row = [&](int s) {
    const double current_value = current_p[s];
    const std::size_t row = static_cast<std::size_t>(s);
    const bool fixed_mu = mu_basis_p == 0;
    for (int g = 0; g < G; ++g) {
      const R_xlen_t cell = static_cast<R_xlen_t>(s) * G + g;
      const double value = values_p[g];
      const double delta = value - current_value;
      bool ok = std::isfinite(value) && std::isfinite(delta);
      if (scale_tau && !(value >= 0.0)) ok = false;
      // The log-tau basis of a logged formula intercept multiplies the
      // difference of the logarithms; only a positive pair has one.
      double tau_delta = delta;
      if (log_delta) {
        if (value > 0.0 && current_value > 0.0) {
          tau_delta = std::log(value) - std::log(current_value);
        } else {
          tau_delta = 0.0;
          ok = false;
        }
      }
      long double mu_sum = 0.0L;
      long double tau_sum = 0.0L;
      long double log_tau_sum = 0.0L;
      bool tau_ok = true;
      double row_sum = 0.0;
      for (int k = 0; k < K; ++k) {
        const int slot = s + S * k;
        const double location = fixed_mu ? mu_p[slot] :
          mu_p[slot] + mu_basis_p[slot] * delta;
        double scale;
        if (has_log_tau) {
          const double log_scale = (scale_tau ? log_base[static_cast<std::size_t>(g)] :
            log_base[static_cast<std::size_t>(slot)]) + log_tau_p[slot] * tau_delta;
          log_tau_sum += static_cast<long double>(log_scale);
          scale = std::exp(log_scale);
          if (!std::isfinite(scale) || scale < 0.0) tau_ok = false;
          tau_sum += static_cast<long double>(scale);
        } else {
          scale = scale_tau ? value : tau_p[slot];
        }
        if (!fixed_mu) mu_sum += static_cast<long double>(location);
        row_sum += norm_predictor_cell(yi_p[k], location, scale, sei_p[k],
          weights_p == 0 ? 1.0 : weights_p[k]);
      }
      if (has_log_tau) {
        if (!std::isfinite(static_cast<double>(log_tau_sum))) ok = false;
        if (!tau_ok) ok = false;
        if (!std::isfinite(static_cast<double>(tau_sum))) ok = false;
      } else {
        const std::size_t fixed_index = scale_tau ?
          static_cast<std::size_t>(g) : row;
        if (tau_row_ok[fixed_index] == 0) ok = false;
        if (!std::isfinite(tau_row[fixed_index])) ok = false;
      }
      const double location_sum = fixed_mu ? mu_row[row] :
        static_cast<double>(mu_sum);
      if (!std::isfinite(location_sum)) ok = false;
      log_lik_p[cell] = row_sum;
      valid_p[cell] = ok ? TRUE : FALSE;
    }
  };

  // One density evaluation per grid value and observation.
  const RobmaRowSchedule schedule = robma_row_schedule(
    S, static_cast<double>(G) * K);
  if (schedule.threads == 1) {
    for (int s = 0; s < S; ++s) compute_row(s);
  } else {
    robma_parallel::for_rows(S, schedule, compute_row);
  }

  UNPROTECT(4);
  return out;
}

extern "C" SEXP RoBMA_norm_loglik_row_sum(
  SEXP yi, SEXP mu_samples, SEXP tau_within, SEXP sei, SEXP weights)
{
  check_real(yi, "yi");
  check_real(mu_samples, "mu_samples");
  check_real(tau_within, "tau_within");
  check_real(sei, "sei");
  if (weights != R_NilValue) {
    check_real(weights, "weights");
  }

  const int S = matrix_nrow(mu_samples, "mu_samples");
  const int K = matrix_ncol(mu_samples, "mu_samples");
  validate_matrix_dims_named(tau_within, "tau_within", S, K, "mu_samples");

  if (Rf_length(yi) != K || Rf_length(sei) != K) {
    Rf_error("'yi' and 'sei' must have one value per observation.");
  }
  if (weights != R_NilValue && Rf_length(weights) != K) {
    Rf_error("'weights' must have one value per observation.");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, S));

  const double *yi_p     = REAL(yi);
  const double *mu_p     = REAL(mu_samples);
  const double *tau_p    = REAL(tau_within);
  const double *sei_p    = REAL(sei);
  const double *weights_p = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p          = REAL(out);

  for (int s = 0; s < S; ++s) {
    double row_sum = 0.0;
    for (int k = 0; k < K; ++k) {
      const int cell = s + S * k;
      const double sd = std::sqrt(tau_p[cell] * tau_p[cell] + sei_p[k] * sei_p[k]);
      const double weight = weights_p == 0 ? 1.0 : weights_p[k];
      row_sum += weight * Rf_dnorm4(yi_p[k], mu_p[cell], sd, 1);
    }
    out_p[s] = row_sum;
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_norm_cluster_analytic_loglik(
  SEXP yi, SEXP vi, SEXP mu_samples, SEXP tau_within, SEXP tau_between,
  SEXP cluster_index, SEXP cluster_size)
{
  NormClusterAnalyticArgs args;
  validate_cluster_common_data(
    yi, vi, "vi", mu_samples, tau_within, tau_between, cluster_index,
    cluster_size, R_NilValue, &args.S, &args.K, &args.G
  );
  args.yi            = REAL(yi);
  args.vi            = REAL(vi);
  args.mu            = REAL(mu_samples);
  args.tau_within    = REAL(tau_within);
  args.tau_between   = REAL(tau_between);
  args.tau_total     = 0;
  args.cluster_index = INTEGER(cluster_index);
  args.cluster_size  = INTEGER(cluster_size);

  return norm_cluster_analytic_accumulate(args, false);
}

extern "C" SEXP RoBMA_norm_cluster_analytic_loglik_row_sum(
  SEXP yi, SEXP vi, SEXP mu_samples, SEXP tau_within, SEXP tau_between,
  SEXP cluster_index, SEXP cluster_size)
{
  NormClusterAnalyticArgs args;
  validate_cluster_common_data(
    yi, vi, "vi", mu_samples, tau_within, tau_between, cluster_index,
    cluster_size, R_NilValue, &args.S, &args.K, &args.G
  );
  args.yi            = REAL(yi);
  args.vi            = REAL(vi);
  args.mu            = REAL(mu_samples);
  args.tau_within    = REAL(tau_within);
  args.tau_between   = REAL(tau_between);
  args.tau_total     = 0;
  args.cluster_index = INTEGER(cluster_index);
  args.cluster_size  = INTEGER(cluster_size);

  return norm_cluster_analytic_accumulate(args, true);
}

extern "C" SEXP RoBMA_norm_cluster_analytic_rho_grid_loglik(
  SEXP yi, SEXP vi, SEXP mu_samples, SEXP tau_total, SEXP rho,
  SEXP cluster_index, SEXP cluster_size)
{
  NormClusterAnalyticArgs args;
  validate_cluster_common_data(
    yi, vi, "vi", mu_samples, tau_total, tau_total, cluster_index,
    cluster_size, R_NilValue, &args.S, &args.K, &args.G
  );
  check_real(rho, "rho");
  const int n_rho = Rf_length(rho);
  if (n_rho < 1) {
    Rf_error("'rho' must contain at least one value.");
  }
  const double *rho_p = REAL(rho);
  for (int r = 0; r < n_rho; ++r) {
    if (!std::isfinite(rho_p[r]) || rho_p[r] < 0.0 || rho_p[r] > 1.0) {
      Rf_error("'rho' must contain only finite values within [0, 1].");
    }
  }
  args.yi            = REAL(yi);
  args.vi            = REAL(vi);
  args.mu            = REAL(mu_samples);
  args.tau_within    = 0;
  args.tau_between   = 0;
  args.tau_total     = REAL(tau_total);
  args.cluster_index = INTEGER(cluster_index);
  args.cluster_size  = INTEGER(cluster_size);

  return norm_cluster_analytic_rho_grid_accumulate(args, rho_p, n_rho);
}

extern "C" SEXP RoBMA_norm_cluster_loglik(
  SEXP yi, SEXP sei, SEXP mu_samples, SEXP tau_within, SEXP tau_between,
  SEXP cluster_index, SEXP cluster_size, SEXP weights, SEXP gamma_grid,
  SEXP log_gamma_weights)
{
  NormClusterArgs args;
  validate_cluster_common_norm(
    yi, sei, mu_samples, tau_within, tau_between, cluster_index, cluster_size,
    weights, gamma_grid, log_gamma_weights, &args.S, &args.K, &args.G,
    &args.n_gamma
  );
  args.selected          = false;
  args.yi                = REAL(yi);
  args.sei               = REAL(sei);
  args.mu                = REAL(mu_samples);
  args.tau_within        = REAL(tau_within);
  args.tau_between       = REAL(tau_between);
  args.cluster_index     = INTEGER(cluster_index);
  args.cluster_size      = INTEGER(cluster_size);
  args.weights           = weights == R_NilValue ? 0 : REAL(weights);
  args.gamma             = REAL(gamma_grid);
  args.log_gamma         = REAL(log_gamma_weights);
  args.omega             = 0;
  args.alpha             = 0;
  args.phack_kind        = 0;
  args.kernel_mode       = 0;
  args.obs_bin           = 0;
  args.alpha_len         = 0;
  args.phack_kind_len    = 0;
  args.kernel_mode_len   = 0;

  return norm_cluster_accumulate(args);
}

extern "C" SEXP RoBMA_norm_cluster_loglik_row_sum(
  SEXP yi, SEXP sei, SEXP mu_samples, SEXP tau_within, SEXP tau_between,
  SEXP cluster_index, SEXP cluster_size, SEXP weights, SEXP gamma_grid,
  SEXP log_gamma_weights)
{
  NormClusterArgs args;
  validate_cluster_common_norm(
    yi, sei, mu_samples, tau_within, tau_between, cluster_index, cluster_size,
    weights, gamma_grid, log_gamma_weights, &args.S, &args.K, &args.G,
    &args.n_gamma
  );
  args.selected          = false;
  args.yi                = REAL(yi);
  args.sei               = REAL(sei);
  args.mu                = REAL(mu_samples);
  args.tau_within        = REAL(tau_within);
  args.tau_between       = REAL(tau_between);
  args.cluster_index     = INTEGER(cluster_index);
  args.cluster_size      = INTEGER(cluster_size);
  args.weights           = weights == R_NilValue ? 0 : REAL(weights);
  args.gamma             = REAL(gamma_grid);
  args.log_gamma         = REAL(log_gamma_weights);
  args.omega             = 0;
  args.alpha             = 0;
  args.phack_kind        = 0;
  args.kernel_mode       = 0;
  args.obs_bin           = 0;
  args.alpha_len         = 0;
  args.phack_kind_len    = 0;
  args.kernel_mode_len   = 0;

  return norm_cluster_accumulate_row_sum(args);
}

extern "C" SEXP RoBMA_selnorm_cluster_loglik(
  SEXP yi, SEXP sei, SEXP mu_samples, SEXP tau_within, SEXP tau_between,
  SEXP cluster_index, SEXP cluster_size, SEXP weights, SEXP gamma_grid,
  SEXP log_gamma_weights, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP q, SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  NormClusterArgs args;
  validate_cluster_common_norm(
    yi, sei, mu_samples, tau_within, tau_between, cluster_index, cluster_size,
    weights, gamma_grid, log_gamma_weights, &args.S, &args.K, &args.G,
    &args.n_gamma
  );
  validate_selnorm_cluster_quadrature_rules(args.n_gamma);

  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(obs_bin, "obs_bin");
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
  if (matrix_nrow(omega, "omega") != args.S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(obs_bin) != args.K) {
    Rf_error("'obs_bin' must have one value per observation.");
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

  args.selected          = true;
  args.yi                = REAL(yi);
  args.sei               = REAL(sei);
  args.mu                = REAL(mu_samples);
  args.tau_within        = REAL(tau_within);
  args.tau_between       = REAL(tau_between);
  args.cluster_index     = INTEGER(cluster_index);
  args.cluster_size      = INTEGER(cluster_size);
  args.weights           = weights == R_NilValue ? 0 : REAL(weights);
  args.gamma             = REAL(gamma_grid);
  args.log_gamma         = REAL(log_gamma_weights);
  args.omega             = REAL(omega);
  args.alpha             = REAL(alpha);
  args.phack_kind        = INTEGER(phack_kind);
  args.kernel_mode       = INTEGER(kernel_mode);
  args.obs_bin           = INTEGER(obs_bin);
  args.alpha_len         = Rf_length(alpha);
  args.phack_kind_len    = Rf_length(phack_kind);
  args.kernel_mode_len   = Rf_length(kernel_mode);

  validate_omega_matrix(args.omega, args.S, B);
  set_kernel_data(
    &args.selection_data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, LOGICAL(telescope_probabilities)[0] == TRUE
  );

  return norm_cluster_accumulate(args);
}

extern "C" SEXP RoBMA_selnorm_cluster_loglik_row_sum(
  SEXP yi, SEXP sei, SEXP mu_samples, SEXP tau_within, SEXP tau_between,
  SEXP cluster_index, SEXP cluster_size, SEXP weights, SEXP gamma_grid,
  SEXP log_gamma_weights, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP q, SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  NormClusterArgs args;
  validate_cluster_common_norm(
    yi, sei, mu_samples, tau_within, tau_between, cluster_index, cluster_size,
    weights, gamma_grid, log_gamma_weights, &args.S, &args.K, &args.G,
    &args.n_gamma
  );
  validate_selnorm_cluster_quadrature_rules(args.n_gamma);

  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(obs_bin, "obs_bin");
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
  if (matrix_nrow(omega, "omega") != args.S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(obs_bin) != args.K) {
    Rf_error("'obs_bin' must have one value per observation.");
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

  args.selected          = true;
  args.yi                = REAL(yi);
  args.sei               = REAL(sei);
  args.mu                = REAL(mu_samples);
  args.tau_within        = REAL(tau_within);
  args.tau_between       = REAL(tau_between);
  args.cluster_index     = INTEGER(cluster_index);
  args.cluster_size      = INTEGER(cluster_size);
  args.weights           = weights == R_NilValue ? 0 : REAL(weights);
  args.gamma             = REAL(gamma_grid);
  args.log_gamma         = REAL(log_gamma_weights);
  args.omega             = REAL(omega);
  args.alpha             = REAL(alpha);
  args.phack_kind        = INTEGER(phack_kind);
  args.kernel_mode       = INTEGER(kernel_mode);
  args.obs_bin           = INTEGER(obs_bin);
  args.alpha_len         = Rf_length(alpha);
  args.phack_kind_len    = Rf_length(phack_kind);
  args.kernel_mode_len   = Rf_length(kernel_mode);

  validate_omega_matrix(args.omega, args.S, B);
  set_kernel_data(
    &args.selection_data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, LOGICAL(telescope_probabilities)[0] == TRUE
  );

  return norm_cluster_accumulate_row_sum(args);
}

extern "C" SEXP RoBMA_selnorm_cluster_location_grid(
  SEXP yi, SEXP sei, SEXP mu_samples, SEXP mu_basis, SEXP current,
  SEXP values, SEXP tau_within, SEXP tau_between, SEXP cluster_index,
  SEXP cluster_size, SEXP weights, SEXP gamma_grid, SEXP log_gamma_weights,
  SEXP omega, SEXP alpha, SEXP phack_kind, SEXP kernel_mode, SEXP z_lower,
  SEXP z_upper, SEXP obs_bin, SEXP sign, SEXP q, SEXP phack_z_source,
  SEXP phack_z_dest, SEXP segment_bounds, SEXP segment_step_bin,
  SEXP segment_phack_region, SEXP telescope_probabilities)
{
  NormClusterArgs args;
  validate_cluster_common_norm(
    yi, sei, mu_samples, tau_within, tau_between, cluster_index, cluster_size,
    weights, gamma_grid, log_gamma_weights, &args.S, &args.K, &args.G,
    &args.n_gamma
  );
  validate_selnorm_cluster_quadrature_rules(args.n_gamma);

  check_real(mu_basis, "mu_basis");
  check_real(current, "current");
  check_real(values, "values");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(obs_bin, "obs_bin");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(telescope_probabilities, "telescope_probabilities");

  validate_matrix_dims_named(mu_basis, "mu_basis", args.S, args.K, "mu_samples");
  if (Rf_length(current) != args.S) {
    Rf_error("'current' must have one value per posterior sample.");
  }

  const int B = matrix_ncol(omega, "omega");
  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }
  if (matrix_nrow(omega, "omega") != args.S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(obs_bin) != args.K) {
    Rf_error("'obs_bin' must have one value per observation.");
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

  const int G_values = Rf_length(values);
  SEXP out = PROTECT(allocate_selnorm_matrix(REALSXP, G_values, args.S));
  SEXP sensitivity = PROTECT(allocate_selnorm_matrix(REALSXP, G_values, args.S));
  SEXP quadrature_order = PROTECT(allocate_selnorm_matrix(INTSXP, G_values, args.S));

  args.selected          = true;
  args.yi                = REAL(yi);
  args.sei               = REAL(sei);
  args.mu                = REAL(mu_samples);
  args.tau_within        = REAL(tau_within);
  args.tau_between       = REAL(tau_between);
  args.cluster_index     = INTEGER(cluster_index);
  args.cluster_size      = INTEGER(cluster_size);
  args.weights           = weights == R_NilValue ? 0 : REAL(weights);
  args.gamma             = REAL(gamma_grid);
  args.log_gamma         = REAL(log_gamma_weights);
  args.omega             = REAL(omega);
  args.alpha             = REAL(alpha);
  args.phack_kind        = INTEGER(phack_kind);
  args.kernel_mode       = INTEGER(kernel_mode);
  args.obs_bin           = INTEGER(obs_bin);
  args.alpha_len         = Rf_length(alpha);
  args.phack_kind_len    = Rf_length(phack_kind);
  args.kernel_mode_len   = Rf_length(kernel_mode);

  validate_omega_matrix(args.omega, args.S, B);
  set_kernel_data(
    &args.selection_data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *current_p        = REAL(current);
  const double *values_p         = REAL(values);
  const double *basis_p          = REAL(mu_basis);
  double *out_p                  = REAL(out);
  double *sensitivity_p          = REAL(sensitivity);
  int *quadrature_order_p        = INTEGER(quadrature_order);

  std::vector<double> gamma_terms(31);

  for (int s = 0; s < args.S; ++s) {
    for (int g_value = 0; g_value < G_values; ++g_value) {
      const double delta = values_p[g_value] - current_p[s];
      if (!std::isfinite(delta)) {
        out_p[g_value + G_values * s] = NA_REAL;
        sensitivity_p[g_value + G_values * s] = NA_REAL;
        quadrature_order_p[g_value + G_values * s] = NA_INTEGER;
        continue;
      }

      double row_sum = 0.0;
      double row_sensitivity = 0.0;
      int row_order = 15;
      int cluster_offset = 0;
      for (int g = 0; g < args.G; ++g) {
        const int cluster_n = args.cluster_size[g];
        const NormClusterIntegral integral = norm_cluster_adaptive_log_integral(
          args, s, g, cluster_offset, cluster_n, basis_p, delta, gamma_terms
        );
        row_sum += integral.value;
        row_sensitivity = std::max(
          row_sensitivity, integral.relative_change
        );
        row_order = std::max(row_order, integral.order);
        cluster_offset += cluster_n;
      }

      out_p[g_value + G_values * s] = row_sum;
      sensitivity_p[g_value + G_values * s] = row_sensitivity;
      quadrature_order_p[g_value + G_values * s] = row_order;
    }
  }

  Rf_setAttrib(out, Rf_install("quadrature_relative_change"), sensitivity);
  Rf_setAttrib(out, Rf_install("quadrature_order"), quadrature_order);

  UNPROTECT(3);
  return out;
}
