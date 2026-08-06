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
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, G_values, args.S));
  SEXP sensitivity = PROTECT(Rf_allocMatrix(REALSXP, G_values, args.S));
  SEXP quadrature_order = PROTECT(Rf_allocMatrix(INTSXP, G_values, args.S));

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
