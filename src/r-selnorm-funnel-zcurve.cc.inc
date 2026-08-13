extern "C" SEXP RoBMA_funnel_model_averaged_quantiles(
  SEXP se_sequence, SEXP mu, SEXP tau, SEXP PET, SEXP PEESE,
  SEXP is_weightfunction, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region, SEXP direction,
  SEXP telescope_probabilities)
{
  check_real(se_sequence, "se_sequence");
  check_real(mu, "mu");
  check_real(tau, "tau");
  check_real(PET, "PET");
  check_real(PEESE, "PEESE");
  check_integer(is_weightfunction, "is_weightfunction");
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
  check_integer(direction, "direction");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = Rf_length(mu);
  const int L = Rf_length(se_sequence);
  const int B = matrix_ncol(omega, "omega");

  if (S < 1) {
    Rf_error("'mu' must contain at least one posterior sample.");
  }
  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }
  if (Rf_length(tau) != S || Rf_length(PET) != S ||
      Rf_length(PEESE) != S || Rf_length(is_weightfunction) != S) {
    Rf_error("'mu', 'tau', 'PET', 'PEESE', and 'is_weightfunction' must have equal lengths.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(alpha) != 1 && Rf_length(alpha) != S) {
    Rf_error("'alpha' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(phack_kind) != 1 && Rf_length(phack_kind) != S) {
    Rf_error("'phack_kind' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(kernel_mode) != 1 && Rf_length(kernel_mode) != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }
  if (Rf_length(direction) != 1 ||
      (INTEGER(direction)[0] != 1 && INTEGER(direction)[0] != -1)) {
    Rf_error("'direction' must be 1 or -1.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out_lower = PROTECT(Rf_allocVector(REALSXP, L));
  SEXP out_upper = PROTECT(Rf_allocVector(REALSXP, L));
  SEXP out_mid = PROTECT(Rf_allocVector(REALSXP, L));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));

  SET_STRING_ELT(names, 0, Rf_mkChar("lower"));
  SET_STRING_ELT(names, 1, Rf_mkChar("upper"));
  SET_STRING_ELT(names, 2, Rf_mkChar("mid"));
  SET_VECTOR_ELT(out, 0, out_lower);
  SET_VECTOR_ELT(out, 1, out_upper);
  SET_VECTOR_ELT(out, 2, out_mid);
  Rf_setAttrib(out, R_NamesSymbol, names);

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *se_p      = REAL(se_sequence);
  const double *mu_p      = REAL(mu);
  const double *tau_p     = REAL(tau);
  const double *PET_p     = REAL(PET);
  const double *PEESE_p   = REAL(PEESE);
  const int *weight_p     = INTEGER(is_weightfunction);
  const double *omega_p   = REAL(omega);
  const double *alpha_p   = REAL(alpha);
  const int *phack_p      = INTEGER(phack_kind);
  const int *mode_p       = INTEGER(kernel_mode);
  const int direction_p   = INTEGER(direction)[0];
  double *lower_p         = REAL(out_lower);
  double *upper_p         = REAL(out_upper);
  double *mid_p           = REAL(out_mid);

  validate_omega_matrix(omega_p, S, B);

  std::vector<double> row_alpha(static_cast<size_t>(S));
  std::vector<int> row_phack(static_cast<size_t>(S));
  std::vector<int> row_mode(static_cast<size_t>(S));
  std::vector<bool> row_active_phack(static_cast<size_t>(S));

  for (int s = 0; s < S; ++s) {
    row_alpha[static_cast<size_t>(s)] = scalar_or_row_real(
      alpha_p, Rf_length(alpha), s, "alpha"
    );
    row_phack[static_cast<size_t>(s)] = scalar_or_row_int(
      phack_p, Rf_length(phack_kind), s, "phack_kind"
    );
    row_mode[static_cast<size_t>(s)] = scalar_or_row_int(
      mode_p, Rf_length(kernel_mode), s, "kernel_mode"
    );
    row_active_phack[static_cast<size_t>(s)] = row_has_active_phack(
      row_mode[static_cast<size_t>(s)],
      row_alpha[static_cast<size_t>(s)],
      row_phack[static_cast<size_t>(s)]
    );
  }

  FunnelSeContext ctx;
  ctx.S = S;
  ctx.direction = direction_p;
  ctx.mu = mu_p;
  ctx.tau = tau_p;
  ctx.PET = PET_p;
  ctx.PEESE = PEESE_p;
  ctx.is_weightfunction = weight_p;
  ctx.omega = omega_p;
  ctx.row_alpha = &row_alpha;
  ctx.row_phack = &row_phack;
  ctx.row_mode = &row_mode;
  ctx.row_active_phack = &row_active_phack;
  ctx.data = &data;

  for (int ell = 0; ell < L; ++ell) {
    ctx.se = se_p[ell];
    funnel_prepare_se_context(&ctx);

    const bool has_previous = ell > 0;
    lower_p[ell] = funnel_model_averaged_quantile_native(
      0.025, ctx, has_previous ? lower_p[ell - 1] : NA_REAL, has_previous
    );
    upper_p[ell] = funnel_model_averaged_quantile_native(
      0.975, ctx, has_previous ? upper_p[ell - 1] : NA_REAL, has_previous
    );
    mid_p[ell] = funnel_model_averaged_quantile_native(
      0.5, ctx, has_previous ? mid_p[ell - 1] : NA_REAL, has_previous
    );
    if (!std::isfinite(lower_p[ell]) || !std::isfinite(upper_p[ell]) ||
        !std::isfinite(mid_p[ell]) || lower_p[ell] > mid_p[ell] ||
        mid_p[ell] > upper_p[ell]) {
      Rf_error(
        "Funnel quantiles could not be computed from a valid bracketed CDF."
      );
    }
  }

  UNPROTECT(5);
  return out;
}

extern "C" SEXP RoBMA_zcurve_normal_density_matrix(
  SEXP z_sequence, SEXP mean, SEXP sd, SEXP sei)
{
  check_real(z_sequence, "z_sequence");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int L = Rf_length(z_sequence);

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, L));

  const double *z_p    = REAL(z_sequence);
  const double *mean_p = REAL(mean);
  const double *sd_p   = REAL(sd);
  const double *sei_p  = REAL(sei);
  double *out_p        = REAL(out);

  std::vector<double> log_scale(static_cast<size_t>(S) * static_cast<size_t>(K));
  for (int k = 0; k < K; ++k) {
    const double log_sei = std::log(sei_p[k]);
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      log_scale[static_cast<size_t>(cell)] =
        log_sei - std::log(sd_p[cell]);
    }
  }

  for (int ell = 0; ell < L; ++ell) {
    for (int s = 0; s < S; ++s) {
      double dens_sum = 0;
      for (int k = 0; k < K; ++k) {
        const int cell = s + S * k;
        dens_sum += std::exp(cpp_selnorm_affine_normal_lpdf_log_scale(
          z_p[ell], sei_p[k], mean_p[cell], sd_p[cell],
          log_scale[static_cast<size_t>(cell)]
        ));
      }
      out_p[s + S * ell] = dens_sum / K;
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_zcurve_threshold_summary(
  SEXP z_threshold, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP extrapolate, SEXP telescope_probabilities)
{
  check_real(z_threshold, "z_threshold");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
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
  check_logical(extrapolate, "extrapolate");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (Rf_length(z_threshold) != 1) {
    Rf_error("'z_threshold' must have length 1.");
  }
  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
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
  if (Rf_length(alpha) != 1 && Rf_length(alpha) != S) {
    Rf_error("'alpha' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(phack_kind) != 1 && Rf_length(phack_kind) != S) {
    Rf_error("'phack_kind' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(kernel_mode) != 1 && Rf_length(kernel_mode) != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }
  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out_edr = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP out_weights = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));

  SET_STRING_ELT(names, 0, Rf_mkChar("EDR"));
  SET_STRING_ELT(names, 1, Rf_mkChar("weights"));
  SET_VECTOR_ELT(out, 0, out_edr);
  SET_VECTOR_ELT(out, 1, out_weights);
  Rf_setAttrib(out, R_NamesSymbol, names);

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double threshold = REAL(z_threshold)[0];
  const double *mean_p   = REAL(mean);
  const double *sd_p     = REAL(sd);
  const double *sei_p    = REAL(sei);
  const double *omega_p  = REAL(omega);
  const double *alpha_p  = REAL(alpha);
  const int *phack_p     = INTEGER(phack_kind);
  const int *mode_p      = INTEGER(kernel_mode);
  const bool extrapolate_b = LOGICAL(extrapolate)[0] == TRUE;
  double *edr_p          = REAL(out_edr);
  double *weights_p      = REAL(out_weights);

  validate_omega_matrix(omega_p, S, B);

  std::vector<double> row_alpha(S);
  std::vector<int> row_phack(S);
  std::vector<int> row_mode(S);
  std::vector<bool> row_active_phack(S);
  std::vector<bool> row_use_step(S);

  for (int s = 0; s < S; ++s) {
    row_alpha[static_cast<size_t>(s)] = scalar_or_row_real(
      alpha_p, Rf_length(alpha), s, "alpha"
    );
    row_phack[static_cast<size_t>(s)] = scalar_or_row_int(
      phack_p, Rf_length(phack_kind), s, "phack_kind"
    );
    row_mode[static_cast<size_t>(s)] = scalar_or_row_int(
      mode_p, Rf_length(kernel_mode), s, "kernel_mode"
    );
    row_active_phack[static_cast<size_t>(s)] = row_has_active_phack(
      row_mode[static_cast<size_t>(s)],
      row_alpha[static_cast<size_t>(s)],
      row_phack[static_cast<size_t>(s)]
    );
    row_use_step[static_cast<size_t>(s)] =
      row_mode[static_cast<size_t>(s)] == SELKERNEL_STEP ||
      row_mode[static_cast<size_t>(s)] == SELKERNEL_STEP_PHACK_POWER;
  }

  for (int s = 0; s < S; ++s) {
    double threshold_sum = 0;
    double weight_sum = 0;

    for (int k = 0; k < K; ++k) {
      const int cell = s + S * k;
      double inverse_weight = 1;
      double threshold_prob;

      if (row_active_phack[static_cast<size_t>(s)]) {
        threshold_prob = NA_REAL;
        inverse_weight = NA_REAL;
      } else if (row_use_step[static_cast<size_t>(s)] && extrapolate_b) {
        double unused_inverse;
        threshold_prob = cpp_selnorm_kernel_threshold(
          threshold,
          mean_p[cell],
          sd_p[cell],
          sei_p[k],
          omega_p + s,
          row_alpha[static_cast<size_t>(s)],
          row_phack[static_cast<size_t>(s)],
          SELKERNEL_NORMAL,
          data,
          &unused_inverse,
          S,
          false
        );

        const double log_norm = cpp_selnorm_kernel_log_norm(
          mean_p[cell],
          sd_p[cell],
          sei_p[k],
          omega_p + s,
          row_alpha[static_cast<size_t>(s)],
          row_phack[static_cast<size_t>(s)],
          row_mode[static_cast<size_t>(s)],
          data,
          S,
          false
        );
        inverse_weight = std::isfinite(log_norm) ?
          std::exp(-log_norm) : NA_REAL;
      } else {
        threshold_prob = cpp_selnorm_kernel_threshold(
          threshold,
          mean_p[cell],
          sd_p[cell],
          sei_p[k],
          omega_p + s,
          row_alpha[static_cast<size_t>(s)],
          row_phack[static_cast<size_t>(s)],
          row_mode[static_cast<size_t>(s)],
          data,
          &inverse_weight,
          S,
          false
        );
      }

      if (extrapolate_b) {
        threshold_sum += threshold_prob * inverse_weight;
      } else {
        threshold_sum += threshold_prob;
      }
      weight_sum += inverse_weight;
    }

    edr_p[s] = extrapolate_b ? threshold_sum / weight_sum : threshold_sum / K;
    weights_p[s] = weight_sum / K;
  }

  UNPROTECT(4);
  return out;
}

extern "C" SEXP RoBMA_selnorm_zcurve_density_matrix(
  SEXP z_sequence, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP extrapolate, SEXP telescope_probabilities)
{
  check_real(z_sequence, "z_sequence");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
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
  if (TYPEOF(extrapolate) != LGLSXP) {
    Rf_error("'extrapolate' must be logical.");
  }
  check_logical(telescope_probabilities, "telescope_probabilities");

  const bool paired = Rf_length(extrapolate) == 2;
  if (Rf_length(extrapolate) != 1 && !paired) {
    Rf_error("'extrapolate' must have length one or request the fitted/extrapolated pair.");
  }
  if (paired &&
      (LOGICAL(extrapolate)[0] != FALSE || LOGICAL(extrapolate)[1] != TRUE)) {
    Rf_error("A paired density request must be c(FALSE, TRUE).");
  }
  if (!paired && LOGICAL(extrapolate)[0] == NA_LOGICAL) {
    Rf_error("'extrapolate' must not be missing.");
  }

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int L = Rf_length(z_sequence);
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
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
  if (Rf_length(alpha) != 1 && Rf_length(alpha) != S) {
    Rf_error("'alpha' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(phack_kind) != 1 && Rf_length(phack_kind) != S) {
    Rf_error("'phack_kind' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(kernel_mode) != 1 && Rf_length(kernel_mode) != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  int n_protect = 0;
  SEXP out          = R_NilValue;
  SEXP fitted       = R_NilValue;
  SEXP extrapolated = R_NilValue;
  double *out_p          = NULL;
  double *fitted_p       = NULL;
  double *extrapolated_p = NULL;

  if (paired) {
    fitted = PROTECT(Rf_allocMatrix(REALSXP, S, L));
    ++n_protect;
    extrapolated = PROTECT(Rf_allocMatrix(REALSXP, S, L));
    ++n_protect;
    out = PROTECT(Rf_allocVector(VECSXP, 2));
    ++n_protect;
    SET_VECTOR_ELT(out, 0, fitted);
    SET_VECTOR_ELT(out, 1, extrapolated);

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    ++n_protect;
    SET_STRING_ELT(names, 0, Rf_mkChar("fitted"));
    SET_STRING_ELT(names, 1, Rf_mkChar("extrapolated"));
    Rf_setAttrib(out, R_NamesSymbol, names);

    fitted_p       = REAL(fitted);
    extrapolated_p = REAL(extrapolated);
  } else {
    out = PROTECT(Rf_allocMatrix(REALSXP, S, L));
    ++n_protect;
    out_p = REAL(out);
  }

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *z_p        = REAL(z_sequence);
  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  const double *omega_p    = REAL(omega);
  const double *alpha_p    = REAL(alpha);
  const int *phack_kind_p  = INTEGER(phack_kind);
  const int *kernel_mode_p = INTEGER(kernel_mode);
  const bool extrapolate_b = LOGICAL(extrapolate)[0] == TRUE;

  validate_omega_matrix(omega_p, S, B);

  std::vector<double> row_alpha(S);
  std::vector<int> row_phack(S);
  std::vector<int> row_mode(S);
  std::vector<bool> row_active_phack(S);
  std::vector<bool> row_use_step(S);
  bool any_step = false;

  for (int s = 0; s < S; ++s) {
    row_alpha[static_cast<size_t>(s)] = scalar_or_row_real(
      alpha_p, Rf_length(alpha), s, "alpha"
    );
    row_phack[static_cast<size_t>(s)] = scalar_or_row_int(
      phack_kind_p, Rf_length(phack_kind), s, "phack_kind"
    );
    row_mode[static_cast<size_t>(s)] = scalar_or_row_int(
      kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"
    );
    row_active_phack[static_cast<size_t>(s)] = row_has_active_phack(
      row_mode[static_cast<size_t>(s)],
      row_alpha[static_cast<size_t>(s)],
      row_phack[static_cast<size_t>(s)]
    );
    row_use_step[static_cast<size_t>(s)] =
      row_mode[static_cast<size_t>(s)] == SELKERNEL_STEP ||
      row_mode[static_cast<size_t>(s)] == SELKERNEL_STEP_PHACK_POWER;
    any_step = any_step || row_use_step[static_cast<size_t>(s)];
  }

  std::vector<double> log_norm(S * K, 0);
  std::vector<double> log_scale(S * K);
  for (int k = 0; k < K; ++k) {
    const double log_sei = std::log(sei_p[k]);
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      log_scale[static_cast<size_t>(cell)] =
        log_sei - std::log(sd_p[cell]);

      if (row_active_phack[static_cast<size_t>(s)]) {
        log_norm[cell] = std::numeric_limits<double>::quiet_NaN();
      } else if (row_use_step[static_cast<size_t>(s)]) {
        log_norm[cell] = cpp_selnorm_kernel_log_norm(
          mean_p[cell],
          sd_p[cell],
          sei_p[k],
          omega_p + s,
          row_alpha[static_cast<size_t>(s)],
          row_phack[static_cast<size_t>(s)],
          row_mode[static_cast<size_t>(s)],
          data,
          S,
          false
        );
      }
    }
  }

  std::vector<int> z_bin(L, 1);
  if (any_step && (!extrapolate_b || paired)) {
    for (int ell = 0; ell < L; ++ell) {
      z_bin[static_cast<size_t>(ell)] = z_step_bin(
        data.effect_sign * z_p[ell], data
      );
    }
  }

  for (int ell = 0; ell < L; ++ell) {
    for (int s = 0; s < S; ++s) {
      if (row_active_phack[static_cast<size_t>(s)]) {
        const double missing = std::numeric_limits<double>::quiet_NaN();
        if (paired) {
          fitted_p[s + S * ell]       = missing;
          extrapolated_p[s + S * ell] = missing;
        } else {
          out_p[s + S * ell] = missing;
        }
        continue;
      }

      if (paired) {
        double fitted_sum       = 0;
        double extrapolated_sum = 0;

        if (row_use_step[static_cast<size_t>(s)]) {
          const int bin  = z_bin[static_cast<size_t>(ell)];
          const double w = omega_p[s + S * (bin - 1)];
          const bool positive_weight = w > 0;
          const double log_observed_weight = positive_weight ? std::log(w) : 0;

          for (int k = 0; k < K; ++k) {
            const int cell = s + S * k;
            double log_density = cpp_selnorm_affine_normal_lpdf_log_scale(
              z_p[ell], sei_p[k], mean_p[cell], sd_p[cell],
              log_scale[static_cast<size_t>(cell)]
            );
            log_density -= log_norm[cell];
            extrapolated_sum += std::exp(log_density);
            if (positive_weight) {
              fitted_sum += std::exp(log_density + log_observed_weight);
            }
          }
        } else {
          for (int k = 0; k < K; ++k) {
            const int cell = s + S * k;
            fitted_sum += std::exp(cpp_selnorm_affine_normal_lpdf_log_scale(
              z_p[ell], sei_p[k], mean_p[cell], sd_p[cell],
              log_scale[static_cast<size_t>(cell)]
            ));
          }
          extrapolated_sum = fitted_sum;
        }

        fitted_p[s + S * ell]       = fitted_sum / K;
        extrapolated_p[s + S * ell] = extrapolated_sum / K;
        continue;
      }

      double dens_sum = 0;
      if (row_use_step[static_cast<size_t>(s)]) {
        double log_observed_weight = 0;

        if (!extrapolate_b) {
          const int bin = z_bin[static_cast<size_t>(ell)];
          const double w = omega_p[s + S * (bin - 1)];
          if (!(w > 0)) {
            out_p[s + S * ell] = 0;
            continue;
          }
          log_observed_weight = std::log(w);
        }

        for (int k = 0; k < K; ++k) {
          const int cell = s + S * k;
          double log_density = cpp_selnorm_affine_normal_lpdf_log_scale(
            z_p[ell], sei_p[k], mean_p[cell], sd_p[cell],
            log_scale[static_cast<size_t>(cell)]
          );

          log_density -= log_norm[cell];
          if (!extrapolate_b) {
            log_density += log_observed_weight;
          }

          dens_sum += std::exp(log_density);
        }
      } else {
        for (int k = 0; k < K; ++k) {
          const int cell = s + S * k;
          dens_sum += std::exp(cpp_selnorm_affine_normal_lpdf_log_scale(
            z_p[ell], sei_p[k], mean_p[cell], sd_p[cell],
            log_scale[static_cast<size_t>(cell)]
          ));
        }
      }

      out_p[s + S * ell] = dens_sum / K;
    }
  }

  UNPROTECT(n_protect);
  return out;
}
