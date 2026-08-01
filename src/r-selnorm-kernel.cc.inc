extern "C" SEXP RoBMA_selnorm_kernel_log_norm_delta_grid(
  SEXP mean, SEXP sd, SEXP basis, SEXP current_log_norm, SEXP current,
  SEXP values, SEXP sei, SEXP weights, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(basis, "basis");
  check_real(current_log_norm, "current_log_norm");
  check_real(current, "current");
  check_real(values, "values");
  check_real(sei, "sei");
  check_real(weights, "weights");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int G = Rf_length(values);
  const SelNormNativeArgs selnorm = validate_selnorm_native_args(
    omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, telescope_probabilities, S
  );

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  validate_matrix_dims_named(basis, "basis", S, K, "mean");
  validate_matrix_dims_named(current_log_norm, "current_log_norm", S, K, "mean");

  if (Rf_length(current) != S) {
    Rf_error("'current' must have one value per posterior sample.");
  }
  if (Rf_length(sei) != K || Rf_length(weights) != K) {
    Rf_error("'sei' and 'weights' must have one value per observation.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, G, S));

  const double *mean_p             = REAL(mean);
  const double *sd_p               = REAL(sd);
  const double *basis_p            = REAL(basis);
  const double *current_log_norm_p = REAL(current_log_norm);
  const double *current_p          = REAL(current);
  const double *values_p           = REAL(values);
  const double *sei_p              = REAL(sei);
  const double *weights_p          = REAL(weights);
  double *out_p                    = REAL(out);

  for (int s = 0; s < S; ++s) {
    for (int g = 0; g < G; ++g) {
      const double delta = values_p[g] - current_p[s];
      double row_sum = 0;

      if (!std::isfinite(delta)) {
        out_p[g + G * s] = NA_REAL;
        continue;
      }

      for (int k = 0; k < K; ++k) {
        const int cell = s + S * k;
        const double candidate_log_norm = cpp_selnorm_kernel_log_norm(
          mean_p[cell] + basis_p[cell] * delta,
          sd_p[cell],
          sei_p[k],
          selnorm.omega_row(s),
          selnorm.alpha_at(s),
          selnorm.phack_kind_at(s),
          selnorm.kernel_mode_at(s),
          selnorm.data,
          S,
          false
        );

        row_sum += weights_p[k] * (candidate_log_norm - current_log_norm_p[cell]);
      }

      out_p[g + G * s] = row_sum;
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_loglik_matrix(
  SEXP yi, SEXP mu_num, SEXP sigma_num, SEXP mu_norm, SEXP sigma_norm,
  SEXP sei, SEXP weights, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP q, SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(yi, "yi");
  check_real(mu_num, "mu_num");
  check_real(sigma_num, "sigma_num");
  check_real(mu_norm, "mu_norm");
  check_real(sigma_norm, "sigma_norm");
  check_real(sei, "sei");
  check_real(weights, "weights");
  check_integer(obs_bin, "obs_bin");

  const int S = matrix_nrow(mu_num, "mu_num");
  const int K = matrix_ncol(mu_num, "mu_num");
  const SelNormNativeArgs selnorm = validate_selnorm_native_args(
    omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, telescope_probabilities, S
  );

  validate_matrix_dims(sigma_num, "sigma_num", S, K);
  validate_matrix_dims(mu_norm, "mu_norm", S, K);
  validate_matrix_dims(sigma_norm, "sigma_norm", S, K);

  if (Rf_length(yi) != K || Rf_length(sei) != K || Rf_length(weights) != K ||
      Rf_length(obs_bin) != K) {
    Rf_error("'yi', 'sei', 'weights', and 'obs_bin' must have one value per observation.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const double *yi_p         = REAL(yi);
  const double *mu_num_p     = REAL(mu_num);
  const double *sigma_num_p  = REAL(sigma_num);
  const double *mu_norm_p    = REAL(mu_norm);
  const double *sigma_norm_p = REAL(sigma_norm);
  const double *sei_p        = REAL(sei);
  const double *weights_p    = REAL(weights);
  const int *obs_bin_p       = INTEGER(obs_bin);
  double *out_p              = REAL(out);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_lpdf(
        yi_p[k],
        mu_num_p[cell],
        sigma_num_p[cell],
        mu_norm_p[cell],
        sigma_norm_p[cell],
        sei_p[k],
        weights_p[k],
        selnorm.omega_row(s),
        obs_bin_p[k],
        selnorm.alpha_at(s),
        selnorm.phack_kind_at(s),
        selnorm.kernel_mode_at(s),
        selnorm.data,
        S,
        false
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_loglik_row_sum(
  SEXP yi, SEXP mu_num, SEXP sigma_num, SEXP mu_norm, SEXP sigma_norm,
  SEXP sei, SEXP weights, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP q, SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(yi, "yi");
  check_real(mu_num, "mu_num");
  check_real(sigma_num, "sigma_num");
  check_real(mu_norm, "mu_norm");
  check_real(sigma_norm, "sigma_norm");
  check_real(sei, "sei");
  check_real(weights, "weights");
  check_integer(obs_bin, "obs_bin");

  const int S = matrix_nrow(mu_num, "mu_num");
  const int K = matrix_ncol(mu_num, "mu_num");
  const SelNormNativeArgs selnorm = validate_selnorm_native_args(
    omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, telescope_probabilities, S
  );

  validate_matrix_dims(sigma_num, "sigma_num", S, K);
  validate_matrix_dims(mu_norm, "mu_norm", S, K);
  validate_matrix_dims(sigma_norm, "sigma_norm", S, K);

  if (Rf_length(yi) != K || Rf_length(sei) != K || Rf_length(weights) != K ||
      Rf_length(obs_bin) != K) {
    Rf_error("'yi', 'sei', 'weights', and 'obs_bin' must have one value per observation.");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, S));

  const double *yi_p         = REAL(yi);
  const double *mu_num_p     = REAL(mu_num);
  const double *sigma_num_p  = REAL(sigma_num);
  const double *mu_norm_p    = REAL(mu_norm);
  const double *sigma_norm_p = REAL(sigma_norm);
  const double *sei_p        = REAL(sei);
  const double *weights_p    = REAL(weights);
  const int *obs_bin_p       = INTEGER(obs_bin);
  double *out_p              = REAL(out);

  for (int s = 0; s < S; ++s) {
    double row_sum = 0.0;
    for (int k = 0; k < K; ++k) {
      const int cell = s + S * k;
      row_sum += cpp_selnorm_kernel_lpdf(
        yi_p[k],
        mu_num_p[cell],
        sigma_num_p[cell],
        mu_norm_p[cell],
        sigma_norm_p[cell],
        sei_p[k],
        weights_p[k],
        selnorm.omega_row(s),
        obs_bin_p[k],
        selnorm.alpha_at(s),
        selnorm.phack_kind_at(s),
        selnorm.kernel_mode_at(s),
        selnorm.data,
        S,
        false
      );
    }
    out_p[s] = row_sum;
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_log_norm_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const SelNormNativeArgs selnorm = validate_selnorm_native_args(
    omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, telescope_probabilities, S
  );

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *sei_p       = REAL(sei);
  double *out_p             = REAL(out);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_log_norm(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        selnorm.omega_row(s),
        selnorm.alpha_at(s),
        selnorm.phack_kind_at(s),
        selnorm.kernel_mode_at(s),
        selnorm.data,
        S,
        false
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_cdf_matrix(
  SEXP x, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP lower_tail, SEXP telescope_probabilities)
{
  check_real(x, "q");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_logical(lower_tail, "lower.tail");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const SelNormNativeArgs selnorm = validate_selnorm_native_args(
    omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, telescope_probabilities, S
  );

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if ((Rf_length(x) != 1 && Rf_length(x) != K) || Rf_length(sei) != K) {
    Rf_error("'q' must have length 1 or K and 'sei' must have one value per observation.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const double *x_p          = REAL(x);
  const double *mean_p       = REAL(mean);
  const double *sd_p         = REAL(sd);
  const double *sei_p        = REAL(sei);
  const bool lower_tail_bool = LOGICAL(lower_tail)[0] == TRUE;
  double *out_p              = REAL(out);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_cdf(
        scalar_or_col_real(x_p, Rf_length(x), k, "q"),
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        selnorm.omega_row(s),
        selnorm.alpha_at(s),
        selnorm.phack_kind_at(s),
        selnorm.kernel_mode_at(s),
        selnorm.data,
        S,
        lower_tail_bool,
        false
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_moments_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const SelNormNativeArgs selnorm = validate_selnorm_native_args(
    omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, telescope_probabilities, S
  );

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }

  SEXP out_mean = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  SEXP out_second = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));

  SET_STRING_ELT(names, 0, Rf_mkChar("mean"));
  SET_STRING_ELT(names, 1, Rf_mkChar("second"));
  SET_VECTOR_ELT(out, 0, out_mean);
  SET_VECTOR_ELT(out, 1, out_second);
  Rf_setAttrib(out, R_NamesSymbol, names);

  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  double *mean_out_p       = REAL(out_mean);
  double *second_out_p     = REAL(out_second);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      cpp_selnorm_kernel_moments(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        selnorm.omega_row(s),
        selnorm.alpha_at(s),
        selnorm.phack_kind_at(s),
        selnorm.kernel_mode_at(s),
        selnorm.data,
        mean_out_p + cell,
        second_out_p + cell,
        S,
        false
      );
    }
  }

  UNPROTECT(4);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_rng_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const SelNormNativeArgs selnorm = validate_selnorm_native_args(
    omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, telescope_probabilities, S
  );

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  double *out_p            = REAL(out);

  GetRNGstate();
  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_rng(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        selnorm.omega_row(s),
        unif_rand(),
        unif_rand(),
        selnorm.alpha_at(s),
        selnorm.phack_kind_at(s),
        selnorm.kernel_mode_at(s),
        selnorm.data,
        S,
        false
      );
    }
  }
  PutRNGstate();

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_weighted_summary(
  SEXP yi, SEXP mean, SEXP sd, SEXP sei, SEXP psis_weights, SEXP omega,
  SEXP alpha, SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(yi, "yi");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(psis_weights, "psis_weights");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const SelNormNativeArgs selnorm = validate_selnorm_native_args(
    omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q,
    phack_z_source, phack_z_dest, segment_bounds, segment_step_bin,
    segment_phack_region, telescope_probabilities, S
  );

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  validate_matrix_dims_named(psis_weights, "psis_weights", S, K, "mean");

  if (Rf_length(yi) != K || Rf_length(sei) != K) {
    Rf_error("'yi' and 'sei' must have one value per observation.");
  }

  SEXP out_log_lower = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP out_log_upper = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP out_mean = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP out_variance = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));

  SET_STRING_ELT(names, 0, Rf_mkChar("log_lower"));
  SET_STRING_ELT(names, 1, Rf_mkChar("log_upper"));
  SET_STRING_ELT(names, 2, Rf_mkChar("mean"));
  SET_STRING_ELT(names, 3, Rf_mkChar("variance"));
  SET_VECTOR_ELT(out, 0, out_log_lower);
  SET_VECTOR_ELT(out, 1, out_log_upper);
  SET_VECTOR_ELT(out, 2, out_mean);
  SET_VECTOR_ELT(out, 3, out_variance);
  Rf_setAttrib(out, R_NamesSymbol, names);

  const double *yi_p        = REAL(yi);
  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *sei_p       = REAL(sei);
  const double *weights_p   = REAL(psis_weights);
  double *log_lower_out_p   = REAL(out_log_lower);
  double *log_upper_out_p   = REAL(out_log_upper);
  double *mean_out_p        = REAL(out_mean);
  double *variance_out_p    = REAL(out_variance);
  double *moment_mean_work = reinterpret_cast<double *>(
    R_alloc(S, sizeof(double))
  );
  double *moment_second_work = reinterpret_cast<double *>(
    R_alloc(S, sizeof(double))
  );

  for (int k = 0; k < K; ++k) {
    double log_lower_sum = -INFINITY;
    double log_upper_sum = -INFINITY;
    long double weight_sum = 0;
    long double mean_sum = 0;

    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      const double w = weights_p[cell];
      if (!std::isfinite(w) || w < 0) {
        Rf_error("'psis_weights' must be finite and non-negative.");
      }
      const double alpha_s = selnorm.alpha_at(s);
      const int phack_s = selnorm.phack_kind_at(s);
      const int mode_s = selnorm.kernel_mode_at(s);
      double moment_mean = NA_REAL;
      double moment_second = NA_REAL;

      cpp_selnorm_kernel_moments(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        selnorm.omega_row(s),
        alpha_s,
        phack_s,
        mode_s,
        selnorm.data,
        &moment_mean,
        &moment_second,
        S,
        false
      );
      if (!std::isfinite(moment_mean) || !std::isfinite(moment_second)) {
        Rf_error("Selected-normal predictive moments are not finite.");
      }
      moment_mean_work[s]   = moment_mean;
      moment_second_work[s] = moment_second;

      if (w > 0) {
        const double log_weight = std::log(w);
        const double log_lower = cpp_selnorm_kernel_log_cdf(
          yi_p[k], mean_p[cell], sd_p[cell], sei_p[k],
          selnorm.omega_row(s), alpha_s, phack_s, mode_s, selnorm.data,
          S, true, false
        );
        const double log_upper = cpp_selnorm_kernel_log_cdf(
          yi_p[k], mean_p[cell], sd_p[cell], sei_p[k],
          selnorm.omega_row(s), alpha_s, phack_s, mode_s, selnorm.data,
          S, false, false
        );

        if (std::isnan(log_lower) || std::isnan(log_upper)) {
          log_lower_sum = NA_REAL;
          log_upper_sum = NA_REAL;
        } else if (!std::isnan(log_lower_sum)) {
          const double lower_term = log_weight + log_lower;
          const double upper_term = log_weight + log_upper;
          const double lower_max  = std::max(log_lower_sum, lower_term);
          const double upper_max  = std::max(log_upper_sum, upper_term);
          log_lower_sum = lower_max == -INFINITY ? -INFINITY :
            lower_max + std::log(
              std::exp(log_lower_sum - lower_max) +
              std::exp(lower_term - lower_max)
            );
          log_upper_sum = upper_max == -INFINITY ? -INFINITY :
            upper_max + std::log(
              std::exp(log_upper_sum - upper_max) +
              std::exp(upper_term - upper_max)
            );
        }
      }
      weight_sum += static_cast<long double>(w);
      mean_sum += static_cast<long double>(w) * moment_mean;
    }

    if (!std::isfinite(weight_sum) || weight_sum <= 0) {
      Rf_error("'psis_weights' must have positive finite column sums.");
    }
    const long double mixture_mean = mean_sum / weight_sum;
    long double variance_sum = 0;
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      const long double component_mean = moment_mean_work[s];
      const long double centered = component_mean - mixture_mean;
      const long double component_variance =
        static_cast<long double>(moment_second_work[s]) -
        component_mean * component_mean;
      variance_sum += static_cast<long double>(weights_p[cell]) *
        (component_variance + centered * centered);
    }
    const long double mixture_variance = variance_sum / weight_sum;
    if (!std::isfinite(mixture_variance) || mixture_variance < 0) {
      Rf_error(
        "Selected-normal predictive variance is not finite and non-negative."
      );
    }

    const double log_weight_sum = std::log(static_cast<double>(weight_sum));
    log_lower_out_p[k] = log_lower_sum - log_weight_sum;
    log_upper_out_p[k] = log_upper_sum - log_weight_sum;
    mean_out_p[k]      = static_cast<double>(mixture_mean);
    variance_out_p[k]  = static_cast<double>(mixture_variance);
  }

  UNPROTECT(6);
  return out;
}
