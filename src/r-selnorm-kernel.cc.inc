namespace {

// Row-wise selection arguments materialized before a parallel region. The row
// accessors of SelNormNativeArgs re-check their lengths and can raise an R
// error, which a worker thread must never do; the lengths are already
// validated by validate_selnorm_native_args(), so this only moves the lookup.
struct SelNormKernelRowArgs {
  std::vector<double> alpha;
  std::vector<int>    phack_kind;
  std::vector<int>    kernel_mode;

  SelNormKernelRowArgs(const SelNormNativeArgs &args, int S)
    : alpha(static_cast<std::size_t>(S)),
      phack_kind(static_cast<std::size_t>(S)),
      kernel_mode(static_cast<std::size_t>(S))
  {
    robma_parallel::SerialWork work;
    for (int s = 0; s < S; ++s) {
      work.add();
      alpha[static_cast<std::size_t>(s)]       = args.alpha_at(s);
      phack_kind[static_cast<std::size_t>(s)]  = args.phack_kind_at(s);
      kernel_mode[static_cast<std::size_t>(s)] = args.kernel_mode_at(s);
    }
  }
};

}

extern "C" SEXP RoBMA_selnorm_kernel_log_norm_delta_grid(
  SEXP mean, SEXP sd, SEXP basis, SEXP current_log_norm, SEXP current,
  SEXP values, SEXP sei, SEXP weights, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  ROBMA_NATIVE_BEGIN(mean, sd, basis, current_log_norm, current, values, sei, weights, omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q, phack_z_source, phack_z_dest, segment_bounds, segment_step_bin, segment_phack_region, telescope_probabilities)
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

  SEXP out = PROTECT(allocate_selnorm_matrix(REALSXP, G, S));

  const double *mean_p             = REAL(mean);
  const double *sd_p               = REAL(sd);
  const double *basis_p            = REAL(basis);
  const double *current_log_norm_p = REAL(current_log_norm);
  const double *current_p          = REAL(current);
  const double *values_p           = REAL(values);
  const double *sei_p              = REAL(sei);
  const double *weights_p          = REAL(weights);
  double *out_p                    = REAL(out);

  const SelNormKernelRowArgs row_args(selnorm, S);
  // Posterior rows are independent: row s reads call-common inputs and its own
  // slice of the inputs, and writes only its own G output slots. The per-row
  // arithmetic and its order are the ones the serial loop uses, so the result
  // is identical at any thread count.
  const RobmaRowSchedule schedule = robma_row_schedule(
    S, static_cast<double>(G) * K * selnorm.B);
  const auto compute_row = [&](int s) {
    robma_parallel::SerialWork row_work(schedule.threads <= 1);
    const std::size_t row = static_cast<std::size_t>(s);
    for (int g = 0; g < G; ++g) {
      row_work.add();
      const double delta = values_p[g] - current_p[s];
      double row_sum = 0;

      if (!std::isfinite(delta)) {
        out_p[g + G * s] = NA_REAL;
        continue;
      }

      for (int k = 0; k < K; ++k) {
        row_work.add(selnorm.B);
        const int cell = s + S * k;
        const double candidate_log_norm = cpp_selnorm_kernel_log_norm(
          mean_p[cell] + basis_p[cell] * delta,
          sd_p[cell],
          sei_p[k],
          selnorm.omega_row(s),
          row_args.alpha[row],
          row_args.phack_kind[row],
          row_args.kernel_mode[row],
          selnorm.data,
          S,
          false
        );

        row_sum += weights_p[k] * (candidate_log_norm - current_log_norm_p[cell]);
      }

      out_p[g + G * s] = row_sum;
    }
  };

  // One log-normalizer evaluation per grid value, observation and bin.
  robma_parallel::for_rows(S, schedule, compute_row);

  UNPROTECT(1);
  return out;
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_kernel_loglik_matrix(
  SEXP yi, SEXP mu_num, SEXP sigma_num, SEXP mu_norm, SEXP sigma_norm,
  SEXP sei, SEXP weights, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP q, SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  ROBMA_NATIVE_BEGIN(yi, mu_num, sigma_num, mu_norm, sigma_norm, sei, weights, omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, obs_bin, sign, q, phack_z_source, phack_z_dest, segment_bounds, segment_step_bin, segment_phack_region, telescope_probabilities)
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

  SEXP out = PROTECT(allocate_selnorm_matrix(REALSXP, S, K));

  const double *yi_p         = REAL(yi);
  const double *mu_num_p     = REAL(mu_num);
  const double *sigma_num_p  = REAL(sigma_num);
  const double *mu_norm_p    = REAL(mu_norm);
  const double *sigma_norm_p = REAL(sigma_norm);
  const double *sei_p        = REAL(sei);
  const double *weights_p    = REAL(weights);
  const int *obs_bin_p       = INTEGER(obs_bin);
  double *out_p              = REAL(out);

  const SelNormKernelRowArgs row_args(selnorm, S);
  const auto compute_cell = [&](int s, int k) {
    const int cell = s + S * k;
    const std::size_t row = static_cast<std::size_t>(s);
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
      row_args.alpha[row],
      row_args.phack_kind[row],
      row_args.kernel_mode[row],
      selnorm.data,
      S,
      false
    );
  };

  // Every (row, observation) cell is an independent evaluation, so visiting the
  // cells row-first gives each posterior row its own disjoint output slots
  // without changing any value.
  const RobmaRowSchedule schedule = robma_row_schedule(
    S, static_cast<double>(K) * selnorm.B);
  if (schedule.threads == 1) {
    robma_parallel::SerialWork work;
    for (int k = 0; k < K; ++k) {
      for (int s = 0; s < S; ++s) {
        compute_cell(s, k);
        work.add(selnorm.B);
      }
    }
  } else {
    robma_parallel::for_rows(S, schedule, [&](int s) {
      for (int k = 0; k < K; ++k) compute_cell(s, k);
    });
  }

  UNPROTECT(1);
  return out;
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_kernel_loglik_row_sum(
  SEXP yi, SEXP mu_num, SEXP sigma_num, SEXP mu_norm, SEXP sigma_norm,
  SEXP sei, SEXP weights, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP q, SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  ROBMA_NATIVE_BEGIN(yi, mu_num, sigma_num, mu_norm, sigma_norm, sei, weights, omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, obs_bin, sign, q, phack_z_source, phack_z_dest, segment_bounds, segment_step_bin, segment_phack_region, telescope_probabilities)
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

  const SelNormKernelRowArgs row_args(selnorm, S);
  // The row sum keeps its serial accumulation order; only whole rows are
  // distributed, and each row writes one output slot of its own.
  const RobmaRowSchedule schedule = robma_row_schedule(
    S, static_cast<double>(K) * selnorm.B);
  const auto compute_row = [&](int s) {
    robma_parallel::SerialWork row_work(schedule.threads <= 1);
    const std::size_t row = static_cast<std::size_t>(s);
    double row_sum = 0.0;
    for (int k = 0; k < K; ++k) {
      row_work.add(selnorm.B);
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
        row_args.alpha[row],
        row_args.phack_kind[row],
        row_args.kernel_mode[row],
        selnorm.data,
        S,
        false
      );
    }
    out_p[s] = row_sum;
  };

  robma_parallel::for_rows(S, schedule, compute_row);

  UNPROTECT(1);
  return out;
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_kernel_log_norm_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  ROBMA_NATIVE_BEGIN(mean, sd, sei, omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q, phack_z_source, phack_z_dest, segment_bounds, segment_step_bin, segment_phack_region, telescope_probabilities)
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

  SEXP out = PROTECT(allocate_selnorm_matrix(REALSXP, S, K));

  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *sei_p       = REAL(sei);
  double *out_p             = REAL(out);

  const SelNormKernelRowArgs row_args(selnorm, S);
  const auto compute_cell = [&](int s, int k) {
    const int cell = s + S * k;
    const std::size_t row = static_cast<std::size_t>(s);
    out_p[cell] = cpp_selnorm_kernel_log_norm(
      mean_p[cell],
      sd_p[cell],
      sei_p[k],
      selnorm.omega_row(s),
      row_args.alpha[row],
      row_args.phack_kind[row],
      row_args.kernel_mode[row],
      selnorm.data,
      S,
      false
    );
  };

  const RobmaRowSchedule schedule = robma_row_schedule(
    S, static_cast<double>(K) * selnorm.B);
  if (schedule.threads == 1) {
    robma_parallel::SerialWork work;
    for (int k = 0; k < K; ++k) {
      for (int s = 0; s < S; ++s) {
        compute_cell(s, k);
        work.add(selnorm.B);
      }
    }
  } else {
    robma_parallel::for_rows(S, schedule, [&](int s) {
      for (int k = 0; k < K; ++k) compute_cell(s, k);
    });
  }

  UNPROTECT(1);
  return out;
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_kernel_cdf_matrix(
  SEXP x, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP lower_tail, SEXP telescope_probabilities)
{
  ROBMA_NATIVE_BEGIN(x, mean, sd, sei, omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q, phack_z_source, phack_z_dest, segment_bounds, segment_step_bin, segment_phack_region, lower_tail, telescope_probabilities)
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

  SEXP out = PROTECT(allocate_selnorm_matrix(REALSXP, S, K));

  const double *x_p          = REAL(x);
  const double *mean_p       = REAL(mean);
  const double *sd_p         = REAL(sd);
  const double *sei_p        = REAL(sei);
  const bool lower_tail_bool = LOGICAL(lower_tail)[0] == TRUE;
  double *out_p              = REAL(out);

  const SelNormKernelRowArgs row_args(selnorm, S);
  // Recycling 'q' over observations is resolved here so no length lookup, and
  // therefore no possible R error, remains inside the row body.
  std::vector<double> q_by_observation(static_cast<std::size_t>(K));
  for (int k = 0; k < K; ++k) {
    q_by_observation[static_cast<std::size_t>(k)] =
      scalar_or_col_real(x_p, Rf_length(x), k, "q");
  }
  const auto compute_cell = [&](int s, int k) {
    const int cell = s + S * k;
    const std::size_t row = static_cast<std::size_t>(s);
    out_p[cell] = cpp_selnorm_kernel_cdf(
      q_by_observation[static_cast<std::size_t>(k)],
      mean_p[cell],
      sd_p[cell],
      sei_p[k],
      selnorm.omega_row(s),
      row_args.alpha[row],
      row_args.phack_kind[row],
      row_args.kernel_mode[row],
      selnorm.data,
      S,
      lower_tail_bool,
      false
    );
  };

  const RobmaRowSchedule schedule = robma_row_schedule(
    S, static_cast<double>(K) * selnorm.B);
  if (schedule.threads == 1) {
    robma_parallel::SerialWork work;
    for (int k = 0; k < K; ++k) {
      for (int s = 0; s < S; ++s) {
        compute_cell(s, k);
        work.add(selnorm.B);
      }
    }
  } else {
    robma_parallel::for_rows(S, schedule, [&](int s) {
      for (int k = 0; k < K; ++k) compute_cell(s, k);
    });
  }

  UNPROTECT(1);
  return out;
  ROBMA_NATIVE_END
}

// Raised when a cell's moments are not representable in double precision
// (cpp_selnorm_kernel_predictive_summary()); the R product-weight path raises
// the same message for the NaN moments of the matrix summary.
static const char *const selnorm_unrepresentable_moments_message =
  "Selected-normal predictive moments are unavailable because the selection "
  "weights span beyond double precision.";

// Paired log tails log P(Y <= q) and log P(Y >= q) with the mean and
// variance of every cell, from cpp_selnorm_kernel_predictive_summary(). Both
// tails are kept so that callers resolve whichever is small without a
// complement.
extern "C" SEXP RoBMA_selnorm_kernel_predictive_summary_matrix(
  SEXP x, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  ROBMA_NATIVE_BEGIN(x, mean, sd, sei, omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q, phack_z_source, phack_z_dest, segment_bounds, segment_step_bin, segment_phack_region, telescope_probabilities)
  check_real(x, "q");
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

  if ((Rf_length(x) != 1 && Rf_length(x) != K) || Rf_length(sei) != K) {
    Rf_error("'q' must have length 1 or K and 'sei' must have one value per observation.");
  }

  const char *labels[] = {"log_lower", "log_upper", "mean", "variance"};
  SEXP out   = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
  double *out_p[4];
  for (int i = 0; i < 4; ++i) {
    SET_VECTOR_ELT(out, i, allocate_selnorm_matrix(REALSXP, S, K));
    SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
    out_p[i] = REAL(VECTOR_ELT(out, i));
  }
  Rf_setAttrib(out, R_NamesSymbol, names);

  const double *x_p    = REAL(x);
  const double *mean_p = REAL(mean);
  const double *sd_p   = REAL(sd);
  const double *sei_p  = REAL(sei);

  const SelNormKernelRowArgs row_args(selnorm, S);
  std::vector<double> q_by_observation(static_cast<std::size_t>(K));
  for (int k = 0; k < K; ++k) {
    q_by_observation[static_cast<std::size_t>(k)] =
      scalar_or_col_real(x_p, Rf_length(x), k, "q");
  }
  const auto compute_cell = [&](int s, int k) {
    const int cell = s + S * k;
    const std::size_t row = static_cast<std::size_t>(s);
    cpp_selnorm_kernel_predictive_summary(
      q_by_observation[static_cast<std::size_t>(k)],
      mean_p[cell],
      sd_p[cell],
      sei_p[k],
      selnorm.omega_row(s),
      row_args.alpha[row],
      row_args.phack_kind[row],
      row_args.kernel_mode[row],
      selnorm.data,
      out_p[0] + cell,
      out_p[1] + cell,
      out_p[2] + cell,
      out_p[3] + cell,
      S
    );
  };

  const RobmaRowSchedule schedule = robma_row_schedule(
    S, static_cast<double>(K) * selnorm.B);
  if (schedule.threads == 1) {
    robma_parallel::SerialWork work;
    for (int k = 0; k < K; ++k) {
      for (int s = 0; s < S; ++s) {
        compute_cell(s, k);
        work.add(selnorm.B);
      }
    }
  } else {
    robma_parallel::for_rows(S, schedule, [&](int s) {
      for (int k = 0; k < K; ++k) compute_cell(s, k);
    });
  }

  UNPROTECT(2);
  return out;
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_kernel_moments_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  ROBMA_NATIVE_BEGIN(mean, sd, sei, omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q, phack_z_source, phack_z_dest, segment_bounds, segment_step_bin, segment_phack_region, telescope_probabilities)
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

  SEXP out_mean = PROTECT(allocate_selnorm_matrix(REALSXP, S, K));
  SEXP out_second = PROTECT(allocate_selnorm_matrix(REALSXP, S, K));
  SEXP out_variance = PROTECT(allocate_selnorm_matrix(REALSXP, S, K));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));

  SET_STRING_ELT(names, 0, Rf_mkChar("mean"));
  SET_STRING_ELT(names, 1, Rf_mkChar("second"));
  SET_STRING_ELT(names, 2, Rf_mkChar("variance"));
  SET_VECTOR_ELT(out, 0, out_mean);
  SET_VECTOR_ELT(out, 1, out_second);
  SET_VECTOR_ELT(out, 2, out_variance);
  Rf_setAttrib(out, R_NamesSymbol, names);

  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  double *mean_out_p       = REAL(out_mean);
  double *second_out_p     = REAL(out_second);
  double *variance_out_p   = REAL(out_variance);

  const SelNormKernelRowArgs row_args(selnorm, S);
  const auto compute_cell = [&](int s, int k) {
    const int cell = s + S * k;
    const std::size_t row = static_cast<std::size_t>(s);
    cpp_selnorm_kernel_moments(
      mean_p[cell],
      sd_p[cell],
      sei_p[k],
      selnorm.omega_row(s),
      row_args.alpha[row],
      row_args.phack_kind[row],
      row_args.kernel_mode[row],
      selnorm.data,
      mean_out_p + cell,
      second_out_p + cell,
      S,
      false,
      variance_out_p + cell
    );
  };

  const RobmaRowSchedule schedule = robma_row_schedule(
    S, static_cast<double>(K) * selnorm.B);
  if (schedule.threads == 1) {
    robma_parallel::SerialWork work;
    for (int k = 0; k < K; ++k) {
      for (int s = 0; s < S; ++s) {
        compute_cell(s, k);
        work.add(selnorm.B);
      }
    }
  } else {
    robma_parallel::for_rows(S, schedule, [&](int s) {
      for (int k = 0; k < K; ++k) compute_cell(s, k);
    });
  }

  UNPROTECT(5);
  return out;
  ROBMA_NATIVE_END
}

extern "C" SEXP RoBMA_selnorm_kernel_rng_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  ROBMA_NATIVE_BEGIN(mean, sd, sei, omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q, phack_z_source, phack_z_dest, segment_bounds, segment_step_bin, segment_phack_region, telescope_probabilities)
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

  SEXP out = PROTECT(allocate_selnorm_matrix(REALSXP, S, K));

  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  double *out_p            = REAL(out);

  // Serial: the draws consume R's random number stream in a fixed order, so
  // distributing rows would change the values and touch the R API in workers.
  robma_parallel::SerialWork work;
  GetRNGstate();
  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      const double u_bin = unif_rand();
      const double u_interval = unif_rand();
      out_p[cell] = cpp_selnorm_kernel_rng(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        selnorm.omega_row(s),
        u_bin,
        u_interval,
        selnorm.alpha_at(s),
        selnorm.phack_kind_at(s),
        selnorm.kernel_mode_at(s),
        selnorm.data,
        S,
        false
      );
      work.add(selnorm.B);
    }
  }
  PutRNGstate();

  UNPROTECT(1);
  return out;
  ROBMA_NATIVE_END
}

// Log of a sum of exponentials accumulated one term at a time against the
// running maximum: one exponential per term instead of a log and two.
struct StreamingLogSum {
  double maximum = -INFINITY;
  double scaled  = 0;

  void add(double term)
  {
    if (term == -INFINITY) {
      return;
    }
    if (term > maximum) {
      scaled  = maximum == -INFINITY ? 1 : scaled * std::exp(maximum - term) + 1;
      maximum = term;
    } else {
      scaled += std::exp(term - maximum);
    }
  }

  double value() const
  {
    return maximum == -INFINITY ? -INFINITY : maximum + std::log(scaled);
  }
};

extern "C" SEXP RoBMA_selnorm_kernel_weighted_summary(
  SEXP yi, SEXP mean, SEXP sd, SEXP sei, SEXP psis_weights, SEXP omega,
  SEXP alpha, SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  ROBMA_NATIVE_BEGIN(yi, mean, sd, sei, psis_weights, omega, alpha, phack_kind, kernel_mode, z_lower, z_upper, sign, q, phack_z_source, phack_z_dest, segment_bounds, segment_step_bin, segment_phack_region, telescope_probabilities)
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

  // Observations are independent reductions over the posterior rows, so they
  // run in parallel; each keeps the serial row order, and a failure is
  // recorded where the serial loop would stop and raised afterwards on the
  // main thread, in observation order.
  enum SummaryStatus {
    summary_ok = 0,
    summary_invalid_weight,
    summary_unavailable_moments,
    summary_invalid_moments,
    summary_invalid_weight_sum,
    summary_invalid_variance,
    summary_undefined_tails
  };
  std::vector<int> status(static_cast<std::size_t>(K), summary_ok);

  const SelNormKernelRowArgs row_args(selnorm, S);
  const RobmaRowSchedule schedule = robma_row_schedule(
    K, static_cast<double>(S) * selnorm.B);
  const auto summarize = [&](int k) {
    // Serial execution polls interrupts within an observation's rows too.
    robma_parallel::SerialWork row_work(schedule.threads <= 1);
    std::vector<double> moment_mean_work(static_cast<std::size_t>(S));
    std::vector<double> moment_variance_work(static_cast<std::size_t>(S));
    StreamingLogSum lower_sum;
    StreamingLogSum upper_sum;
    bool tails_unavailable = false;
    long double weight_sum = 0;
    long double mean_sum = 0;

    for (int s = 0; s < S; ++s) {
      row_work.add(selnorm.B);
      const int cell = s + S * k;
      const double w = weights_p[cell];
      if (!std::isfinite(w) || w < 0) {
        status[static_cast<std::size_t>(k)] = summary_invalid_weight;
        return;
      }
      if (w == 0) {
        continue;
      }
      double log_lower;
      double log_upper;
      double moment_mean;
      double moment_variance;
      bool moments_unavailable;
      cpp_selnorm_kernel_predictive_summary(
        yi_p[k], mean_p[cell], sd_p[cell], sei_p[k], selnorm.omega_row(s),
        row_args.alpha[static_cast<std::size_t>(s)],
        row_args.phack_kind[static_cast<std::size_t>(s)],
        row_args.kernel_mode[static_cast<std::size_t>(s)], selnorm.data,
        &log_lower, &log_upper, &moment_mean, &moment_variance, S,
        &moments_unavailable
      );
      if (moments_unavailable) {
        status[static_cast<std::size_t>(k)] = summary_unavailable_moments;
        return;
      }
      if (!std::isfinite(moment_mean) || !std::isfinite(moment_variance) ||
          moment_variance < 0) {
        status[static_cast<std::size_t>(k)] = summary_invalid_moments;
        return;
      }
      moment_mean_work[static_cast<std::size_t>(s)]     = moment_mean;
      moment_variance_work[static_cast<std::size_t>(s)] = moment_variance;

      if (std::isnan(log_lower) || std::isnan(log_upper)) {
        tails_unavailable = true;
      } else if (!tails_unavailable) {
        const double log_weight = std::log(w);
        lower_sum.add(log_weight + log_lower);
        upper_sum.add(log_weight + log_upper);
      }
      weight_sum += static_cast<long double>(w);
      mean_sum += static_cast<long double>(w) * moment_mean;
    }

    if (!std::isfinite(weight_sum) || weight_sum <= 0) {
      status[static_cast<std::size_t>(k)] = summary_invalid_weight_sum;
      return;
    }
    const long double mixture_mean = mean_sum / weight_sum;
    long double variance_sum = 0;
    for (int s = 0; s < S; ++s) {
      row_work.add();
      const int cell = s + S * k;
      if (weights_p[cell] == 0) {
        continue;
      }
      const long double component_mean =
        moment_mean_work[static_cast<std::size_t>(s)];
      const long double centered = component_mean - mixture_mean;
      const long double component_variance =
        moment_variance_work[static_cast<std::size_t>(s)];
      variance_sum += static_cast<long double>(weights_p[cell]) *
        (component_variance + centered * centered);
    }
    const long double mixture_variance = variance_sum / weight_sum;
    if (!std::isfinite(mixture_variance) || mixture_variance < 0) {
      status[static_cast<std::size_t>(k)] = summary_invalid_variance;
      return;
    }

    const double log_lower_sum = lower_sum.value();
    const double log_upper_sum = upper_sum.value();
    if (tails_unavailable) {
      log_lower_out_p[k] = NA_REAL;
      log_upper_out_p[k] = NA_REAL;
    } else if (log_lower_sum == -INFINITY && log_upper_sum == -INFINITY) {
      status[static_cast<std::size_t>(k)] = summary_undefined_tails;
      return;
    } else if (log_lower_sum == -INFINITY) {
      log_lower_out_p[k] = -INFINITY;
      log_upper_out_p[k] = 0;
    } else if (log_upper_sum == -INFINITY) {
      log_lower_out_p[k] = 0;
      log_upper_out_p[k] = -INFINITY;
    } else {
      const double difference = log_upper_sum - log_lower_sum;
      log_lower_out_p[k] = difference > 0 ?
        -difference - std::log1p(std::exp(-difference)) :
        -std::log1p(std::exp(difference));
      log_upper_out_p[k] = difference < 0 ?
        difference - std::log1p(std::exp(difference)) :
        -std::log1p(std::exp(-difference));
    }
    mean_out_p[k]      = static_cast<double>(mixture_mean);
    variance_out_p[k]  = static_cast<double>(mixture_variance);
  };

  robma_parallel::for_rows(K, schedule, summarize);

  for (int k = 0; k < K; ++k) {
    switch (status[static_cast<std::size_t>(k)]) {
    case summary_invalid_weight:
      Rf_error("'psis_weights' must be finite and non-negative.");
    case summary_unavailable_moments:
      Rf_error("%s", selnorm_unrepresentable_moments_message);
    case summary_invalid_moments:
      Rf_error("Selected-normal predictive moments are not finite.");
    case summary_invalid_weight_sum:
      Rf_error("'psis_weights' must have positive finite column sums.");
    case summary_invalid_variance:
      Rf_error(
        "Selected-normal predictive variance is not finite and non-negative."
      );
    case summary_undefined_tails:
      Rf_error("Selected-normal predictive tail probabilities are undefined.");
    default:
      break;
    }
  }

  UNPROTECT(6);
  return out;
  ROBMA_NATIVE_END
}
