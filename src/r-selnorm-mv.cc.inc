// A posterior branch can use only a subset of the shared prior-union cuts.
// This call-owned view removes exact adjacent repeats after validation, without
// changing the context, branch routing, or the weight at an observed boundary.
struct CompactStepPartition {
  std::vector<double> omega, lower, upper;
  std::vector<int> obs_bin, bin_map;

  bool prepare(const double *weights, const double *z_lower,
      const double *z_upper, int bins, const int *observed = nullptr,
      int observations = 0)
  {
    if (bins <= 2) return false;
    bool repeated = false;
    for (int bin = 1; bin < bins; ++bin)
      repeated = repeated || weights[bin] == weights[bin - 1];
    if (!repeated) return false;
    omega.clear(); lower.clear(); upper.clear();
    if (observed != nullptr) bin_map.resize(bins);
    for (int first = 0; first < bins;) {
      int last = first;
      while (last + 1 < bins && weights[last + 1] == weights[first]) ++last;
      omega.push_back(weights[first]);
      lower.push_back(z_lower[last]);
      upper.push_back(z_upper[first]);
      if (observed != nullptr)
        for (int bin = first; bin <= last; ++bin)
          bin_map[bin] = static_cast<int>(omega.size());
      first = last + 1;
    }
    if (observed != nullptr) {
      obs_bin.resize(observations);
      for (int row = 0; row < observations; ++row)
        obs_bin[row] = bin_map[observed[row] - 1];
    }
    return true;
  }
};

static SelNormDenseIntegration validate_dense_integration(
  SEXP normalizer_quadrature, SEXP relative_tolerance,
  std::vector<double> &normalizer_orders)
{
  check_real(relative_tolerance, "relative_tolerance");
  if (Rf_length(relative_tolerance) != 1 ||
      !std::isfinite(REAL(relative_tolerance)[0]) ||
      REAL(relative_tolerance)[0] <= 0.0) {
    Rf_error("'relative_tolerance' must be one finite positive value.");
  }

  SelNormDenseIntegration integration_rules;
  if (normalizer_quadrature != R_NilValue) {
    if (TYPEOF(normalizer_quadrature) != VECSXP || Rf_length(normalizer_quadrature) != 3) {
      Rf_error("'normalizer_quadrature' must contain orders, nodes, and log weights.");
    }
    SEXP orders = VECTOR_ELT(normalizer_quadrature, 0);
    SEXP nodes = VECTOR_ELT(normalizer_quadrature, 1);
    SEXP log_weights = VECTOR_ELT(normalizer_quadrature, 2);
    if (TYPEOF(orders) != INTSXP && TYPEOF(orders) != REALSXP) {
      Rf_error("'normalizer_quadrature' orders must be numeric.");
    }
    check_real(nodes, "normalizer_quadrature nodes");
    check_real(log_weights, "normalizer_quadrature log weights");
    if (Rf_length(orders) < 3 || Rf_length(nodes) != Rf_length(log_weights)) {
      Rf_error("'normalizer_quadrature' dimensions are invalid.");
    }
  }
  if (normalizer_quadrature != R_NilValue) {
    SEXP orders = VECTOR_ELT(normalizer_quadrature, 0);
    SEXP nodes = VECTOR_ELT(normalizer_quadrature, 1);
    SEXP log_weights = VECTOR_ELT(normalizer_quadrature, 2);
    int previous = 0;
    int total = 0;
    for (int i = 0; i < Rf_length(orders); ++i) {
      const double raw = TYPEOF(orders) == INTSXP ? INTEGER(orders)[i] : REAL(orders)[i];
      if (!std::isfinite(raw) || raw != std::floor(raw) || raw <= previous ||
          raw > Rf_length(nodes) - total) {
        Rf_error("'normalizer_quadrature' orders are invalid.");
      }
      previous = static_cast<int>(raw);
      total += previous;
      normalizer_orders.push_back(raw);
    }
    if (total != Rf_length(nodes)) Rf_error("'normalizer_quadrature' orders do not match its nodes.");
    for (int i = 0; i < total; ++i) {
      if (!std::isfinite(REAL(nodes)[i]) || std::isnan(REAL(log_weights)[i]) ||
          REAL(log_weights)[i] == R_PosInf) {
        Rf_error("'normalizer_quadrature' nodes or log weights are invalid.");
      }
    }
    integration_rules.nodes = REAL(nodes);
    integration_rules.log_weights = REAL(log_weights);
    integration_rules.orders = normalizer_orders.data();
    integration_rules.rule_count = normalizer_orders.size();
    integration_rules.tolerance = REAL(relative_tolerance)[0];
    SEXP bounded_cdf = Rf_getAttrib(normalizer_quadrature, Rf_install("bounded_cdf"));
    if (bounded_cdf != R_NilValue) {
      if (TYPEOF(bounded_cdf) != LGLSXP || Rf_length(bounded_cdf) != 1 ||
          LOGICAL(bounded_cdf)[0] == NA_LOGICAL) {
        Rf_error("'bounded_cdf' must be one nonmissing logical value.");
      }
      integration_rules.bounded_cdf = LOGICAL(bounded_cdf)[0] == TRUE;
    }
  }

  return integration_rules;
}

static int validate_vector_rule(SEXP vector_rule, int samples)
{
  check_integer(vector_rule, "vector_rule");
  const int length = Rf_length(vector_rule);
  if (length != 1 && length != samples) {
    Rf_error("'vector_rule' must have length 1 or one value per posterior sample.");
  }
  for (int i = 0; i < length; ++i) {
    if (INTEGER(vector_rule)[i] < SELVECTOR_PRODUCT ||
        INTEGER(vector_rule)[i] > SELVECTOR_BEST_TWO_SIDED) {
      Rf_error("'vector_rule' must identify product, one-sided best, or two-sided best weights.");
    }
  }
  return length;
}

extern "C" SEXP RoBMA_selnorm_gaussian_event_mass_batch(
  SEXP means, SEXP covariance_lower, SEXP sei, SEXP omega,
  SEXP z_lower, SEXP z_upper, SEXP sign, SEXP kernel_mode,
  SEXP vector_rule, SEXP lower, SEXP upper, SEXP qmc,
  SEXP points, SEXP scrambles, SEXP rank_one_loading,
  SEXP normalizer_quadrature, SEXP relative_tolerance)
{
  check_real(means, "means");
  if (rank_one_loading == R_NilValue) {
    check_real(covariance_lower, "covariance_lower");
  } else {
    check_real(rank_one_loading, "rank_one_loading");
    if (covariance_lower != R_NilValue) {
      Rf_error("'rank_one_loading' and 'covariance_lower' are alternative Gaussian covariance representations.");
    }
  }
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(kernel_mode, "kernel_mode");
  check_real(lower, "lower");
  check_real(upper, "upper");
  check_real(qmc, "qmc");
  check_integer(points, "points");
  check_integer(scrambles, "scrambles");
  const int S = matrix_nrow(means, "means");
  const int K = matrix_ncol(means, "means");
  const int B = matrix_ncol(omega, "omega");
  const int L = K * (K + 1) / 2;
  const int rule_length = validate_vector_rule(vector_rule, S);
  const int mode_length = Rf_length(kernel_mode);
  if (S < 1 || K < 1 || B < 1 || Rf_length(sei) != K ||
      Rf_length(z_lower) != B || Rf_length(z_upper) != B ||
      (mode_length != 1 && mode_length != S)) {
    Rf_error("Gaussian selection event inputs have invalid dimensions.");
  }
  if (rank_one_loading == R_NilValue) {
    validate_matrix_dims_named(covariance_lower, "covariance_lower", S, L, "means");
  } else {
    validate_matrix_dims_named(rank_one_loading, "rank_one_loading", S, K, "means");
    for (R_xlen_t i = 0; i < XLENGTH(rank_one_loading); ++i) {
      if (!std::isfinite(REAL(rank_one_loading)[i])) {
        Rf_error("'rank_one_loading' must contain finite values.");
      }
    }
  }
  validate_matrix_dims_named(omega, "omega", S, B, "means");
  validate_matrix_dims_named(lower, "lower", S, K, "means");
  validate_matrix_dims_named(upper, "upper", S, K, "means");
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1) ||
      Rf_length(points) != 1 || INTEGER(points)[0] < 1 ||
      Rf_length(scrambles) != 1 || INTEGER(scrambles)[0] < 2) {
    Rf_error("Gaussian selection event controls are invalid.");
  }
  if (!selnorm_is_descending_step_partition(REAL(z_lower), REAL(z_upper), B)) {
    Rf_error("Selection bounds do not define a descending step partition.");
  }
  validate_omega_matrix(REAL(omega), S, B);
  const int n_points = INTEGER(points)[0];
  const int n_scrambles = INTEGER(scrambles)[0];
  if (rank_one_loading == R_NilValue && K > 1 &&
      XLENGTH(qmc) != static_cast<R_xlen_t>(2 * K) * n_points * n_scrambles) {
    Rf_error("'qmc' dimensions do not match the integration settings.");
  }
  for (R_xlen_t i = 0; i < XLENGTH(qmc); ++i) {
    if (!(REAL(qmc)[i] > 0.0 && REAL(qmc)[i] < 1.0)) {
      Rf_error("'qmc' must contain finite values strictly between zero and one.");
    }
  }
  for (int i = 0; i < K; ++i) {
    if (!(REAL(sei)[i] > 0.0) || !std::isfinite(REAL(sei)[i])) {
      Rf_error("'sei' must contain finite positive values.");
    }
  }
  for (R_xlen_t i = 0; i < XLENGTH(means); ++i) {
    if (!std::isfinite(REAL(means)[i]) || std::isnan(REAL(lower)[i]) ||
        std::isnan(REAL(upper)[i])) {
      Rf_error("Gaussian selection event locations or limits are invalid.");
    }
  }
  std::vector<double> normalizer_orders;
  SelNormDenseIntegration integration_rules;
  if (normalizer_quadrature != R_NilValue) {
    integration_rules = validate_dense_integration(
      normalizer_quadrature, relative_tolerance, normalizer_orders);
  }
  SEXP quadrature_error = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP log_mass = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP mcse = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP cdf_error = PROTECT(integration_rules.bounded_cdf ? Rf_allocVector(REALSXP, S) : R_NilValue);
  SEXP cdf_log_error = PROTECT(integration_rules.bounded_cdf ? Rf_allocVector(REALSXP, S) : R_NilValue);
  // Kernel modes are validated before the row loop so parallel workers never
  // raise R errors from a worker thread.
  for (int s = 0; s < S; ++s) {
    const int mode = scalar_or_row_int(INTEGER(kernel_mode), mode_length, s, "kernel_mode");
    if (mode != SELKERNEL_NORMAL && mode != SELKERNEL_STEP) {
      Rf_error("Gaussian selection events require a normal or step selection kernel.");
    }
  }
  const int threads = robma_batch_threads(S);
  robma_parallel::for_rows(S, threads, [&](int s) {
    std::vector<double> row_mean(static_cast<std::size_t>(K));
    std::vector<double> row_covariance(static_cast<std::size_t>(L));
    std::vector<double> row_weights(static_cast<std::size_t>(B));
    std::vector<double> row_lo(static_cast<std::size_t>(K));
    std::vector<double> row_hi(static_cast<std::size_t>(K));
    std::vector<double> row_loading(static_cast<std::size_t>(K));
    const int mode = scalar_or_row_int(INTEGER(kernel_mode), mode_length, s, "kernel_mode");
    for (int i = 0; i < K; ++i) {
      row_mean[i] = REAL(means)[s + S * i];
      row_lo[i] = REAL(lower)[s + S * i];
      row_hi[i] = REAL(upper)[s + S * i];
    }
    if (rank_one_loading == R_NilValue) {
      for (int i = 0; i < L; ++i) row_covariance[i] = REAL(covariance_lower)[s + S * i];
    } else {
      for (int i = 0; i < K; ++i) row_loading[i] = REAL(rank_one_loading)[s + S * i];
    }
    for (int i = 0; i < B; ++i) row_weights[i] = REAL(omega)[s + S * i];
    SelNormDenseIntegration integration = integration_rules;
    REAL(log_mass)[s] = cpp_selnorm_gaussian_event_log_mass(row_mean.data(),
      row_covariance.data(), K, REAL(sei), row_weights.data(), B, REAL(z_lower),
      REAL(z_upper), INTEGER(sign)[0], mode,
      scalar_or_row_int(INTEGER(vector_rule), rule_length, s, "vector_rule"),
      row_lo.data(), row_hi.data(), REAL(qmc), n_points, n_scrambles, REAL(mcse) + s,
      rank_one_loading == R_NilValue ? nullptr : row_loading.data(),
      normalizer_quadrature == R_NilValue ? nullptr : &integration);
    // This combines covariance/exterior bounds with estimated interior error;
    // it is not a Monte Carlo standard error or a rigorous total error bound.
    REAL(quadrature_error)[s] = integration.used_envelope ?
      integration.covariance_width + 2.0 * integration.quadrature_change +
      integration.tail_bound + integration.cdf_relative_error : 0.0;
    if (integration_rules.bounded_cdf) {
      REAL(cdf_error)[s] = integration.cdf_relative_error;
      REAL(cdf_log_error)[s] = integration.cdf_log_error;
    }
  });
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_VECTOR_ELT(out, 0, log_mass);
  SET_VECTOR_ELT(out, 1, mcse);
  SET_VECTOR_ELT(out, 2, quadrature_error);
  SET_STRING_ELT(names, 0, Rf_mkChar("log_mass"));
  SET_STRING_ELT(names, 1, Rf_mkChar("relative_mcse"));
  SET_STRING_ELT(names, 2, Rf_mkChar("relative_quadrature_error"));
  if (integration_rules.bounded_cdf) {
    Rf_setAttrib(quadrature_error, Rf_install("cdf_relative_error"), cdf_error);
    Rf_setAttrib(quadrature_error, Rf_install("cdf_log_error"), cdf_log_error);
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(7);
  return out;
}

static SEXP selnorm_mnorm_step_batch(
  SEXP yi, SEXP means, SEXP covariance_lower, SEXP sei, SEXP omega,
  SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP telescope_probabilities, SEXP kernel_mode, SEXP qmc,
  SEXP points, SEXP scrambles, SEXP relative_tolerance,
  SEXP z, SEXP probability, SEXP factors, bool return_normalizer,
  SEXP vector_rule, SEXP normalizer_quadrature)
{
  check_real(yi, "yi");
  check_real(means, "means");
  check_real(covariance_lower, "covariance_lower");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(obs_bin, "obs_bin");
  check_integer(sign, "sign");
  check_logical(telescope_probabilities, "telescope_probabilities");
  check_integer(kernel_mode, "kernel_mode");
  check_real(qmc, "qmc");
  check_integer(points, "points");
  check_integer(scrambles, "scrambles");

  const int S = matrix_nrow(means, "means");
  const int vector_length = validate_vector_rule(vector_rule, S);
  const int K = matrix_ncol(means, "means");
  const int L = K * (K + 1) / 2;
  const int B = matrix_ncol(omega, "omega");
  if (S < 1 || K < 1) {
    Rf_error("'means' must have positive dimensions.");
  }
  int rank = 0;
  if (factors != R_NilValue) {
    if (TYPEOF(factors) != VECSXP ||
        (Rf_length(factors) != 2 && Rf_length(factors) != 5)) {
      Rf_error("'factors' must contain residual SDs and loadings.");
    }
    check_real(VECTOR_ELT(factors, 0), "residual_sd");
    check_real(VECTOR_ELT(factors, 1), "loading");
    validate_matrix_dims_named(VECTOR_ELT(factors, 0), "residual_sd", S, K, "means");
    const int columns = matrix_ncol(VECTOR_ELT(factors, 1), "loading");
    if (columns % K != 0 || columns == 0 || columns / K > K ||
        matrix_nrow(VECTOR_ELT(factors, 1), "loading") != S) {
      Rf_error("'loading' has invalid dimensions.");
    }
    rank = columns / K;
  }
  validate_matrix_dims_named(
    covariance_lower, "covariance_lower", S, L, "means"
  );
  if (matrix_nrow(omega, "omega") != S || B < 1) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(yi) != K || Rf_length(sei) != K ||
      Rf_length(obs_bin) != K) {
    Rf_error("'yi', 'sei', and 'obs_bin' must have one value per observation.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("Selection bounds must have one value per selection bin.");
  }
  if (Rf_length(sign) != 1 ||
      (INTEGER(sign)[0] != -1 && INTEGER(sign)[0] != 1)) {
    Rf_error("'sign' must be either -1 or 1.");
  }
  const int kernel_length = Rf_length(kernel_mode);
  if (kernel_length != 1 && kernel_length != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(points) != 1 || INTEGER(points)[0] < 1 ||
      Rf_length(scrambles) != 1 || INTEGER(scrambles)[0] < 2) {
    Rf_error("'points' and 'scrambles' contain invalid integration settings.");
  }
  const int n_points    = INTEGER(points)[0];
  const int n_scrambles = INTEGER(scrambles)[0];
  const bool quadrature = factors != R_NilValue && Rf_length(factors) == 5;
  const int *loading_support = nullptr;
  if (quadrature) {
    check_real(VECTOR_ELT(factors, 2), "nodes");
    check_real(VECTOR_ELT(factors, 3), "log_weights");
    const int order = Rf_length(VECTOR_ELT(factors, 2));
    double expected_points = 1.0;
    for (int factor = 0; factor < rank; ++factor) expected_points *= order;
    if (order < 1 || rank > 4 || Rf_length(VECTOR_ELT(factors, 3)) != order ||
        n_points != expected_points) {
      Rf_error("'factors' quadrature dimensions are invalid.");
    }
    SEXP support = VECTOR_ELT(factors, 4);
    if (rank > 2 || support != R_NilValue) {
      if (TYPEOF(support) != LGLSXP) {
        Rf_error("'loading_support' must be a logical matrix.");
      }
      if (matrix_nrow(support, "loading_support") != K ||
          matrix_ncol(support, "loading_support") != rank) {
        Rf_error("'loading_support' has invalid dimensions.");
      }
      loading_support = LOGICAL(support);
      for (int index = 0; index < K * rank; ++index) {
        if (loading_support[index] == NA_LOGICAL) {
          Rf_error("'loading_support' must not contain missing values.");
        }
      }
    }
  }
  const R_xlen_t expected_qmc = static_cast<R_xlen_t>(n_scrambles) *
    static_cast<R_xlen_t>(n_points) * static_cast<R_xlen_t>(2 * K);
  if (!quadrature && K > 1 && XLENGTH(qmc) != expected_qmc) {
    Rf_error("'qmc' dimensions do not match the integration settings.");
  }
  std::vector<double> normalizer_orders;
  SelNormDenseIntegration integration_rules = validate_dense_integration(
    normalizer_quadrature, relative_tolerance, normalizer_orders);
  // The dense R interface approximates only post-fit log densities.
  // Its existing projected-outcome path retains the original primitive.
  if (z != R_NilValue) integration_rules.bounded_cdf = false;
  if (integration_rules.bounded_cdf &&
      !selnorm_is_descending_step_partition(REAL(z_lower), REAL(z_upper), B)) {
    Rf_error("Selection bounds do not define a descending step partition.");
  }

  validate_omega_matrix(REAL(omega), S, B);
  for (int k = 0; k < K; ++k) {
    if (!std::isfinite(REAL(yi)[k]) || !std::isfinite(REAL(sei)[k]) ||
        REAL(sei)[k] <= 0.0 || INTEGER(obs_bin)[k] < 1 ||
        INTEGER(obs_bin)[k] > B) {
      Rf_error("Observed selection inputs are invalid.");
    }
  }

  SEXP log_density = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP relative_mcse = PROTECT(Rf_allocVector(REALSXP, S));
  const bool project = z != R_NilValue;
  SEXP log_normalizer = PROTECT(Rf_allocVector(REALSXP, project ? 0 : S));
  std::fill_n(REAL(log_normalizer), project ? 0 : S, NA_REAL);
  const int Z = project ? Rf_length(z) : 0;
  if (project) {
    check_real(z, "z");
    check_logical(probability, "probability");
    if (Rf_length(probability) != 1 || LOGICAL(probability)[0] == NA_LOGICAL) {
      Rf_error("'probability' must be one logical value.");
    }
    for (int j = 0; j < Z; ++j) {
      if (!std::isfinite(REAL(z)[j])) Rf_error("'z' must be finite.");
      if (LOGICAL(probability)[0] && REAL(z)[j] < 0) {
        Rf_error("Z thresholds must be non-negative.");
      }
    }
  }
  SEXP density = PROTECT(Rf_allocMatrix(REALSXP, S, Z));
  SEXP integration_diagnostics = PROTECT(project ? R_NilValue : Rf_allocMatrix(REALSXP, S, 4));
  SEXP cdf_error = PROTECT(integration_rules.bounded_cdf ? Rf_allocVector(REALSXP, S) : R_NilValue);
  SEXP cdf_log_error = PROTECT(integration_rules.bounded_cdf ? Rf_allocVector(REALSXP, S) : R_NilValue);
  if (!project) {
    SEXP diagnostic_names = PROTECT(Rf_allocVector(STRSXP, 4));
    SEXP diagnostic_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_STRING_ELT(diagnostic_names, 0, Rf_mkChar("used_covariance_envelope"));
    SET_STRING_ELT(diagnostic_names, 1, Rf_mkChar("quadrature_change"));
    SET_STRING_ELT(diagnostic_names, 2, Rf_mkChar("covariance_width"));
    SET_STRING_ELT(diagnostic_names, 3, Rf_mkChar("tail_bound"));
    SET_VECTOR_ELT(diagnostic_dimnames, 1, diagnostic_names);
    Rf_setAttrib(integration_diagnostics, R_DimNamesSymbol, diagnostic_dimnames);
    UNPROTECT(2);
  }
  std::vector<double> projected(static_cast<std::size_t>(Z));
  std::vector<double> mean(static_cast<std::size_t>(K));
  std::vector<double> lower(static_cast<std::size_t>(L));
  std::vector<double> omega_row(static_cast<std::size_t>(B));
  std::vector<double> residual_sd(static_cast<std::size_t>(K));
  std::vector<double> loading(static_cast<std::size_t>(K * rank));
  CompactStepPartition partition;

  // Per-row factor inputs are validated before the row loop so the parallel
  // workers below never need to raise R errors from a worker thread.
  if (rank) {
    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < K; ++i) {
        const double value = REAL(VECTOR_ELT(factors, 0))[s + S * i];
        if (!(value > 0) || !std::isfinite(value)) {
          Rf_error("'residual_sd' must be finite and positive.");
        }
      }
      for (int j = 0; j < K * rank; ++j) {
        const double value = REAL(VECTOR_ELT(factors, 1))[s + S * j];
        if (!std::isfinite(value)) Rf_error("'loading' must be finite.");
        if (loading_support != nullptr && !loading_support[j] && value != 0.0) {
          Rf_error("'loading' contains nonzero values outside 'loading_support'.");
        }
      }
    }
  }

  // The existing owned row buffers are also the previous physical key.
  // Only this call's non-projected dense likelihood can reuse a completed row;
  // static observations, selection controls and quadrature/QMC are call-common.
  const int threads = robma_batch_threads(S);
  const bool adjacent_dense = !project && factors == R_NilValue && threads == 1;
  bool previous_accepted = false;
  int previous_kernel = 0, previous_rule = 0;
  const auto same_double = [](const double &left, const double &right) {
    return std::memcmp(&left, &right, sizeof(double)) == 0;
  };
  const auto compute_row = [&](int s, std::vector<double> &row_mean,
      std::vector<double> &row_lower, std::vector<double> &row_omega,
      std::vector<double> &row_residual_sd, std::vector<double> &row_loading,
      std::vector<double> &row_projected, CompactStepPartition &row_partition) {
    const int row_kernel = scalar_or_row_int(INTEGER(kernel_mode), kernel_length, s, "kernel_mode");
    const int row_rule = scalar_or_row_int(INTEGER(vector_rule), vector_length, s, "vector_rule");
    for (int k = 0; k < K; ++k) {
      row_mean[static_cast<std::size_t>(k)] = REAL(means)[s + S * k];
    }
    for (int l = 0; l < L; ++l) {
      row_lower[static_cast<std::size_t>(l)] = REAL(covariance_lower)[s + S * l];
    }
    for (int b = 0; b < B; ++b) {
      row_omega[static_cast<std::size_t>(b)] = REAL(omega)[s + S * b];
    }
    if (rank) {
      for (int i = 0; i < K; ++i) {
        row_residual_sd[i] = REAL(VECTOR_ELT(factors, 0))[s + S * i];
      }
      for (int j = 0; j < K * rank; ++j) {
        row_loading[j] = REAL(VECTOR_ELT(factors, 1))[s + S * j];
      }
    }

    double error = 0.0;
    SelNormZProjection projection;
    projection.z = project ? REAL(z) : nullptr;
    projection.size = Z;
    projection.density = row_projected.data();
    projection.relative_error = R_PosInf;
    projection.probability = project && LOGICAL(probability)[0];
    projection.rank = rank;
    projection.factor_loading = rank ? row_loading.data() : nullptr;
    projection.loading_support = loading_support;
    projection.residual_sd = rank ? row_residual_sd.data() : nullptr;
    projection.order = quadrature ? Rf_length(VECTOR_ELT(factors, 2)) : 0;
    projection.nodes = quadrature ? REAL(VECTOR_ELT(factors, 2)) : nullptr;
    projection.log_weights = quadrature ? REAL(VECTOR_ELT(factors, 3)) : nullptr;
    SelNormDenseIntegration integration = integration_rules;
    const bool compact = !project && integration_rules.bounded_cdf &&
      row_kernel == SELKERNEL_STEP && row_rule == SELVECTOR_PRODUCT &&
      LOGICAL(telescope_probabilities)[0] && row_partition.prepare(row_omega.data(),
        REAL(z_lower), REAL(z_upper), B, INTEGER(obs_bin), K);
    REAL(log_density)[s] = cpp_selnorm_mnorm_step_lpdf(
      REAL(yi), row_mean.data(), row_lower.data(), K, REAL(sei),
      compact ? row_partition.omega.data() : row_omega.data(),
      compact ? static_cast<int>(row_partition.omega.size()) : B,
      compact ? row_partition.lower.data() : REAL(z_lower),
      compact ? row_partition.upper.data() : REAL(z_upper),
      compact ? row_partition.obs_bin.data() : INTEGER(obs_bin), INTEGER(sign)[0],
      LOGICAL(telescope_probabilities)[0],
      row_kernel,
      REAL(qmc), n_points, n_scrambles, &error,
      project ? &projection : nullptr,
      return_normalizer ? REAL(log_normalizer) + s : nullptr,
      row_rule,
      normalizer_quadrature != R_NilValue ? &integration : nullptr
    );
    REAL(relative_mcse)[s] = project ? std::max(error, projection.relative_error) : error;
    if (!project) {
      REAL(integration_diagnostics)[s] = integration.used_envelope;
      REAL(integration_diagnostics)[s + S] = integration.quadrature_change;
      REAL(integration_diagnostics)[s + 2 * S] = integration.covariance_width;
      REAL(integration_diagnostics)[s + 3 * S] = integration.tail_bound;
      if (integration_rules.bounded_cdf) {
        REAL(cdf_error)[s] = integration.cdf_relative_error;
        REAL(cdf_log_error)[s] = integration.cdf_log_error;
      }
    }
    for (int j = 0; j < Z; ++j) REAL(density)[s + S * j] = row_projected[j];
    return std::pair<int, SelNormDenseIntegration>(row_kernel, integration);
  };
  if (threads == 1) {
    for (int s = 0; s < S; ++s) {
      const int row_kernel = scalar_or_row_int(INTEGER(kernel_mode), kernel_length, s, "kernel_mode");
      const int row_rule = scalar_or_row_int(INTEGER(vector_rule), vector_length, s, "vector_rule");
      bool same_row = adjacent_dense && previous_accepted &&
        row_kernel == previous_kernel && row_rule == previous_rule;
      for (int k = 0; k < K; ++k) {
        same_row = same_row && same_double(mean[static_cast<std::size_t>(k)], REAL(means)[s + S * k]);
        mean[static_cast<std::size_t>(k)] = REAL(means)[s + S * k];
      }
      for (int l = 0; l < L; ++l) {
        same_row = same_row && same_double(lower[static_cast<std::size_t>(l)], REAL(covariance_lower)[s + S * l]);
        lower[static_cast<std::size_t>(l)] =
          REAL(covariance_lower)[s + S * l];
      }
      for (int b = 0; b < B; ++b) {
        same_row = same_row && same_double(omega_row[static_cast<std::size_t>(b)], REAL(omega)[s + S * b]);
        omega_row[static_cast<std::size_t>(b)] = REAL(omega)[s + S * b];
      }
      if (same_row) {
        REAL(log_density)[s] = REAL(log_density)[s - 1];
        REAL(relative_mcse)[s] = REAL(relative_mcse)[s - 1];
        REAL(log_normalizer)[s] = REAL(log_normalizer)[s - 1];
        for (int column = 0; column < 4; ++column)
          REAL(integration_diagnostics)[s + S * column] = REAL(integration_diagnostics)[s - 1 + S * column];
        if (integration_rules.bounded_cdf) {
          REAL(cdf_error)[s] = REAL(cdf_error)[s - 1];
          REAL(cdf_log_error)[s] = REAL(cdf_log_error)[s - 1];
        }
        R_CheckUserInterrupt();
        continue;
      }

      auto row_result = compute_row(s, mean, lower, omega_row, residual_sd,
                                    loading, projected, partition);
      previous_kernel = row_result.first;
      previous_rule = scalar_or_row_int(INTEGER(vector_rule), vector_length, s, "vector_rule");
      const SelNormDenseIntegration &integration = row_result.second;
      previous_accepted = adjacent_dense && std::isfinite(REAL(log_density)[s]) &&
        std::isfinite(REAL(relative_mcse)[s]) && REAL(relative_mcse)[s] >= 0.0 &&
        REAL(relative_mcse)[s] <= REAL(relative_tolerance)[0] &&
        (!return_normalizer || std::isfinite(REAL(log_normalizer)[s]));
      if (previous_accepted) {
        for (int column = 0; column < 4; ++column) {
          const double diagnostic = REAL(integration_diagnostics)[s + S * column];
          if (!std::isfinite(diagnostic) || diagnostic < 0.0) previous_accepted = false;
        }
        if (integration_rules.bounded_cdf &&
            (!std::isfinite(integration.cdf_relative_error) || integration.cdf_relative_error < 0.0 ||
             !std::isfinite(integration.cdf_log_error) || integration.cdf_log_error < 0.0)) previous_accepted = false;
      }
      R_CheckUserInterrupt();
    }
  } else {
    // Rows are independent: each writes only its own output slots. Workers
    // own their buffers and partitions, so results match the serial rows.
    robma_parallel::for_rows(S, threads, [&](int s) {
      std::vector<double> row_mean(static_cast<std::size_t>(K));
      std::vector<double> row_lower(static_cast<std::size_t>(L));
      std::vector<double> row_omega(static_cast<std::size_t>(B));
      std::vector<double> row_residual_sd(static_cast<std::size_t>(K));
      std::vector<double> row_loading(static_cast<std::size_t>(K * rank));
      std::vector<double> row_projected(static_cast<std::size_t>(Z));
      CompactStepPartition row_partition;
      compute_row(s, row_mean, row_lower, row_omega, row_residual_sd,
                  row_loading, row_projected, row_partition);
    });
  }

  SEXP output = PROTECT(Rf_allocVector(VECSXP, project ? 3 : 4));
  SEXP names  = PROTECT(Rf_allocVector(STRSXP, project ? 3 : 4));
  SET_VECTOR_ELT(output, 0, log_density);
  SET_VECTOR_ELT(output, 1, relative_mcse);
  SET_STRING_ELT(names, 0, Rf_mkChar("log_density"));
  SET_STRING_ELT(names, 1, Rf_mkChar("relative_mcse"));
  if (project) {
    SET_VECTOR_ELT(output, 2, density);
    SET_STRING_ELT(names, 2, Rf_mkChar("density"));
  } else {
    SET_VECTOR_ELT(output, 2, log_normalizer);
    SET_STRING_ELT(names, 2, Rf_mkChar("log_normalizer"));
    SET_VECTOR_ELT(output, 3, integration_diagnostics);
    SET_STRING_ELT(names, 3, Rf_mkChar("integration_diagnostics"));
    if (integration_rules.bounded_cdf) {
      Rf_setAttrib(integration_diagnostics, Rf_install("cdf_relative_error"), cdf_error);
      Rf_setAttrib(integration_diagnostics, Rf_install("cdf_log_error"), cdf_log_error);
    }
  }
  Rf_setAttrib(output, R_NamesSymbol, names);

  UNPROTECT(9);
  return output;
}

extern "C" SEXP RoBMA_selnorm_mnorm_step_loglik_batch(
  SEXP yi, SEXP means, SEXP covariance_lower, SEXP sei, SEXP omega,
  SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP telescope_probabilities, SEXP kernel_mode, SEXP qmc,
  SEXP points, SEXP scrambles, SEXP relative_tolerance,
  SEXP return_normalizer, SEXP vector_rule, SEXP normalizer_quadrature)
{
  check_logical(return_normalizer, "return_normalizer");
  if (LOGICAL(return_normalizer)[0] == NA_LOGICAL) {
    Rf_error("'return_normalizer' must not be missing.");
  }
  return selnorm_mnorm_step_batch(yi, means, covariance_lower, sei, omega,
    z_lower, z_upper, obs_bin, sign, telescope_probabilities, kernel_mode,
    qmc, points, scrambles, relative_tolerance, R_NilValue, R_NilValue,
    R_NilValue, LOGICAL(return_normalizer)[0], vector_rule, normalizer_quadrature);
}

extern "C" SEXP RoBMA_selnorm_mnorm_zplot_batch(
  SEXP yi, SEXP means, SEXP covariance_lower, SEXP sei, SEXP omega,
  SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP telescope_probabilities, SEXP kernel_mode, SEXP qmc,
  SEXP points, SEXP scrambles, SEXP relative_tolerance,
  SEXP z, SEXP probability, SEXP factors, SEXP vector_rule)
{
  return selnorm_mnorm_step_batch(yi, means, covariance_lower, sei, omega,
    z_lower, z_upper, obs_bin, sign, telescope_probabilities, kernel_mode,
    qmc, points, scrambles, relative_tolerance, z, probability, factors, false,
    vector_rule, R_NilValue);
}

extern "C" SEXP RoBMA_selnorm_cluster_step_loglik_batch(
  SEXP yi, SEXP means, SEXP residual_sd, SEXP loading, SEXP sei, SEXP omega,
  SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP telescope_probabilities, SEXP kernel_mode, SEXP quadrature_nodes,
  SEXP quadrature_log_weights, SEXP quadrature_orders, SEXP qmc,
  SEXP initial_points, SEXP max_points, SEXP scrambles,
  SEXP relative_tolerance, SEXP return_normalizer, SEXP vector_rule)
{
  check_logical(return_normalizer, "return_normalizer");
  if (LOGICAL(return_normalizer)[0] == NA_LOGICAL) {
    Rf_error("'return_normalizer' must not be missing.");
  }
  check_real(yi, "yi");
  check_real(means, "means");
  check_real(residual_sd, "residual_sd");
  check_real(loading, "loading");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(obs_bin, "obs_bin");
  check_integer(sign, "sign");
  check_logical(telescope_probabilities, "telescope_probabilities");
  check_integer(kernel_mode, "kernel_mode");
  check_real(quadrature_nodes, "quadrature_nodes");
  check_real(quadrature_log_weights, "quadrature_log_weights");
  check_real(quadrature_orders, "quadrature_orders");
  check_real(qmc, "qmc");
  check_integer(initial_points, "initial_points");
  check_integer(max_points, "max_points");
  check_integer(scrambles, "scrambles");
  check_real(relative_tolerance, "relative_tolerance");

  const int S = matrix_nrow(means, "means");
  const int vector_length = validate_vector_rule(vector_rule, S);
  const int K = matrix_ncol(means, "means");
  const int B = matrix_ncol(omega, "omega");
  validate_matrix_dims_named(residual_sd, "residual_sd", S, K, "means");
  validate_matrix_dims_named(loading, "loading", S, K, "means");
  if (S < 1 || K < 1 || B < 1 || matrix_nrow(omega, "omega") != S ||
      Rf_length(yi) != K || Rf_length(sei) != K || Rf_length(obs_bin) != K ||
      Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("Cluster selection inputs have invalid dimensions.");
  }
  if (Rf_length(sign) != 1 ||
      (INTEGER(sign)[0] != -1 && INTEGER(sign)[0] != 1) ||
      Rf_length(telescope_probabilities) != 1 ||
      Rf_length(quadrature_nodes) < 1 ||
      Rf_length(quadrature_log_weights) != Rf_length(quadrature_nodes) ||
      Rf_length(quadrature_orders) < 2 ||
      Rf_length(initial_points) != 1 || INTEGER(initial_points)[0] < 4 ||
      Rf_length(max_points) != 1 ||
      INTEGER(max_points)[0] < INTEGER(initial_points)[0] ||
      Rf_length(scrambles) != 1 || INTEGER(scrambles)[0] < 2 ||
      Rf_length(relative_tolerance) != 1 ||
      !std::isfinite(REAL(relative_tolerance)[0]) ||
      REAL(relative_tolerance)[0] <= 0.0) {
    Rf_error("Cluster selection controls are invalid.");
  }
  int quadrature_length = 0;
  int previous_order = 0;
  for (int rule = 0; rule < Rf_length(quadrature_orders); ++rule) {
    const double raw_order = REAL(quadrature_orders)[rule];
    const int order = static_cast<int>(raw_order);
    if (!std::isfinite(raw_order) || raw_order != static_cast<double>(order) ||
        order <= previous_order) {
      Rf_error("'quadrature_orders' must contain increasing positive integers.");
    }
    quadrature_length += order;
    previous_order = order;
  }
  if (quadrature_length != Rf_length(quadrature_nodes)) {
    Rf_error("'quadrature_orders' do not match the quadrature rules.");
  }
  const int n_initial_points = INTEGER(initial_points)[0];
  const int n_max_points = INTEGER(max_points)[0];
  const int n_scrambles = INTEGER(scrambles)[0];
  const R_xlen_t expected_qmc = static_cast<R_xlen_t>(2) *
    static_cast<R_xlen_t>(n_max_points) * n_scrambles;
  if (XLENGTH(qmc) != expected_qmc) {
    Rf_error("'qmc' dimensions do not match the cluster integration settings.");
  }
  const int kernel_length = Rf_length(kernel_mode);
  if (kernel_length != 1 && kernel_length != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }
  validate_omega_matrix(REAL(omega), S, B);

  SEXP log_density = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP relative_mcse = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP relative_change = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP log_normalizer = PROTECT(Rf_allocVector(REALSXP, S));
  std::fill_n(REAL(log_normalizer), S, NA_REAL);
  // These rule weights are invariant across posterior draws and observations.
  std::vector<long double> quadrature_weights(quadrature_length);
  for (int i = 0; i < quadrature_length; ++i) {
    quadrature_weights[i] =
      std::exp(static_cast<long double>(REAL(quadrature_log_weights)[i]));
  }
  const int threads = robma_batch_threads(S);
  robma_parallel::for_rows(S, threads, [&](int s) {
    std::vector<double> row_mean(static_cast<std::size_t>(K));
    std::vector<double> row_sd(static_cast<std::size_t>(K));
    std::vector<double> row_coefficient(static_cast<std::size_t>(K));
    std::vector<double> row_omega(static_cast<std::size_t>(B));
    for (int k = 0; k < K; ++k) {
      row_mean[static_cast<std::size_t>(k)] = REAL(means)[s + S * k];
      row_sd[static_cast<std::size_t>(k)] = REAL(residual_sd)[s + S * k];
      row_coefficient[static_cast<std::size_t>(k)] = REAL(loading)[s + S * k];
    }
    for (int b = 0; b < B; ++b) {
      row_omega[static_cast<std::size_t>(b)] = REAL(omega)[s + S * b];
    }
    double mcse = 0.0;
    double error = 0.0;
    REAL(log_density)[s] = cpp_selnorm_cluster_step_lpdf(
      REAL(yi), row_mean.data(), row_sd.data(), row_coefficient.data(), K, REAL(sei),
      row_omega.data(), B, REAL(z_lower), REAL(z_upper), INTEGER(obs_bin),
      INTEGER(sign)[0], LOGICAL(telescope_probabilities)[0],
      scalar_or_row_int(INTEGER(kernel_mode), kernel_length, s, "kernel_mode"),
      REAL(quadrature_nodes), REAL(quadrature_log_weights),
      REAL(quadrature_orders), Rf_length(quadrature_orders),
      REAL(qmc), n_initial_points, n_max_points, n_scrambles,
      REAL(relative_tolerance)[0], &mcse, &error, quadrature_weights.data(),
      LOGICAL(return_normalizer)[0] ? REAL(log_normalizer) + s : nullptr,
      scalar_or_row_int(INTEGER(vector_rule), vector_length, s, "vector_rule")
    );
    REAL(relative_mcse)[s] = mcse;
    REAL(relative_change)[s] = error;
  });

  SEXP output = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_VECTOR_ELT(output, 0, log_density);
  SET_VECTOR_ELT(output, 1, relative_mcse);
  SET_VECTOR_ELT(output, 2, relative_change);
  SET_VECTOR_ELT(output, 3, log_normalizer);
  SET_STRING_ELT(names, 0, Rf_mkChar("log_density"));
  SET_STRING_ELT(names, 1, Rf_mkChar("relative_mcse"));
  SET_STRING_ELT(names, 2, Rf_mkChar("relative_change"));
  SET_STRING_ELT(names, 3, Rf_mkChar("log_normalizer"));
  Rf_setAttrib(output, R_NamesSymbol, names);
  UNPROTECT(6);
  return output;
}

extern "C" SEXP RoBMA_selnorm_factor_step_loglik_batch(
  SEXP yi, SEXP means, SEXP residual_sd, SEXP loading, SEXP sei, SEXP omega,
  SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP telescope_probabilities, SEXP kernel_mode, SEXP quadrature_nodes,
  SEXP quadrature_log_weights, SEXP quadrature_orders,
  SEXP quadrature_rule_counts, SEXP qmc,
  SEXP initial_points, SEXP max_points, SEXP scrambles,
  SEXP relative_tolerance, SEXP return_normalizer, SEXP vector_rule)
{
  check_logical(return_normalizer, "return_normalizer");
  if (LOGICAL(return_normalizer)[0] == NA_LOGICAL) {
    Rf_error("'return_normalizer' must not be missing.");
  }
  check_real(yi, "yi");
  check_real(means, "means");
  check_real(residual_sd, "residual_sd");
  check_real(loading, "loading");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(obs_bin, "obs_bin");
  check_integer(sign, "sign");
  check_logical(telescope_probabilities, "telescope_probabilities");
  check_integer(kernel_mode, "kernel_mode");
  check_real(quadrature_nodes, "quadrature_nodes");
  check_real(quadrature_log_weights, "quadrature_log_weights");
  check_real(quadrature_orders, "quadrature_orders");
  check_real(quadrature_rule_counts, "quadrature_rule_counts");
  check_real(qmc, "qmc");
  check_integer(initial_points, "initial_points");
  check_integer(max_points, "max_points");
  check_integer(scrambles, "scrambles");
  check_real(relative_tolerance, "relative_tolerance");

  const int S = matrix_nrow(means, "means");
  const int vector_length = validate_vector_rule(vector_rule, S);
  const int K = matrix_ncol(means, "means");
  const int loading_columns = matrix_ncol(loading, "loading");
  const int B = matrix_ncol(omega, "omega");
  if (S < 1 || K < 1 || loading_columns % K != 0) {
    Rf_error("Exact factor selection inputs have invalid dimensions.");
  }
  const int rank = loading_columns / K;
  validate_matrix_dims_named(residual_sd, "residual_sd", S, K, "means");
  if (matrix_nrow(loading, "loading") != S || rank < 2 || rank > 4 ||
      B < 1 || matrix_nrow(omega, "omega") != S ||
      Rf_length(yi) != K || Rf_length(sei) != K ||
      Rf_length(obs_bin) != K || Rf_length(z_lower) != B ||
      Rf_length(z_upper) != B) {
    Rf_error("Exact factor selection inputs have invalid dimensions.");
  }
  if (Rf_length(sign) != 1 ||
      (INTEGER(sign)[0] != -1 && INTEGER(sign)[0] != 1) ||
      Rf_length(telescope_probabilities) != 1 ||
      Rf_length(quadrature_nodes) < 1 ||
      Rf_length(quadrature_log_weights) != Rf_length(quadrature_nodes) ||
      Rf_length(quadrature_orders) < 3 ||
      Rf_length(quadrature_rule_counts) != rank - 1 ||
      Rf_length(initial_points) != 1 || INTEGER(initial_points)[0] < 4 ||
      Rf_length(max_points) != 1 ||
      INTEGER(max_points)[0] < INTEGER(initial_points)[0] ||
      Rf_length(scrambles) != 1 || INTEGER(scrambles)[0] < 2 ||
      Rf_length(relative_tolerance) != 1 ||
      !std::isfinite(REAL(relative_tolerance)[0]) ||
      REAL(relative_tolerance)[0] <= 0.0) {
    Rf_error("Exact factor selection controls are invalid.");
  }
  int quadrature_length = 0;
  int previous_order = 0;
  for (int rule = 0; rule < Rf_length(quadrature_orders); ++rule) {
    const double raw_order = REAL(quadrature_orders)[rule];
    const int order = static_cast<int>(raw_order);
    if (!std::isfinite(raw_order) || raw_order != static_cast<double>(order) ||
        order <= previous_order) {
      Rf_error("'quadrature_orders' must contain increasing positive integers.");
    }
    quadrature_length += order;
    previous_order = order;
  }
  if (quadrature_length != Rf_length(quadrature_nodes)) {
    Rf_error("'quadrature_orders' do not match the quadrature rules.");
  }
  for (int index = 0; index < rank - 1; ++index) {
    const double count = REAL(quadrature_rule_counts)[index];
    if (!std::isfinite(count) || count != std::floor(count) || count < 3.0 ||
        count > Rf_length(quadrature_orders)) {
      Rf_error("'quadrature_rule_counts' must contain valid rule-prefix lengths.");
    }
  }
  const int n_initial_points = INTEGER(initial_points)[0];
  const int n_max_points = INTEGER(max_points)[0];
  const int n_scrambles = INTEGER(scrambles)[0];
  const R_xlen_t expected_qmc = static_cast<R_xlen_t>(2 * rank) *
    static_cast<R_xlen_t>(n_max_points) *
    static_cast<R_xlen_t>(n_scrambles);
  if (XLENGTH(qmc) != expected_qmc) {
    Rf_error("'qmc' dimensions do not match the factor integration settings.");
  }
  const int kernel_length = Rf_length(kernel_mode);
  if (kernel_length != 1 && kernel_length != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }
  validate_omega_matrix(REAL(omega), S, B);
  for (int i = 0; i < K; ++i) {
    if (!std::isfinite(REAL(yi)[i]) || !std::isfinite(REAL(sei)[i]) ||
        REAL(sei)[i] <= 0.0 || INTEGER(obs_bin)[i] < 1 ||
        INTEGER(obs_bin)[i] > B) {
      Rf_error("Observed selection inputs are invalid.");
    }
  }
  // The execution plan validates this fixed design at construction. As in
  // JAGS, the shared evaluator checks coordinates only if QMC consumes them;
  // quadrature-only evaluations must not rescan the unused maximum budget.

  SEXP log_density = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP relative_mcse = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP relative_change = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP log_normalizer = PROTECT(Rf_allocVector(REALSXP, S));
  std::fill_n(REAL(log_normalizer), S, NA_REAL);
  // Per-row location and loading inputs are validated before the row loop so
  // parallel workers never raise R errors from a worker thread.
  for (int s = 0; s < S; ++s) {
    for (int i = 0; i < K; ++i) {
      if (!std::isfinite(REAL(means)[s + S * i]) ||
          !(REAL(residual_sd)[s + S * i] > 0.0) ||
          !std::isfinite(REAL(residual_sd)[s + S * i])) {
        Rf_error("Exact factor selection location inputs are invalid.");
      }
    }
    for (int column = 0; column < K * rank; ++column) {
      if (!std::isfinite(REAL(loading)[s + S * column])) {
        Rf_error("Exact factor selection loadings must be finite.");
      }
    }
  }
  const int threads = robma_batch_threads(S);
  robma_parallel::for_rows(S, threads, [&](int s) {
    std::vector<double> row_mean(static_cast<std::size_t>(K));
    std::vector<double> row_sd(static_cast<std::size_t>(K));
    std::vector<double> row_coefficient(
      static_cast<std::size_t>(K * rank)
    );
    std::vector<double> row_omega(static_cast<std::size_t>(B));
    for (int i = 0; i < K; ++i) {
      row_mean[i] = REAL(means)[s + S * i];
      row_sd[i] = REAL(residual_sd)[s + S * i];
    }
    for (int column = 0; column < K * rank; ++column) {
      row_coefficient[column] = REAL(loading)[s + S * column];
    }
    for (int bin = 0; bin < B; ++bin) {
      row_omega[bin] = REAL(omega)[s + S * bin];
    }
    double mcse = 0.0;
    double change = 0.0;
    REAL(log_density)[s] = cpp_selnorm_factor_step_lpdf(
      REAL(yi), row_mean.data(), row_sd.data(), row_coefficient.data(), K, rank,
      REAL(sei), row_omega.data(), B, REAL(z_lower), REAL(z_upper),
      INTEGER(obs_bin), INTEGER(sign)[0],
      LOGICAL(telescope_probabilities)[0],
      scalar_or_row_int(INTEGER(kernel_mode), kernel_length, s, "kernel_mode"),
      REAL(quadrature_nodes), REAL(quadrature_log_weights),
      REAL(quadrature_orders), Rf_length(quadrature_orders),
      REAL(quadrature_rule_counts),
      REAL(qmc), n_initial_points, n_max_points, n_scrambles,
      REAL(relative_tolerance)[0], &mcse, &change,
      LOGICAL(return_normalizer)[0] ? REAL(log_normalizer) + s : nullptr,
      scalar_or_row_int(INTEGER(vector_rule), vector_length, s, "vector_rule")
    );
    REAL(relative_mcse)[s] = mcse;
    REAL(relative_change)[s] = change;
  });

  SEXP output = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_VECTOR_ELT(output, 0, log_density);
  SET_VECTOR_ELT(output, 1, relative_mcse);
  SET_VECTOR_ELT(output, 2, relative_change);
  SET_VECTOR_ELT(output, 3, log_normalizer);
  SET_STRING_ELT(names, 0, Rf_mkChar("log_density"));
  SET_STRING_ELT(names, 1, Rf_mkChar("relative_mcse"));
  SET_STRING_ELT(names, 2, Rf_mkChar("relative_change"));
  SET_STRING_ELT(names, 3, Rf_mkChar("log_normalizer"));
  Rf_setAttrib(output, R_NamesSymbol, names);
  UNPROTECT(6);
  return output;
}

namespace {

bool selnorm_mvn_exact_rank_one_factor(const std::vector<double> &covariance,
                                       int n, double *factor,
                                       double *rank_one)
{
  int pivot = 0;
  double largest = -1.0;
  for (int i = 0; i < n; ++i) {
    const double diagonal = covariance[static_cast<std::size_t>(i + n * i)];
    if (!(diagonal > 0.0)) return false;
    const double magnitude = std::sqrt(diagonal);
    if (magnitude > largest) {
      largest = magnitude;
      pivot = i;
    }
  }

  for (int i = 0; i < n; ++i) {
    const double pivot_value = covariance[static_cast<std::size_t>(pivot + n * i)];
    if (pivot_value == 0.0) return false;
    const double direction = i == pivot ? 1.0 : (pivot_value > 0.0 ? 1.0 : -1.0);
    rank_one[i] =
      direction * std::sqrt(covariance[static_cast<std::size_t>(i + n * i)]);
  }
  for (int column = 0; column < n; ++column) {
    for (int row = 0; row < n; ++row) {
      if (rank_one[row] * rank_one[column] !=
          covariance[static_cast<std::size_t>(row + n * column)]) {
        return false;
      }
    }
  }

  std::fill(factor, factor + static_cast<std::size_t>(n * n), 0.0);
  for (int column = 0; column < n; ++column) {
    factor[n * column] = rank_one[column];
  }
  return true;
}

bool selnorm_mvn_sampling_factor(const std::vector<double> &covariance,
                                 int n, double *factor,
                                 double *eigenvalues, double *work,
                                 int work_length, double *rank_one,
                                 double *cholesky, double *eigenvectors,
                                 bool *positive_definite)
{
  *positive_definite = false;
  // These are exact PSD necessities, independent of eigensolver sign error.
  for (int row = 0; row < n; ++row) {
    const double diagonal = covariance[row + n * row];
    if (diagonal < 0.0) return false;
    if (diagonal == 0.0) {
      for (int column = 0; column < n; ++column) {
        if (covariance[row + n * column] != 0.0 ||
            covariance[column + n * row] != 0.0) return false;
      }
    }
  }
  // Preserve the established one-direction representation before a rounded
  // positive Cholesky pivot can add a spurious direction. Singletons use chol.
  if (n > 1 && selnorm_mvn_exact_rank_one_factor(
        covariance, n, factor, rank_one)) {
    return true;
  }

  int info = 0;
  // Apply the existing numerical PSD roundoff policy in correlation units,
  // so another coordinate's large variance cannot erase a small marginal SD.
  // The exact zero rows were checked above. Reuse existing scratch only.
  int positive = 0;
  for (int row = 0; row < n; ++row) {
    rank_one[row] = std::sqrt(covariance[row + n * row]);
    if (rank_one[row] > 0.0) ++positive;
  }
  std::fill(factor, factor + static_cast<std::size_t>(n * n), 0.0);
  if (positive == 0) return true;
  int scaled_column = 0;
  for (int column = 0; column < n; ++column) {
    if (rank_one[column] == 0.0) continue;
    int scaled_row = 0;
    for (int row = 0; row < n; ++row) {
      if (rank_one[row] == 0.0) continue;
      const double denominator = rank_one[row] * rank_one[column];
      if (!(denominator > 0.0) || !std::isfinite(denominator)) return false;
      const double correlation = row == column ? 1.0 :
        covariance[row + n * column] / denominator;
      if (!std::isfinite(correlation)) return false;
      eigenvectors[scaled_row + positive * scaled_column] = correlation;
      ++scaled_row;
    }
    ++scaled_column;
  }
  F77_CALL(dsyev)("V", "U", &positive, eigenvectors, &positive, eigenvalues,
                  work, &work_length, &info FCONE FCONE);
  if (info != 0) return false;

  double scale = 0.0;
  for (int i = 0; i < positive; ++i) {
    if (!std::isfinite(eigenvalues[i])) return false;
    scale = std::max(scale, std::abs(eigenvalues[i]));
  }
  const double operations = 4.0 * static_cast<double>(positive);
  const double scaled_epsilon = operations * std::numeric_limits<double>::epsilon();
  const double tolerance = scaled_epsilon / (1.0 - scaled_epsilon) * scale;
  if (eigenvalues[0] < -tolerance) return false;

  int numerical_rank = 0;
  for (int i = 0; i < positive; ++i) {
    if (std::abs(eigenvalues[i]) <= tolerance) eigenvalues[i] = 0.0;
    if (eigenvalues[i] < 0.0) return false;
    if (eigenvalues[i] > 0.0) ++numerical_rank;
  }
  std::fill(eigenvalues + positive, eigenvalues + n, 0.0);
  // A rounded positive Cholesky pivot cannot establish a new direction.
  // Reuse its original factor only after the existing SRS policy retains all
  // n directions; otherwise the spectral factor below owns the support.
  if (numerical_rank == n) {
    std::copy(covariance.begin(), covariance.end(), cholesky);
    F77_CALL(dpotrf)("U", &n, cholesky, &n, &info FCONE);
    if (info == 0) {
      std::fill(factor, factor + static_cast<std::size_t>(n * n), 0.0);
      for (int column = 0; column < n; ++column) {
        for (int row = 0; row <= column; ++row) {
          factor[row + n * column] =
            cholesky[row + n * column];
        }
      }
      *positive_definite = true;
      return true;
    }
  }
  scaled_column = 0;
  for (int column = 0; column < n; ++column) {
    if (rank_one[column] == 0.0) continue;
    for (int component = 0; component < positive; ++component) {
      factor[component + n * column] = std::sqrt(eigenvalues[component]) *
        eigenvectors[scaled_column + positive * component] * rank_one[column];
    }
    ++scaled_column;
  }
  return true;
}

int selnorm_mvn_step_bin(double value, double sei, int sign,
                         const double *z_lower, int n_bins)
{
  const double score = static_cast<double>(sign) * value / sei;
  for (int bin = 0; bin < n_bins; ++bin) {
    if (score >= z_lower[bin]) return bin;
  }
  return -1;
}

}

extern "C" SEXP RoBMA_selnorm_mnorm_step_rng_batch(
  SEXP means, SEXP covariance, SEXP sei, SEXP omega,
  SEXP z_lower, SEXP z_upper, SEXP sign, SEXP kernel_mode,
  SEXP dependency_blocks, SEXP max_attempts, SEXP vector_rule)
{
  check_real(means, "means");
  check_real(covariance, "covariance");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(kernel_mode, "kernel_mode");
  check_integer(max_attempts, "max_attempts");

  const int S = matrix_nrow(means, "means");
  const int vector_length = validate_vector_rule(vector_rule, S);
  const int K = matrix_ncol(means, "means");
  const int B = matrix_ncol(omega, "omega");
  if (S < 1 || K < 1 || B < 1 || matrix_nrow(omega, "omega") != S) {
    Rf_error("Selected response inputs have invalid dimensions.");
  }
  SEXP covariance_dim = Rf_getAttrib(covariance, R_DimSymbol);
  if (TYPEOF(covariance_dim) != INTSXP || Rf_length(covariance_dim) != 3 ||
      INTEGER(covariance_dim)[0] != S ||
      INTEGER(covariance_dim)[1] != K ||
      INTEGER(covariance_dim)[2] != K) {
    Rf_error("'covariance' dimensions must match 'means'.");
  }
  if (Rf_length(sei) != K || Rf_length(z_lower) != B ||
      Rf_length(z_upper) != B) {
    Rf_error("Selected response vectors have invalid dimensions.");
  }
  if (Rf_length(sign) != 1 ||
      (INTEGER(sign)[0] != -1 && INTEGER(sign)[0] != 1)) {
    Rf_error("'sign' must be either -1 or 1.");
  }
  const int mode_length = Rf_length(kernel_mode);
  if (mode_length != 1 && mode_length != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(max_attempts) != 1 || INTEGER(max_attempts)[0] < 1) {
    Rf_error("'max_attempts' must be one positive integer.");
  }
  if (TYPEOF(dependency_blocks) != VECSXP || Rf_length(dependency_blocks) < 1) {
    Rf_error("'dependency_blocks' must be a non-empty list.");
  }
  if (!selnorm_is_descending_step_partition(
        REAL(z_lower), REAL(z_upper), B)) {
    Rf_error("Selection bounds do not define a descending step partition.");
  }

  validate_omega_matrix(REAL(omega), S, B);
  std::vector<int> seen(static_cast<std::size_t>(K), 0);
  int maximum_block_size = 0;
  for (int block_index = 0;
       block_index < Rf_length(dependency_blocks); ++block_index) {
    SEXP block = VECTOR_ELT(dependency_blocks, block_index);
    if (TYPEOF(block) != INTSXP || Rf_length(block) < 1) {
      Rf_error("'dependency_blocks' contains an invalid block.");
    }
    maximum_block_size = std::max(maximum_block_size, Rf_length(block));
    for (int j = 0; j < Rf_length(block); ++j) {
      const int index = INTEGER(block)[j];
      if (index < 1 || index > K || seen[static_cast<std::size_t>(index - 1)] != 0) {
        Rf_error("'dependency_blocks' must partition the observations.");
      }
      seen[static_cast<std::size_t>(index - 1)] = block_index + 1;
    }
  }
  for (int k = 0; k < K; ++k) {
    if (seen[static_cast<std::size_t>(k)] == 0) {
      Rf_error("'dependency_blocks' must partition the observations.");
    }
    if (!std::isfinite(REAL(sei)[k]) || REAL(sei)[k] <= 0.0) {
      Rf_error("'sei' must contain finite positive values.");
    }
  }
  for (R_xlen_t i = 0; i < XLENGTH(means); ++i) {
    if (!std::isfinite(REAL(means)[i])) {
      Rf_error("'means' must contain finite values.");
    }
  }
  for (int column = 0; column < K; ++column) {
    for (int row = 0; row < K; ++row) {
      const bool cross_block = seen[row] != seen[column];
      const R_xlen_t offset = static_cast<R_xlen_t>(S) *
        (row + static_cast<R_xlen_t>(K) * column);
      for (int s = 0; s < S; ++s) {
        const double value = REAL(covariance)[offset + s];
        if (!std::isfinite(value)) Rf_error("'covariance' must contain finite values.");
        if (cross_block && value != 0.0) {
          Rf_error("'dependency_blocks' must include all observations with nonzero covariance.");
        }
      }
    }
  }
  for (int s = 0; s < S; ++s) {
    const int mode = scalar_or_row_int(
      INTEGER(kernel_mode), mode_length, s, "kernel_mode"
    );
    if (mode != SELKERNEL_NORMAL && mode != SELKERNEL_STEP) {
      Rf_error("Exact selected response simulation requires a step selection kernel.");
    }
  }

  SEXP draws = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  std::fill(REAL(draws), REAL(draws) + static_cast<R_xlen_t>(S) * K,
            NA_REAL);
  std::vector<double> block_covariance(
    static_cast<std::size_t>(maximum_block_size * maximum_block_size)
  );
  std::vector<double> factor(
    static_cast<std::size_t>(maximum_block_size * maximum_block_size)
  );
  std::vector<double> eigenvalues(static_cast<std::size_t>(maximum_block_size));
  const int work_length = std::max(1, 3 * maximum_block_size - 1);
  std::vector<double> work(static_cast<std::size_t>(work_length));
  std::vector<double> rank_one(static_cast<std::size_t>(maximum_block_size));
  std::vector<double> cholesky(
    static_cast<std::size_t>(maximum_block_size * maximum_block_size)
  );
  std::vector<double> eigenvectors(
    static_cast<std::size_t>(maximum_block_size * maximum_block_size)
  );
  std::vector<double> normal(static_cast<std::size_t>(maximum_block_size));
  std::vector<double> candidate(static_cast<std::size_t>(maximum_block_size));
  std::vector<double> log_weights(static_cast<std::size_t>(B));
  std::vector<double> pivot_log_mass(static_cast<std::size_t>(maximum_block_size) * B);
  std::vector<double> pivot_log_normalizer(static_cast<std::size_t>(maximum_block_size));
  std::vector<double> mixture_log_terms(static_cast<std::size_t>(maximum_block_size));

  int failure_code  = 0;
  int failure_size  = 0;

  std::vector<double> scalar_mass(B), scalar_lower(B), scalar_upper(B);
  SelNormKernelData diagonal_selection = {};
  diagonal_selection.n_bins = B;
  diagonal_selection.effect_sign = INTEGER(sign)[0];
  diagonal_selection.q = 1;
  diagonal_selection.z_lower = REAL(z_lower);
  diagonal_selection.z_upper = REAL(z_upper);
  diagonal_selection.trusted_step_partition = true;
  diagonal_selection.telescope_probabilities = true;

  GetRNGstate();
  for (int s = 0; s < S && failure_code == 0; ++s) {
    const int mode = scalar_or_row_int(
      INTEGER(kernel_mode), mode_length, s, "kernel_mode"
    );
    double maximum_weight = 0.0;
    for (int b = 0; b < B; ++b) {
      maximum_weight = std::max(maximum_weight, REAL(omega)[s + S * b]);
    }
    if (mode == SELKERNEL_STEP && !(maximum_weight > 0.0)) {
      failure_code = 3;
      break;
    }
    const int rule = scalar_or_row_int(INTEGER(vector_rule), vector_length, s, "vector_rule");
    bool constant_weights = true;
    for (int b = 0; b < B; ++b) {
      log_weights[b] = mode == SELKERNEL_NORMAL ? 0.0 :
        std::log(REAL(omega)[s + S * b]) - std::log(maximum_weight);
      constant_weights = constant_weights && REAL(omega)[s + S * b] == maximum_weight;
    }

    for (int block_index = 0;
         block_index < Rf_length(dependency_blocks) && failure_code == 0;
         ++block_index) {
      SEXP block = VECTOR_ELT(dependency_blocks, block_index);
      const int n = Rf_length(block);
      block_covariance.resize(static_cast<std::size_t>(n * n));
      for (int column = 0; column < n; ++column) {
        const int source_column = INTEGER(block)[column] - 1;
        for (int row = 0; row < n; ++row) {
          const int source_row = INTEGER(block)[row] - 1;
          const double a = REAL(covariance)[
            s + S * (source_row + K * source_column)
          ];
          const double b = REAL(covariance)[
            s + S * (source_column + K * source_row)
          ];
          if (a != b) {
            failure_code  = 4;
            failure_size  = n;
            break;
          }
          block_covariance[static_cast<std::size_t>(row + n * column)] =
            a;
        }
        if (failure_code != 0) break;
      }
      if (failure_code != 0) break;
      bool independent_positive = mode == SELKERNEL_STEP && rule == SELVECTOR_PRODUCT;
      for (int column = 0; column < n && independent_positive; ++column) {
        for (int row = 0; row < n; ++row) {
          const double value = block_covariance[row + n * column];
          if (row == column ? !(value > 0.0) : value != 0.0) {
            independent_positive = false;
            break;
          }
        }
      }
      if (independent_positive) {
        // An exactly diagonal Gaussian product law consists of independent
        // scalar selected laws. It needs neither Holder nor factor rejection.
        for (int row = 0; row < n; ++row) {
          const int observation = INTEGER(block)[row] - 1;
          const double bin_uniform = unif_rand();
          const double interval_uniform = unif_rand();
          candidate[row] = cpp_selnorm_kernel_rng_workspace(
            REAL(means)[s + S * observation], std::sqrt(block_covariance[row + n * row]),
            REAL(sei)[observation], REAL(omega) + s, bin_uniform, interval_uniform,
            0.0, 0, SELKERNEL_STEP, diagonal_selection,
            scalar_mass.data(), scalar_lower.data(), scalar_upper.data(), S, false);
          if (!std::isfinite(candidate[row])) {
            failure_code = 6;
            failure_size = n;
            break;
          }
        }
        if (failure_code != 0) break;
        for (int row = 0; row < n; ++row) {
          REAL(draws)[s + S * (INTEGER(block)[row] - 1)] = candidate[row];
        }
        continue;
      }
      const bool deterministic = std::all_of(block_covariance.begin(),
        block_covariance.end(), [](double value) { return value == 0.0; });
      bool positive_definite = false;
      if (!deterministic && !selnorm_mvn_sampling_factor(
            block_covariance, n, factor.data(), eigenvalues.data(),
            work.data(), work_length, rank_one.data(), cholesky.data(),
            eigenvectors.data(), &positive_definite)) {
        failure_code  = 1;
        failure_size  = n;
        break;
      }

      // Holder's inequality bounds the product event by a mixture of n
      // powered marginal weights. Its rejection constant is GM(L_i), where
      // L_i=E[(w_i/max(w))^n], no larger than the old constant one.
      bool holder = positive_definite && n > 1 &&
        mode == SELKERNEL_STEP && rule == SELVECTOR_PRODUCT && !constant_weights;
      for (int row = 0; row < n && holder; ++row) {
        const int observation = INTEGER(block)[row] - 1;
        const double mean = REAL(means)[s + S * observation];
        const double sd = std::sqrt(block_covariance[row + n * row]);
        double total = -std::numeric_limits<double>::infinity();
        for (int bin = 0; bin < B; ++bin) {
          const double lower = REAL(sei)[observation] *
            (INTEGER(sign)[0] > 0 ? REAL(z_lower)[bin] : -REAL(z_upper)[bin]);
          const double upper = REAL(sei)[observation] *
            (INTEGER(sign)[0] > 0 ? REAL(z_upper)[bin] : -REAL(z_lower)[bin]);
          const double mass = n * log_weights[bin] +
            cpp_selnorm_normal_interval_log_prob(lower, upper, mean, sd);
          pivot_log_mass[row * B + bin] = mass;
          total = logspace_add(total, mass);
        }
        pivot_log_normalizer[row] = total;
        if (!std::isfinite(total)) holder = false;
      }

      bool accepted = false;
      for (int attempt = 0;
           attempt < (deterministic ? 1 : INTEGER(max_attempts)[0]) && !accepted; ++attempt) {
        if (holder && attempt == 128) {
          // The preceding Holder attempts have all completed. This optional
          // branch owns the remaining attempts without redrawing the context.
          int status = -2;
          try {
            std::vector<double> box_mean(n), box_se(n), box_omega(B);
            std::vector<double> box_covariance(static_cast<std::size_t>(n) * n);
            for (int row = 0; row < n; ++row) {
              const int observation = INTEGER(block)[row] - 1;
              box_mean[row] = REAL(means)[s + S * observation];
              box_se[row] = REAL(sei)[observation];
              for (int column = 0; column < n; ++column) {
                box_covariance[row + n * column] = block_covariance[
                  std::min(row, column) + n * std::max(row, column)];
              }
            }
            for (int bin = 0; bin < B; ++bin) box_omega[bin] = REAL(omega)[s + S * bin];
            status = cpp_selnorm_factor_box_rng(box_mean.data(), box_covariance.data(), n,
              box_se.data(), box_omega.data(), B, REAL(z_lower), REAL(z_upper), INTEGER(sign)[0],
              INTEGER(max_attempts)[0] - attempt, unif_rand, candidate.data());
          } catch (const std::bad_alloc &) {
            status = -3;
          } catch (const std::exception &) {
            status = -2;
          }
          if (status != -1) {
            accepted = status == 1;
            if (accepted) {
              for (int row = 0; row < n; ++row) {
                REAL(draws)[s + S * (INTEGER(block)[row] - 1)] = candidate[row];
              }
            } else {
              failure_code = status == 0 ? 2 : (status == -3 ? 7 : 6);
              failure_size = n;
            }
            break;
          }
        }
        if (!deterministic) {
          for (int row = 0; row < n; ++row) normal[row] = norm_rand();
        }
        for (int column = 0; column < n; ++column) {
          double value = REAL(means)[
            s + S * (INTEGER(block)[column] - 1)
          ];
          if (!deterministic) {
            const int active_rows = positive_definite ? column + 1 : n;
            for (int row = 0; row < active_rows; ++row) {
              value += normal[row] * factor[row + n * column];
            }
          }
          candidate[column] = value;
        }
        if (holder) {
          const int pivot = static_cast<int>(unif_rand() * n);
          const double log_uniform = std::log(unif_rand()) + pivot_log_normalizer[pivot];
          double cumulative = -std::numeric_limits<double>::infinity();
          int chosen_bin = -1;
          for (int bin = 0; bin < B; ++bin) {
            const double mass = pivot_log_mass[pivot * B + bin];
            if (!std::isfinite(mass)) continue;
            chosen_bin = bin;
            cumulative = logspace_add(cumulative, mass);
            if (log_uniform <= cumulative) break;
          }
          const int observation = INTEGER(block)[pivot] - 1;
          const double lower = REAL(sei)[observation] *
            (INTEGER(sign)[0] > 0 ? REAL(z_lower)[chosen_bin] : -REAL(z_upper)[chosen_bin]);
          const double upper = REAL(sei)[observation] *
            (INTEGER(sign)[0] > 0 ? REAL(z_upper)[chosen_bin] : -REAL(z_lower)[chosen_bin]);
          const double selected_pivot = cpp_selnorm_normal_interval_quantile(
            lower, upper, REAL(means)[s + S * observation],
            std::sqrt(block_covariance[pivot + n * pivot]), unif_rand());
          if (!std::isfinite(selected_pivot)) {
            failure_code = 6;
            failure_size = n;
            break;
          }
          const double displacement = selected_pivot - candidate[pivot];
          for (int row = 0; row < n; ++row) {
            // Match the authoritative triangle used by the existing U-Cholesky
            // of this transposed block (the original input's lower triangle).
            const double covariance_with_pivot = block_covariance[
              std::min(row, pivot) + n * std::max(row, pivot)];
            candidate[row] += covariance_with_pivot /
              block_covariance[pivot + n * pivot] * displacement;
            if (!std::isfinite(candidate[row])) {
              failure_code = 6;
              failure_size = n;
              break;
            }
          }
          if (failure_code != 0) break;
          // Keep the sampled marginal value, including its boundary bin.
          candidate[pivot] = selected_pivot;
        }

        accepted = mode == SELKERNEL_NORMAL;
        if (!accepted) {
          double log_acceptance = 0.0;
          int best_bin = B;
          bool valid_best = true;
          for (int row = 0; row < n; ++row) {
            const int observation = INTEGER(block)[row] - 1;
            if (rule != SELVECTOR_PRODUCT) {
              const double value = rule == SELVECTOR_BEST_TWO_SIDED ?
                std::fabs(candidate[row]) : INTEGER(sign)[0] * candidate[row];
              int bin = 0;
              while (bin < B && !(value >= REAL(sei)[observation] * REAL(z_lower)[bin])) ++bin;
              if (bin == B) { valid_best = false; break; }
              best_bin = std::min(best_bin, bin);
              continue;
            }
            const int bin = selnorm_mvn_step_bin(
              candidate[row], REAL(sei)[observation], INTEGER(sign)[0],
              REAL(z_lower), B
            );
            if (bin < 0 || (holder && !std::isfinite(log_weights[bin]))) {
              log_acceptance = -std::numeric_limits<double>::infinity();
              break;
            }
            log_acceptance += log_weights[bin];
            if (holder) mixture_log_terms[row] = n * log_weights[bin] - pivot_log_normalizer[row];
          }
          if (holder && std::isfinite(log_acceptance)) {
            const double maximum = *std::max_element(mixture_log_terms.begin(),
              mixture_log_terms.begin() + n);
            long double sum = 0.0L;
            long double sum_expm1 = 0.0L;
            for (int row = 0; row < n; ++row) {
              const long double centered = static_cast<long double>(mixture_log_terms[row]) - maximum;
              sum += centered;
              sum_expm1 += std::expm1(centered);
            }
            // mean(t)-log(mean(exp(t))) is stable after centering; Jensen's
            // inequality supplies the rejection bound without clamping.
            log_acceptance = static_cast<double>(sum / n - std::log1p(sum_expm1 / n));
          }
          if (rule != SELVECTOR_PRODUCT) {
            log_acceptance = valid_best ?
              std::log(REAL(omega)[s + S * best_bin]) - std::log(maximum_weight) :
              -std::numeric_limits<double>::infinity();
          }
          accepted = deterministic ? std::isfinite(log_acceptance) :
            std::log(unif_rand()) <= log_acceptance;
        }
        if (accepted) {
          for (int row = 0; row < n; ++row) {
            REAL(draws)[s + S * (INTEGER(block)[row] - 1)] = candidate[row];
          }
        }
      }
      if (!accepted && failure_code == 0) {
        failure_code  = deterministic ? 5 : 2;
        failure_size  = n;
      }
    }
  }
  PutRNGstate();
  if (failure_code == 7) {
    Rf_error("Selected response simulation could not allocate its proposal workspace.");
  }

  if (failure_code == 5) {
    Rf_error("Selected response simulation is unavailable because a fully retained outcome has zero acceptance probability. Use strictly positive selection weights or integrate an outcome-generating source.");
  }

  SEXP output = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names  = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_VECTOR_ELT(output, 0, draws);
  SET_VECTOR_ELT(output, 1, Rf_ScalarInteger(failure_code));
  SET_VECTOR_ELT(output, 2, Rf_ScalarInteger(failure_size));
  SET_STRING_ELT(names, 0, Rf_mkChar("draws"));
  SET_STRING_ELT(names, 1, Rf_mkChar("failure_code"));
  SET_STRING_ELT(names, 2, Rf_mkChar("failure_size"));
  Rf_setAttrib(output, R_NamesSymbol, names);

  UNPROTECT(3);
  return output;
}

extern "C" SEXP RoBMA_selnorm_covariance_envelope_components(
    SEXP covariance_lower, SEXP selection_se, SEXP upper)
{
  check_real(covariance_lower, "covariance_lower");
  check_real(selection_se, "selection_se");
  check_logical(upper, "upper");
  const int k = Rf_length(selection_se);
  if (k < 1 || k > 46340 ||
      XLENGTH(covariance_lower) != static_cast<R_xlen_t>(k) * (k + 1) / 2 ||
      XLENGTH(upper) != 1 || LOGICAL(upper)[0] == NA_LOGICAL) {
    Rf_error("Compact envelope dimensions or side are invalid.");
  }
  SEXP result = PROTECT(Rf_allocVector(VECSXP, 10));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 10));
  const char *labels[] = {"available", "reason", "residual_variance", "residual_sd",
    "loading", "covariance_residual", "groups", "type_index", "upper",
    "covariance_residual_radius"};
  for (int i = 0; i < 10; ++i) SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  Rf_setAttrib(result, R_NamesSymbol, names);
  SEXP sd = PROTECT(Rf_allocVector(REALSXP, k));
  SEXP work = PROTECT(Rf_allocMatrix(REALSXP, k, 4));
  SEXP types = PROTECT(Rf_allocVector(INTSXP, k));
  int groups = 0;
  bool available = false;
  bool allocation_failed = false;
  try {
    available = cpp_selnorm_covariance_envelope_components(
      REAL(covariance_lower), k, REAL(selection_se), LOGICAL(upper)[0],
      REAL(sd), REAL(work), &groups, INTEGER(types));
  } catch (const std::bad_alloc &) {
    allocation_failed = true;
  }
  if (allocation_failed) Rf_error("Compact envelope allocation failed.");
  SET_VECTOR_ELT(result, 0, Rf_ScalarLogical(available));
  SET_VECTOR_ELT(result, 8, upper);
  if (!available) {
    SET_VECTOR_ELT(result, 1, Rf_mkString(
      "No admissible positive-residual covariance envelope is available."));
    UNPROTECT(5);
    return result;
  }
  const int rank = groups == 1 ? 1 : groups + 1;
  SEXP variance = PROTECT(Rf_allocVector(REALSXP, k));
  SEXP loading = PROTECT(Rf_allocMatrix(REALSXP, k, rank));
  SEXP residual = PROTECT(Rf_allocMatrix(REALSXP, k, k));
  SEXP radius = PROTECT(Rf_allocMatrix(REALSXP, k, k));
  std::copy_n(REAL(work), static_cast<std::size_t>(k) * rank, REAL(loading));
  for (int row = 0; row < k; ++row) {
    REAL(variance)[row] = static_cast<double>(
      static_cast<long double>(REAL(sd)[row]) * REAL(sd)[row]);
    INTEGER(types)[row] += 1;
  }
  // Ctilde is defined by the actual emitted double SDs and loadings, not by
  // nominal envelope coefficients or a rounded reconstruction of Ctilde.
  const long double infinity = std::numeric_limits<long double>::infinity();
  const auto operation_radius = [infinity](long double value) {
    if (!std::isfinite(value)) return infinity;
    return std::max(std::nextafter(value, infinity) - value,
                    value - std::nextafter(value, -infinity));
  };
  const auto add_radius = [infinity](long double left, long double right) {
    return std::nextafter(left + right, infinity);
  };
  std::size_t position = 0;
  for (int column = 0; column < k; ++column) {
    for (int row = column; row < k; ++row) {
      long double represented = row == column ?
        static_cast<long double>(REAL(sd)[row]) * REAL(sd)[row] : 0.0L;
      long double error = row == column ? operation_radius(represented) : 0.0L;
      for (int factor = 0; factor < rank; ++factor) {
        const long double product = static_cast<long double>(REAL(loading)[row + k * factor]) *
          REAL(loading)[column + k * factor];
        error = add_radius(error, operation_radius(product));
        represented += product;
        error = add_radius(error, operation_radius(represented));
      }
      const long double exact_input_difference =
        static_cast<long double>(REAL(covariance_lower)[position++]) - represented;
      error = add_radius(error, operation_radius(exact_input_difference));
      const double difference = static_cast<double>(exact_input_difference);
      const double double_infinity = std::numeric_limits<double>::infinity();
      const double cast_error = std::max(
        std::nextafter(difference, double_infinity) - difference,
        difference - std::nextafter(difference, -double_infinity));
      error = add_radius(error, static_cast<long double>(cast_error));
      const double returned_radius = std::nextafter(static_cast<double>(error), double_infinity);
      REAL(residual)[row + k * column] = difference;
      REAL(residual)[column + k * row] = difference;
      REAL(radius)[row + k * column] = returned_radius;
      REAL(radius)[column + k * row] = returned_radius;
    }
  }
  SET_VECTOR_ELT(result, 2, variance);
  SET_VECTOR_ELT(result, 3, sd);
  SET_VECTOR_ELT(result, 4, loading);
  SET_VECTOR_ELT(result, 5, residual);
  SET_VECTOR_ELT(result, 6, Rf_ScalarInteger(groups));
  SET_VECTOR_ELT(result, 7, types);
  SET_VECTOR_ELT(result, 9, radius);
  UNPROTECT(9);
  return result;
}

extern "C" SEXP RoBMA_selnorm_context_star_zplot(
  SEXP means, SEXP covariance_lower, SEXP selection_se, SEXP omega,
  SEXP lower, SEXP upper, SEXP sign, SEXP telescope, SEXP log_normalizers,
  SEXP context_log_weights, SEXP nodes, SEXP log_weights, SEXP z, SEXP allowance)
{
  check_real(means, "means"); check_real(covariance_lower, "covariance_lower");
  check_real(selection_se, "selection_se"); check_real(omega, "omega");
  check_real(lower, "lower"); check_real(upper, "upper"); check_integer(sign, "sign");
  check_logical(telescope, "telescope");
  if (log_normalizers != R_NilValue) check_real(log_normalizers, "log_normalizers");
  check_real(context_log_weights, "context_log_weights"); check_real(nodes, "nodes");
  check_real(log_weights, "log_weights"); check_real(z, "z"); check_real(allowance, "allowance");
  const int contexts = matrix_nrow(means, "means");
  const int k = matrix_ncol(means, "means");
  const int bins = Rf_length(omega), order = Rf_length(nodes), size = Rf_length(z);
  if (contexts < 1 || k < 1 || bins < 1 || order < 1 || size < 1 ||
      XLENGTH(covariance_lower) != static_cast<R_xlen_t>(k) * (k + 1) / 2 ||
      Rf_length(selection_se) != k || Rf_length(lower) != bins || Rf_length(upper) != bins ||
      Rf_length(sign) != 1 || (INTEGER(sign)[0] != -1 && INTEGER(sign)[0] != 1) ||
      Rf_length(telescope) != 1 || LOGICAL(telescope)[0] == NA_LOGICAL ||
      (log_normalizers != R_NilValue && Rf_length(log_normalizers) != contexts) ||
      Rf_length(context_log_weights) != contexts ||
      Rf_length(log_weights) != order || Rf_length(allowance) != 1 ||
      !std::isfinite(REAL(allowance)[0]) || REAL(allowance)[0] < 0.0) {
    Rf_error("Deterministic star projection dimensions or controls are invalid.");
  }
  const auto finite = [](SEXP values) {
    for (R_xlen_t i = 0; i < XLENGTH(values); ++i) if (!std::isfinite(REAL(values)[i])) return false;
    return true;
  };
  if (!finite(means) || (log_normalizers != R_NilValue && !finite(log_normalizers)) ||
      !finite(nodes) || !finite(log_weights) || !finite(z)) {
    Rf_error("Deterministic star projection locations and quadrature values must be finite.");
  }
  for (int context = 0; context < contexts; ++context) {
    if (std::isnan(REAL(context_log_weights)[context]) || REAL(context_log_weights)[context] == R_PosInf) {
      Rf_error("Context log weights must be finite or negative infinity.");
    }
  }
  for (int row = 0; row < k; ++row) {
    if (!std::isfinite(REAL(selection_se)[row]) || !(REAL(selection_se)[row] > 0)) {
      Rf_error("'selection_se' must be finite and positive.");
    }
  }
  validate_omega_matrix(REAL(omega), 1, bins);
  if (!selnorm_is_descending_step_partition(REAL(lower), REAL(upper), bins)) {
    Rf_error("Selection bounds do not define a descending step partition.");
  }
  SEXP cdf_policy = Rf_getAttrib(nodes, Rf_install("bounded_cdf"));
  bool bounded_cdf = false;
  if (cdf_policy != R_NilValue) {
    if (TYPEOF(cdf_policy) != LGLSXP || Rf_length(cdf_policy) != 1 ||
        LOGICAL(cdf_policy)[0] == NA_LOGICAL) {
      Rf_error("'bounded_cdf' must be one nonmissing logical value.");
    }
    bounded_cdf = LOGICAL(cdf_policy)[0] == TRUE;
  }
  SEXP result = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 6));
  SEXP density = PROTECT(Rf_allocVector(REALSXP, size));
  SEXP compact = PROTECT(Rf_allocVector(REALSXP, contexts));
  std::fill_n(REAL(compact), contexts, NA_REAL);
  SEXP cdf_factor = PROTECT(Rf_allocVector(REALSXP, contexts));
  SEXP cdf_mass = PROTECT(Rf_allocVector(REALSXP, contexts));
  SEXP cdf_numerator = PROTECT(Rf_allocVector(REALSXP, contexts));
  std::fill_n(REAL(cdf_factor), contexts, 0.0);
  const char *labels[] = {"available", "density", "mass", "compression_error", "mixture_components", "compact_log_normalizers"};
  for (int i = 0; i < 6; ++i) SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  Rf_setAttrib(result, R_NamesSymbol, names);
  double mass = NA_REAL, error = NA_REAL, omitted_mass_error = 0.0;
  std::size_t components = 0;
  bool available = false;
  try {
    CompactStepPartition partition;
    const bool compact_partition = bounded_cdf && LOGICAL(telescope)[0] &&
      partition.prepare(REAL(omega), REAL(lower), REAL(upper), bins);
    available = cpp_selnorm_context_star_projection(REAL(means), contexts, REAL(covariance_lower), k,
      REAL(selection_se), compact_partition ? partition.omega.data() : REAL(omega),
      compact_partition ? static_cast<int>(partition.omega.size()) : bins,
      compact_partition ? partition.lower.data() : REAL(lower),
      compact_partition ? partition.upper.data() : REAL(upper), INTEGER(sign)[0], LOGICAL(telescope)[0],
      log_normalizers == R_NilValue ? nullptr : REAL(log_normalizers),
      REAL(context_log_weights), REAL(nodes), REAL(log_weights), order,
      REAL(z), size, REAL(allowance)[0], REAL(density), &mass, &error, &components, REAL(compact),
      bounded_cdf, REAL(cdf_factor), &omitted_mass_error);
  } catch (const std::exception &exception) {
    Rf_error("Deterministic star projection failed: %s.", exception.what());
  }
  SET_VECTOR_ELT(result, 0, Rf_ScalarLogical(available));
  if (available) SET_VECTOR_ELT(result, 1, density);
  SET_VECTOR_ELT(result, 2, Rf_ScalarReal(mass));
  SEXP compression = PROTECT(Rf_ScalarReal(error));
  SEXP omission = PROTECT(Rf_ScalarReal(omitted_mass_error));
  Rf_setAttrib(compression, Rf_install("omitted_mass_error"), omission);
  SET_VECTOR_ELT(result, 3, compression);
  SET_VECTOR_ELT(result, 4, Rf_ScalarReal(static_cast<double>(components)));
  SET_VECTOR_ELT(result, 5, compact);
  const auto multiple_bound = [](double bound, int count) {
    return bound > 0 && count > 0 ?
      std::nextafter(bound * count, std::numeric_limits<double>::infinity()) : 0.0;
  };
  for (int context = 0; context < contexts; ++context) {
    REAL(cdf_mass)[context] = multiple_bound(REAL(cdf_factor)[context], k);
    REAL(cdf_numerator)[context] = multiple_bound(REAL(cdf_factor)[context], k - 1);
  }
  // These describe this rule's own raw-A and numerator primitives only.
  // Supplied-denominator uncertainty belongs to the R controller. Attributes
  // remain available for complete mass-only tables with available == FALSE.
  Rf_setAttrib(compact, Rf_install("cdf_factor_log_error"), cdf_factor);
  Rf_setAttrib(compact, Rf_install("cdf_log_error"), cdf_mass);
  Rf_setAttrib(compact, Rf_install("cdf_numerator_log_error"), cdf_numerator);
  UNPROTECT(9);
  return result;
}
