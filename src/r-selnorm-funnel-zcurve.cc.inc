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
        dens_sum += cpp_selnorm_affine_normal_pdf_log_scale(
          z_p[ell], sei_p[k], mean_p[cell], sd_p[cell],
          log_scale[static_cast<size_t>(cell)]
        );
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

namespace {

// Inverse-normalizer interpolation on first-kind Chebyshev roots.
// Only H is interpolated. Gaussian densities and original step weights remain
// in the caller. The analytic remainder inherits ordinary floating arithmetic.
struct ZcurveChebyshevGrid {
  int count;
  long double half_log_factorial;
  bool valid = true;
  std::vector<double> nodes;
  std::vector<long double> cardinal;
  bool double_available = true;
  std::vector<double> cardinal_double;
  std::vector<long double> cardinal_cast_sum, cardinal_absolute_sum;
  long double minimum_cardinal = std::numeric_limits<long double>::infinity();
  long double maximum_cardinal = 0;

  ZcurveChebyshevGrid(const double *grid, int size, int n, double lower, double upper)
    : count(n), half_log_factorial(0.5L * std::lgamma(static_cast<long double>(n + 1))),
      nodes(n), cardinal(static_cast<std::size_t>(size) * n)
  {
    const long double pi = std::acos(-1.0L);
    const long double midpoint = (static_cast<long double>(lower) + upper) / 2;
    const long double radius = (static_cast<long double>(upper) - lower) / 2;
    std::vector<long double> weights(n), terms(n);
    for (int j = 0; j < n; ++j) {
      const long double angle = pi * (2 * j + 1) / (2 * n);
      nodes[j] = static_cast<double>(midpoint + radius * std::cos(angle));
      weights[j] = (j % 2 ? -1 : 1) * std::sin(angle);
      if (!std::isfinite(nodes[j]) || (j > 0 && !(nodes[j] < nodes[j - 1]))) {
        valid = false;
        return;
      }
    }
    for (int point = 0; point < size; ++point) {
      int match = -1;
      for (int j = 0; j < n; ++j) if (grid[point] == nodes[j]) { match = j; break; }
      if (match >= 0) {
        cardinal[static_cast<std::size_t>(point) * n + match] = 1.0L;
        continue;
      }
      long double denominator = 0.0L;
      for (int j = 0; j < n; ++j) {
        terms[j] = weights[j] / (static_cast<long double>(grid[point]) - nodes[j]);
        denominator += terms[j];
      }
      if (!std::isfinite(denominator) || denominator == 0.0L) {
        valid = false;
        return;
      }
      for (int j = 0; j < n; ++j) {
        cardinal[static_cast<std::size_t>(point) * n + j] = terms[j] / denominator;
      }
    }
  
    cardinal_double.resize(cardinal.size());
    cardinal_cast_sum.assign(size, 0.0L);
    cardinal_absolute_sum.assign(size, 0.0L);
    const long double operations = (2 * count + 2) * std::numeric_limits<long double>::epsilon();
    const long double table_gamma = operations / (1 - operations);
    for (int point = 0; point < size; ++point) {
      for (int j = 0; j < count; ++j) {
        const std::size_t index = static_cast<std::size_t>(point) * count + j;
        const double coefficient = static_cast<double>(cardinal[index]);
        if (!std::isfinite(coefficient) || (coefficient == 0 && cardinal[index] != 0)) {
          double_available = false;
          return;
        }
        cardinal_double[index] = coefficient;
        const long double magnitude = std::fabs(static_cast<long double>(coefficient));
        cardinal_cast_sum[point] += std::fabs(cardinal[index] - static_cast<long double>(coefficient));
        cardinal_absolute_sum[point] += magnitude;
        if (magnitude > 0) minimum_cardinal = std::min(minimum_cardinal, magnitude);
        maximum_cardinal = std::max(maximum_cardinal, magnitude);
      }
      cardinal_cast_sum[point] /= 1 - table_gamma;
      cardinal_absolute_sum[point] /= 1 - table_gamma;
    }
  }
};

// The working cardinal table is rounded once. Charges below refer to the
// stored long-double polynomial, independently of its analytic H remainder.
struct ZcurveDoubleInterpolation {
  const ZcurveChebyshevGrid &grid;
  std::vector<double> H, error;
  long double maximum_H = 0, H_cast_error = 0, maximum_node_error = 0;
  long double gamma = 0;
  bool valid = true;

  static double upward(long double value)
  {
    double out = static_cast<double>(value);
    if (static_cast<long double>(out) < value) {
      out = std::nextafter(out, std::numeric_limits<double>::infinity());
    }
    return out;
  }

  ZcurveDoubleInterpolation(const ZcurveChebyshevGrid &geometry,
      const std::vector<long double> &values, const std::vector<long double> &errors)
    : grid(geometry), H(values.size()), error(errors.size())
  {
    valid = grid.double_available;
    long double minimum_H = std::numeric_limits<long double>::infinity();
    long double minimum_error = std::numeric_limits<long double>::infinity();
    long double maximum_error = 0;
    for (std::size_t j = 0; j < H.size() && valid; ++j) {
      H[j] = static_cast<double>(values[j]);
      error[j] = upward(errors[j]);
      valid = H[j] > 0 && std::isfinite(H[j]) && error[j] >= 0 && std::isfinite(error[j]);
      maximum_H = std::max(maximum_H, values[j]);
      minimum_H = std::min(minimum_H, static_cast<long double>(H[j]));
      H_cast_error = std::max(H_cast_error, std::fabs(values[j] - static_cast<long double>(H[j])));
      maximum_node_error = std::max(maximum_node_error, errors[j]);
      if (error[j] > 0) minimum_error = std::min(minimum_error, static_cast<long double>(error[j]));
      maximum_error = std::max(maximum_error, static_cast<long double>(error[j]));
    }
    const long double normal_min = std::numeric_limits<double>::min();
    const long double safe_max = std::numeric_limits<double>::max() / (2.0L * grid.count);
    // Exclude product underflow/overflow before entering the vectorizable loop.
    valid = valid && grid.minimum_cardinal * minimum_H >= normal_min &&
      grid.maximum_cardinal * maximum_H <= safe_max;
    if (maximum_error > 0) {
      valid = valid && grid.minimum_cardinal * minimum_error >= normal_min &&
        grid.maximum_cardinal * maximum_error <= safe_max;
    }
    const long double operations = 2 * grid.count + 8;
    const long double bound = operations * std::numeric_limits<double>::epsilon();
    gamma = bound / (1 - bound);
  }

  bool evaluate(int point, long double *value, long double *estimated_error,
                long double *rounding_error) const
  {
    if (!valid) return false;
    const double *coefficient = grid.cardinal_double.data() + static_cast<std::size_t>(point) * grid.count;
    double values[4] = {}, absolutes[4] = {}, estimates[4] = {};
    int j = 0;
    for (; j + 3 < grid.count; j += 4) {
      for (int lane = 0; lane < 4; ++lane) {
        const double product = coefficient[j + lane] * H[j + lane];
        values[lane] += product;
        absolutes[lane] += std::fabs(product);
        estimates[lane] += std::fabs(coefficient[j + lane]) * error[j + lane];
      }
    }
    for (; j < grid.count; ++j) {
      const double product = coefficient[j] * H[j];
      values[0] += product;
      absolutes[0] += std::fabs(product);
      estimates[0] += std::fabs(coefficient[j]) * error[j];
    }
    const double sum = (values[0] + values[1]) + (values[2] + values[3]);
    const double absolute = (absolutes[0] + absolutes[1]) + (absolutes[2] + absolutes[3]);
    const double estimate = (estimates[0] + estimates[1]) + (estimates[2] + estimates[3]);
    *value = sum;
    *rounding_error = grid.cardinal_cast_sum[point] * maximum_H +
      grid.cardinal_absolute_sum[point] * H_cast_error + gamma / (1 - gamma) * absolute;
    *estimated_error = estimate / (1 - gamma) +
      grid.cardinal_cast_sum[point] * maximum_node_error;
    return std::isfinite(sum) && std::isfinite(*rounding_error) && std::isfinite(*estimated_error);
  }
};

}

extern "C" SEXP RoBMA_selnorm_zcurve_density_matrix(
  SEXP z_sequence, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP extrapolate, SEXP telescope_probabilities,
  SEXP latent_sd, SEXP quadrature)
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

  const bool latent_mixture = latent_sd != R_NilValue;
  const double *latent_p = NULL;
  const double *nodes_p = NULL;
  const double *log_weights_p = NULL;
  int n_nodes = 0;
  if (latent_mixture) {
    check_real(latent_sd, "latent_sd");
    validate_matrix_dims_named(latent_sd, "latent_sd", S, K, "mean");
    if (TYPEOF(quadrature) != VECSXP || Rf_length(quadrature) != 3) {
      Rf_error("'quadrature' must contain nodes, weights, and log weights.");
    }
    SEXP nodes = VECTOR_ELT(quadrature, 0);
    SEXP log_weights = VECTOR_ELT(quadrature, 2);
    check_real(nodes, "quadrature nodes");
    check_real(log_weights, "quadrature log weights");
    n_nodes = Rf_length(nodes);
    if (n_nodes < 1 || Rf_length(log_weights) != n_nodes) {
      Rf_error("Quadrature nodes and log weights must have equal positive length.");
    }
    latent_p = REAL(latent_sd);
    nodes_p = REAL(nodes);
    log_weights_p = REAL(log_weights);
    for (int j = 0; j < n_nodes; ++j) {
      if (!std::isfinite(nodes_p[j]) || !std::isfinite(log_weights_p[j])) {
        Rf_error("Quadrature nodes and log weights must be finite.");
      }
    }
  } else if (quadrature != R_NilValue) {
    Rf_error("'quadrature' requires 'latent_sd'.");
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
  std::vector<double> total_sd(latent_mixture ? S * K : 0);
  std::vector<double> conditional_sd(latent_mixture ? S * K : 0);
  std::vector<double> context_fraction(latent_mixture ? S * K : 0);
  std::vector<double> candidate_fraction(latent_mixture ? S * K : 0);
  for (int k = 0; k < K; ++k) {
    const double log_sei = std::log(sei_p[k]);
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      if (latent_mixture) {
        if (!(sd_p[cell] > 0) || !std::isfinite(sd_p[cell]) ||
            !(latent_p[cell] >= 0) || !std::isfinite(latent_p[cell])) {
          Rf_error("Latent z densities require positive integrated SDs and nonnegative retained SDs.");
        }
        const double total = std::hypot(sd_p[cell], latent_p[cell]);
        const double retained_ratio = latent_p[cell] / total;
        const double integrated_ratio = sd_p[cell] / total;
        total_sd[cell] = total;
        conditional_sd[cell] = sd_p[cell] * retained_ratio;
        context_fraction[cell] = retained_ratio * retained_ratio;
        candidate_fraction[cell] = integrated_ratio * integrated_ratio;
      }
      log_scale[static_cast<size_t>(cell)] =
        log_sei - std::log(latent_mixture ? total_sd[cell] : sd_p[cell]);

      if (row_active_phack[static_cast<size_t>(s)]) {
        log_norm[cell] = std::numeric_limits<double>::quiet_NaN();
      } else if (row_use_step[static_cast<size_t>(s)] &&
                 (!latent_mixture || latent_p[cell] == 0)) {
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

  if (latent_mixture) {
    SelNormNormalProjectionWorkspace projection(z_p, L);
    std::vector<double> row_weights(B), node_means(n_nodes), log_coefficients(n_nodes);
    std::vector<long double> normalizers(n_nodes), quadrature_weights(n_nodes);
    std::vector<long double> density(L), cell_density(L);
    SEXP interpolation = Rf_getAttrib(quadrature, Rf_install("density_control"));
    const bool requested = interpolation != R_NilValue;
    bool interpolate = requested && L >= 16;
    double interpolation_tolerance = 0.0;
    const double *coarse_nodes = nullptr, *coarse_logs = nullptr;
    int coarse_count = 0;
    std::vector<long double> coarse_weights;
    SEXP fitted_error = R_NilValue, extrapolated_error = R_NilValue;
    if (requested) {
      if (TYPEOF(interpolation) != VECSXP || Rf_length(interpolation) != 2) {
        Rf_error("Density integration controls must supply a relative tolerance and a coarse quadrature rule.");
      }
      SEXP tolerance = VECTOR_ELT(interpolation, 0), coarse = VECTOR_ELT(interpolation, 1);
      check_real(tolerance, "interpolation tolerance");
      if (Rf_length(tolerance) != 1 || !std::isfinite(REAL(tolerance)[0]) || !(REAL(tolerance)[0] > 0)) {
        Rf_error("Density integration tolerance must be finite and positive.");
      }
      interpolation_tolerance = REAL(tolerance)[0];
      if (coarse == R_NilValue) Rf_error("Density integration requires a coarse quadrature rule.");
      if (coarse != R_NilValue) {
        if (TYPEOF(coarse) != VECSXP || Rf_length(coarse) != 3) Rf_error("Invalid coarse quadrature rule.");
        SEXP cn = VECTOR_ELT(coarse, 0), cw = VECTOR_ELT(coarse, 2);
        check_real(cn, "coarse nodes"); check_real(cw, "coarse log weights");
        coarse_count = Rf_length(cn);
        if (coarse_count < 1 || coarse_count >= n_nodes || Rf_length(cw) != coarse_count) {
          Rf_error("Coarse quadrature dimensions are invalid.");
        }
        coarse_nodes = REAL(cn); coarse_logs = REAL(cw);
        coarse_weights.resize(coarse_count);
        for (int j = 0; j < coarse_count; ++j) {
          if (!std::isfinite(coarse_nodes[j]) || !std::isfinite(coarse_logs[j])) {
            Rf_error("Coarse quadrature values must be finite.");
          }
          coarse_weights[j] = std::exp(static_cast<long double>(coarse_logs[j]));
        }
      }
      fitted_error = PROTECT(Rf_allocVector(REALSXP, S)); ++n_protect;
      extrapolated_error = PROTECT(Rf_allocVector(REALSXP, S)); ++n_protect;
    }
    std::vector<ZcurveChebyshevGrid> chebyshev;
    long double log_interval_quarter = 0.0L;
    if (interpolate) {
      double lower = z_p[0], upper = z_p[0];
      for (int ell = 0; ell < L; ++ell) {
        if (!std::isfinite(z_p[ell])) { interpolate = false; break; }
        lower = std::min(lower, z_p[ell]); upper = std::max(upper, z_p[ell]);
      }
      if (!(upper > lower)) interpolate = false;
      if (interpolate) {
        log_interval_quarter = std::log((static_cast<long double>(upper) - lower) / 4);
        for (int n : {8, 12, 16, 24, 32, 48, 64}) {
          if (2 * n < L) chebyshev.emplace_back(z_p, L, n, lower, upper);
        }
      }
    }
    std::vector<long double> interpolation_error(L), noninterpolated(L), cell_error(L);
    std::vector<long double> coarse_noninterpolated(L), cell_coarse(L), cell_gaussian(L);
    const long double half_log_two_pi = 0.5L * std::log(2 * std::acos(-1.0L));

    for (int j = 0; j < n_nodes; ++j) {
      quadrature_weights[j] = std::exp(static_cast<long double>(log_weights_p[j]));
    }
    for (int s = 0; s < S; ++s) {
      std::fill(density.begin(), density.end(), 0.0L);
      std::fill(interpolation_error.begin(), interpolation_error.end(), 0.0L);
      std::fill(noninterpolated.begin(), noninterpolated.end(), 0.0L);
      std::fill(coarse_noninterpolated.begin(), coarse_noninterpolated.end(), 0.0L);
      long double reference_peak = 0.0L;
      long double fitted_absolute_error = 0.0L, extrapolated_absolute_error = 0.0L;
      long double fitted_peak = 0.0L, extrapolated_peak = 0.0L;
      bool zero_weight_display = true;
      if (interpolate) {
        // Actual requested grid values give a lower bound on the Gaussian
        // grid peak and therefore only tighten the interpolation allocation.
        long double highest = -1.0L, anchor = 0.0L, total_weight = 0.0L, weighted_mean = 0.0L;
        for (int k = 0; k < K; ++k) {
          const int cell = s + S * k;
          const long double weight = static_cast<long double>(sei_p[k]) / total_sd[cell];
          if (weight > highest) {
            highest = weight;
            anchor = static_cast<long double>(mean_p[cell]) / sei_p[k];
          }
          total_weight += weight;
          weighted_mean += static_cast<long double>(mean_p[cell]) / total_sd[cell];
        }
        const long double central = total_weight > 0 ? weighted_mean / total_weight : anchor;
        const auto nearest = [&](long double location) {
          int selected = 0;
          long double distance = std::fabs(static_cast<long double>(z_p[0]) - location);
          for (int ell = 1; ell < L; ++ell) {
            const long double candidate = std::fabs(static_cast<long double>(z_p[ell]) - location);
            if (candidate < distance) { selected = ell; distance = candidate; }
          }
          return selected;
        };
        const auto evaluate_reference = [&](int ell) {
          long double value = 0.0L;
          for (int k = 0; k < K; ++k) {
            const int cell = s + S * k;
            value += std::exp(static_cast<long double>(cpp_selnorm_affine_normal_lpdf_log_scale(
              z_p[ell], sei_p[k], mean_p[cell], total_sd[cell], log_scale[cell])));
          }
          return value / K;
        };
        const int first = nearest(anchor), second = nearest(central);
        reference_peak = evaluate_reference(first);
        if (second != first) reference_peak = std::max(reference_peak, evaluate_reference(second));
        // Extreme separated locations may make both chosen scores vanish.
        // Retain the original grid census in that exceptional case.
        if (!(reference_peak > 0) || !std::isfinite(reference_peak)) {
          reference_peak = 0.0L;
          for (int ell = 0; ell < L; ++ell) reference_peak = std::max(reference_peak, evaluate_reference(ell));
        }
      }

      for (int b = 0; b < B; ++b) row_weights[b] = omega_p[s + S * b];
      const bool selected = row_use_step[static_cast<size_t>(s)];
      const double minimum_weight = *std::min_element(row_weights.begin(), row_weights.end());
      const double maximum_weight = *std::max_element(row_weights.begin(), row_weights.end());
      const double log_maximum_weight = maximum_weight > 0 ? std::log(maximum_weight) : 0;
      const long double log_range = minimum_weight > 0 && maximum_weight > minimum_weight ?
        std::log(static_cast<long double>(maximum_weight) - minimum_weight) -
          std::log(minimum_weight) - log_maximum_weight : 0.0L;
      const long double log_budget = interpolate && reference_peak > 0 ?
        std::log(interpolation_tolerance) + std::log(reference_peak) - std::log(64.0L) : 0.0L;

      if (row_active_phack[static_cast<size_t>(s)]) {
        Rf_error("Latent z densities are unavailable for active p-hacking kernels.");
      }
      for (int k = 0; k < K; ++k) {
        const int cell = s + S * k;
        // One implementation serves original grid evaluation and Cheb nodes.
        const auto inverse_integral = [&](double value, const double *rule_nodes,
            const double *rule_logs, const long double *rule_weights, int count,
            bool logarithm) -> long double {
          const double location = selnorm_fma(candidate_fraction[cell], mean_p[cell],
            context_fraction[cell] * (value * sei_p[k]));
          for (int j = 0; j < count; ++j) {
            node_means[j] = selnorm_fma(conditional_sd[cell], rule_nodes[j], location);
            normalizers[j] = 1.0L;
          }
          const bool ordinary = cpp_selnorm_step_normalizer_product(
            node_means.data(), count, sd_p[cell], sei_p[k], row_weights.data(),
            data, normalizers.data());
          long double inverse = 0.0L;
          if (ordinary) {
            for (int j = 0; j < count; ++j) inverse += rule_weights[j] / normalizers[j];
          }
          if (ordinary && inverse > 0.0L && std::isfinite(inverse)) {
            // Cheb values retain the already accumulated positive integral;
            // direct grid evaluation retains its original double log result.
            return logarithm ? static_cast<double>(std::log(inverse)) : inverse;
          }
          double value_log = -std::numeric_limits<double>::infinity();
          for (int j = 0; j < count; ++j) {
            const double normalizer = cpp_selnorm_step_log_norm(node_means[j],
              sd_p[cell], sei_p[k], row_weights.data(), data, 1, false);
            if (!std::isfinite(normalizer)) {
              Rf_error("The retained-context selection normalizer is unavailable.");
            }
            value_log = logspace_add(value_log, rule_logs[j] - normalizer);
          }
          return logarithm ? static_cast<long double>(value_log) :
            std::exp(static_cast<long double>(value_log));
        };
        const bool population_nodes = selected && latent_p[cell] > 0 &&
          latent_p[cell] <= sd_p[cell];
        if (interpolate && selected && latent_p[cell] > 0.0 &&
            minimum_weight > 0.0 && maximum_weight > minimum_weight && reference_peak > 0) {
          const long double log_a = std::log(static_cast<long double>(latent_p[cell])) -
            std::log(static_cast<long double>(total_sd[cell])) + std::log(sei_p[k]) - std::log(sd_p[cell]);
          const long double log_density_peak = log_scale[cell] - half_log_two_pi + log_maximum_weight;
          const ZcurveChebyshevGrid *chosen = nullptr;
          long double remainder = 0.0L;
          for (ZcurveChebyshevGrid const &grid : chebyshev) {
            // Computational routing only; the same remainder/error gate applies.
            if (population_nodes && grid.count > 16) break;
            if (!grid.valid) continue;
            const long double log_remainder = log_range + grid.count * (log_a + log_interval_quarter) -
              grid.half_log_factorial;
            if (std::isfinite(log_remainder) && log_remainder + log_density_peak <= log_budget) {
              chosen = &grid;
              remainder = std::exp(log_remainder);
              break;
            }
          }
          if (chosen != nullptr) {
            std::vector<long double> H(chosen->count), node_error(chosen->count);
            bool valid = std::isfinite(remainder);
            for (int j = 0; j < chosen->count && valid; ++j) {
              H[j] = inverse_integral(chosen->nodes[j], nodes_p, log_weights_p,
                quadrature_weights.data(), n_nodes, false);
              if (coarse_count > 0) {
                const long double coarse = inverse_integral(chosen->nodes[j],
                  coarse_nodes, coarse_logs, coarse_weights.data(), coarse_count, false);
                node_error[j] = std::fabs(H[j] - coarse);
              }
              valid = H[j] > 0.0L && std::isfinite(H[j]) && std::isfinite(node_error[j]);
            }
            ZcurveDoubleInterpolation double_rule(*chosen, H, node_error);
            const long double rounding_budget = std::exp(log_budget);
            const double unit_log_coefficient = 0.0;
            const bool gaussian_projected = valid && projection.density(mean_p + cell, 1,
              total_sd[cell], sei_p[k], &unit_log_coefficient, cell_gaussian.data());
            for (int ell = 0; ell < L && valid; ++ell) {
              const long double gaussian = gaussian_projected ? cell_gaussian[ell] :
                std::exp(static_cast<long double>(cpp_selnorm_affine_normal_lpdf_log_scale(
                  z_p[ell], sei_p[k], mean_p[cell], total_sd[cell], log_scale[cell])));
              long double value = 0.0L, error = remainder, rounding_error = 0.0L, node_bound = 0.0L;
              bool ordinary = double_rule.evaluate(ell, &value, &node_bound, &rounding_error);
              ordinary = ordinary && value > 0 &&
                gaussian * maximum_weight * rounding_error <= rounding_budget;
              if (ordinary) {
                error += node_bound + rounding_error;
              } else {
                value = 0.0L;
                for (int j = 0; j < chosen->count; ++j) {
                  const long double coefficient = chosen->cardinal[static_cast<std::size_t>(ell) * chosen->count + j];
                  value += coefficient * H[j];
                  error += std::fabs(coefficient) * node_error[j];
                }
              }
              if (!(value > 0.0L) || !std::isfinite(value) || !std::isfinite(error) ||
                  value + error < 1.0L / maximum_weight || value - error > 1.0L / minimum_weight) {
                valid = false;
                break;
              }
              cell_density[ell] = gaussian * value;
              cell_error[ell] = gaussian * error;
              valid = std::isfinite(cell_density[ell]) && std::isfinite(cell_error[ell]);
            }
            if (valid) {
              for (int ell = 0; ell < L; ++ell) {
                density[ell] += cell_density[ell];
                interpolation_error[ell] += cell_error[ell];
              }
              continue;
            }
          }
        }
        if (population_nodes) {
          const auto forward_curve = [&](const double *rule_nodes, const double *rule_logs,
              int count, std::vector<long double> &curve) {
            for (int j = 0; j < count; ++j) {
              node_means[j] = selnorm_fma(latent_p[cell], rule_nodes[j], mean_p[cell]);
              normalizers[j] = 1.0L;
            }
            const bool ordinary = cpp_selnorm_step_normalizer_product(node_means.data(),
              count, sd_p[cell], sei_p[k], row_weights.data(), data, normalizers.data());
            for (int j = 0; j < count; ++j) {
              const double normalizer = ordinary ? static_cast<double>(std::log(normalizers[j])) :
                cpp_selnorm_step_log_norm(node_means[j], sd_p[cell], sei_p[k], row_weights.data(), data, 1, false);
              if (!std::isfinite(normalizer)) Rf_error("The retained-context selection normalizer is unavailable.");
              log_coefficients[j] = rule_logs[j] - normalizer;
            }
            if (projection.density(node_means.data(), count, sd_p[cell], sei_p[k],
                log_coefficients.data(), curve.data())) return;
            const double scale = std::log(sei_p[k]) - std::log(sd_p[cell]);
            for (int ell = 0; ell < L; ++ell) {
              double value = -std::numeric_limits<double>::infinity();
              for (int j = 0; j < count; ++j) {
                value = logspace_add(value, log_coefficients[j] +
                  cpp_selnorm_affine_normal_lpdf_log_scale(z_p[ell], sei_p[k], node_means[j], sd_p[cell], scale));
              }
              curve[ell] = std::exp(static_cast<long double>(value));
            }
          };
          forward_curve(nodes_p, log_weights_p, n_nodes, cell_density);
          if (requested && coarse_count > 0) {
            forward_curve(coarse_nodes, coarse_logs, coarse_count, cell_coarse);
          }
          for (int ell = 0; ell < L; ++ell) {
            density[ell] += cell_density[ell];
            if (requested) {
              noninterpolated[ell] += cell_density[ell];
              if (coarse_count > 0) coarse_noninterpolated[ell] += cell_coarse[ell];
            }
          }
          continue;
        }

        for (int ell = 0; ell < L; ++ell) {
          const double gaussian = cpp_selnorm_affine_normal_lpdf_log_scale(
            z_p[ell], sei_p[k], mean_p[cell], total_sd[cell], log_scale[cell]);
          double fine_inverse = selected ? -log_norm[cell] : 0;
          double coarse_inverse = fine_inverse;
          if (selected && latent_p[cell] > 0) {
            fine_inverse = inverse_integral(z_p[ell], nodes_p, log_weights_p,
              quadrature_weights.data(), n_nodes, true);
            if (requested && coarse_count > 0) {
              coarse_inverse = inverse_integral(z_p[ell], coarse_nodes, coarse_logs,
                coarse_weights.data(), coarse_count, true);
            }
          }
          const long double value = std::exp(static_cast<long double>(gaussian + fine_inverse));
          density[ell] += value;
          if (requested) {
            noninterpolated[ell] += value;
            if (coarse_count > 0) {
              coarse_noninterpolated[ell] += std::exp(static_cast<long double>(gaussian + coarse_inverse));
            }
          }
        }
      }
      for (int ell = 0; ell < L; ++ell) {
        const double weight = selected && (!extrapolate_b || paired) ?
          row_weights[z_bin[static_cast<size_t>(ell)] - 1] : 1.0;
        const long double value = density[ell] / K;
        if (requested) {
          const long double error = (std::fabs(noninterpolated[ell] - coarse_noninterpolated[ell]) +
            interpolation_error[ell]) / K;
          const double fitted_value = weight > 0 ? static_cast<double>(weight * value) : 0.0;
          const double extrapolated_value = static_cast<double>(value);
          const long double unavailable = std::numeric_limits<long double>::infinity();
          if (!std::isfinite(error) || error < 0 || !std::isfinite(extrapolated_value) || extrapolated_value < 0) {
            extrapolated_absolute_error = unavailable;
          } else {
            extrapolated_absolute_error = std::max(extrapolated_absolute_error, error);
            extrapolated_peak = std::max(extrapolated_peak, static_cast<long double>(extrapolated_value));
          }
          if (weight > 0) {
            if (!std::isfinite(error) || error < 0 || !std::isfinite(fitted_value) || fitted_value < 0) {
              fitted_absolute_error = unavailable;
            } else {
              fitted_absolute_error = std::max(fitted_absolute_error, weight * error);
              fitted_peak = std::max(fitted_peak, static_cast<long double>(fitted_value));
            }
          }
          if (weight > 0) zero_weight_display = false;
        }
        if (paired) {
          extrapolated_p[s + S * ell] = static_cast<double>(value);
          fitted_p[s + S * ell] = weight > 0 ? static_cast<double>(value * weight) : 0;
        } else {
          out_p[s + S * ell] = weight > 0 ? static_cast<double>(value * weight) : 0;
        }
      }
      if (requested) {
        const auto relative_error = [](long double error, long double peak, bool structural_zero) {
          if (structural_zero && error == 0 && peak == 0) return 0.0;
          if (!std::isfinite(error) || !std::isfinite(peak) || !(peak > error)) {
            return std::numeric_limits<double>::infinity();
          }
          return ZcurveDoubleInterpolation::upward(error / (peak - error));
        };
        REAL(fitted_error)[s] = relative_error(fitted_absolute_error, fitted_peak, zero_weight_display);
        REAL(extrapolated_error)[s] = relative_error(extrapolated_absolute_error, extrapolated_peak, false);
      }
    }
    if (requested) {
      if (paired) {
        Rf_setAttrib(fitted, Rf_install("relative_integration_error"), fitted_error);
        Rf_setAttrib(extrapolated, Rf_install("relative_integration_error"), extrapolated_error);
      } else {
        Rf_setAttrib(out, Rf_install("relative_integration_error"),
          extrapolate_b ? extrapolated_error : fitted_error);
      }
    }
    UNPROTECT(n_protect);
    return out;
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
            fitted_sum += cpp_selnorm_affine_normal_pdf_log_scale(
              z_p[ell], sei_p[k], mean_p[cell], sd_p[cell],
              log_scale[static_cast<size_t>(cell)]
            );
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
          dens_sum += cpp_selnorm_affine_normal_pdf_log_scale(
            z_p[ell], sei_p[k], mean_p[cell], sd_p[cell],
            log_scale[static_cast<size_t>(cell)]
          );
        }
      }

      out_p[s + S * ell] = dens_sum / K;
    }
  }

  UNPROTECT(n_protect);
  return out;
}
