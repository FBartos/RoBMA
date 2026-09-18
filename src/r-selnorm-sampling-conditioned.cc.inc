// Independent row deletion integrates inverse acceptance under E | Y. Inputs
// are already on the positive-selection scale. An unconverged cell stays NA so
// the R caller retains its existing adaptive-integration fallback and messages.
extern "C" SEXP RoBMA_selnorm_sampling_deletion_loglik_batch(
    SEXP yi, SEXP means, SEXP sampling_variances, SEXP integrated_variances,
    SEXP total_variances, SEXP sei, SEXP omega, SEXP z_lower, SEXP z_upper,
    SEXP obs_bin, SEXP kernel_mode, SEXP telescope_probabilities,
    SEXP nodes, SEXP log_weights, SEXP orders, SEXP relative_tolerance)
{
  const SEXP numeric_inputs[] = {yi, means, sampling_variances, integrated_variances,
    total_variances, sei, omega, z_lower, z_upper, nodes, log_weights, relative_tolerance};
  const char *numeric_labels[] = {"yi", "means", "sampling_variances", "integrated_variances",
    "total_variances", "sei", "omega", "z_lower", "z_upper", "nodes", "log_weights", "relative_tolerance"};
  for (int i = 0; i < 12; ++i) check_real(numeric_inputs[i], numeric_labels[i]);
  check_integer(obs_bin, "obs_bin");
  check_integer(kernel_mode, "kernel_mode");
  check_integer(orders, "orders");
  check_logical(telescope_probabilities, "telescope_probabilities");
  const int S = matrix_nrow(means, "means"), K = matrix_ncol(means, "means");
  const int B = matrix_ncol(omega, "omega"), rule_count = Rf_length(orders);
  const int mode_length = Rf_length(kernel_mode);
  if (S < 1 || K < 1 || B < 1 || Rf_length(yi) != K ||
      Rf_length(sampling_variances) != K || Rf_length(sei) != K ||
      Rf_length(obs_bin) != K || Rf_length(z_lower) != B || Rf_length(z_upper) != B ||
      (mode_length != 1 && mode_length != S) || rule_count < 2 ||
      Rf_length(relative_tolerance) != 1 || !(REAL(relative_tolerance)[0] > 0.0) ||
      !std::isfinite(REAL(relative_tolerance)[0]) ||
      LOGICAL(telescope_probabilities)[0] == NA_LOGICAL) {
    Rf_error("Sampling-conditioned scalar deletion inputs have invalid dimensions or controls.");
  }
  validate_matrix_dims_named(integrated_variances, "integrated_variances", S, K, "means");
  validate_matrix_dims_named(total_variances, "total_variances", S, K, "means");
  validate_matrix_dims_named(omega, "omega", S, B, "means");
  validate_omega_matrix(REAL(omega), S, B);
  for (R_xlen_t i = 0; i < XLENGTH(omega); ++i) {
    if (!(REAL(omega)[i] > 0.0)) {
      Rf_error("Sampling-conditioned scalar deletion requires strictly positive selection weights.");
    }
  }
  if (!selnorm_is_descending_step_partition(REAL(z_lower), REAL(z_upper), B)) {
    Rf_error("Selection bounds do not define a descending step partition.");
  }
  for (int k = 0; k < K; ++k) {
    if (!std::isfinite(REAL(yi)[k]) || !std::isfinite(REAL(sei)[k]) ||
        !(REAL(sei)[k] > 0.0) || !std::isfinite(REAL(sampling_variances)[k]) ||
        !(REAL(sampling_variances)[k] > 0.0) ||
        INTEGER(obs_bin)[k] < 1 || INTEGER(obs_bin)[k] > B) {
      Rf_error("Sampling-conditioned scalar deletion outcome inputs are invalid.");
    }
  }
  for (R_xlen_t i = 0; i < XLENGTH(means); ++i) {
    if (!std::isfinite(REAL(means)[i]) || !std::isfinite(REAL(integrated_variances)[i]) ||
        REAL(integrated_variances)[i] < 0.0 || !std::isfinite(REAL(total_variances)[i]) ||
        !(REAL(total_variances)[i] > 0.0)) {
      Rf_error("Sampling-conditioned scalar deletion Gaussian inputs are invalid.");
    }
  }
  int total_nodes = 0, maximum_order = 0;
  for (int r = 0; r < rule_count; ++r) {
    const int order = INTEGER(orders)[r];
    if (order <= maximum_order || total_nodes > Rf_length(nodes) - order) {
      Rf_error("Sampling-conditioned scalar deletion quadrature orders are invalid.");
    }
    maximum_order = order;
    total_nodes += order;
  }
  if (Rf_length(nodes) != total_nodes || Rf_length(log_weights) != total_nodes) {
    Rf_error("Sampling-conditioned scalar deletion quadrature arrays do not match their orders.");
  }
  for (int i = 0; i < total_nodes; ++i) {
    if (!std::isfinite(REAL(nodes)[i]) || !std::isfinite(REAL(log_weights)[i])) {
      Rf_error("Sampling-conditioned scalar deletion quadrature must contain finite values.");
    }
  }
  for (int s = 0; s < S; ++s) {
    const int mode = scalar_or_row_int(INTEGER(kernel_mode), mode_length, s, "kernel_mode");
    if (mode != SELKERNEL_NORMAL && mode != SELKERNEL_STEP) {
      Rf_error("Sampling-conditioned scalar deletion requires a normal or step kernel.");
    }
  }

  SelNormKernelData data = {};
  data.n_bins = B;
  data.effect_sign = 1;
  data.z_lower = REAL(z_lower);
  data.z_upper = REAL(z_upper);
  data.trusted_step_partition = true;
  data.telescope_probabilities = LOGICAL(telescope_probabilities)[0] == TRUE;
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  std::fill_n(REAL(out), static_cast<R_xlen_t>(S) * K, NA_REAL);
  std::vector<double> weights(B), node_means(maximum_order);
  std::vector<long double> normalizers(maximum_order);
  // The rule weights depend on nothing but the rule; see selnorm_rule_weights.
  const SelNormRuleWeights quadrature_weights = selnorm_rule_weights(
    REAL(log_weights), total_nodes);
  const double tolerance = REAL(relative_tolerance)[0];
  std::string failure;
  try {
    for (int s = 0; s < S; ++s) {
      const int mode = scalar_or_row_int(INTEGER(kernel_mode), mode_length, s, "kernel_mode");
      bool constant_weights = true;
      for (int b = 0; b < B; ++b) {
        weights[b] = REAL(omega)[s + S * b];
        constant_weights = constant_weights && weights[b] == weights[0];
      }
      for (int k = 0; k < K; ++k) {
        const int cell = s + S * k;
        const double mean = REAL(means)[cell];
        const double variance = REAL(integrated_variances)[cell];
        if (variance == 0.0 || mode == SELKERNEL_NORMAL || constant_weights) {
          REAL(out)[cell] = Rf_dnorm4(REAL(yi)[k], mean, std::sqrt(REAL(total_variances)[cell]), 1);
          continue;
        }
        const double sampling = REAL(sampling_variances)[k];
        const double total = sampling + variance;
        const double conditional_mean = mean + sampling / total * (REAL(yi)[k] - mean);
        const double conditional_sd = std::sqrt(sampling * variance / total);
        const double sd = std::sqrt(variance);
        const double gaussian = Rf_dnorm4(REAL(yi)[k], mean, std::sqrt(total), 1);
        const double log_weight = std::log(weights[INTEGER(obs_bin)[k] - 1]);
        double previous = NA_REAL;
        int offset = 0;
        for (int r = 0; r < rule_count; ++r) {
          const int order = INTEGER(orders)[r];
          for (int j = 0; j < order; ++j) {
            node_means[j] = conditional_mean + conditional_sd * REAL(nodes)[offset + j];
            normalizers[j] = 1.0L;
          }
          const bool direct = cpp_selnorm_step_normalizer_product(
            node_means.data(), order, sd, REAL(sei)[k], weights.data(), data,
            normalizers.data());
          long double inverse = 0.0L;
          if (direct) {
            for (int j = 0; j < order; ++j) {
              inverse += (*quadrature_weights)[offset + j] / normalizers[j];
            }
          }
          double log_inverse;
          if (direct && inverse > 0.0L && std::isfinite(inverse)) {
            log_inverse = static_cast<double>(std::log(inverse));
          } else {
            log_inverse = -std::numeric_limits<double>::infinity();
            for (int j = 0; j < order; ++j) {
              const double normalizer = cpp_selnorm_step_log_norm(
                node_means[j], sd, REAL(sei)[k], weights.data(), data, 1, false);
              if (!std::isfinite(normalizer)) {
                log_inverse = NA_REAL;
                break;
              }
              log_inverse = logspace_add(log_inverse, REAL(log_weights)[offset + j] - normalizer);
            }
          }
          const double current = gaussian + log_weight + log_inverse;
          const double change = current == previous ? 0.0 : std::abs(std::expm1(current - previous));
          if (r > 0 && !std::isnan(current) && current != std::numeric_limits<double>::infinity() &&
              std::isfinite(change) && change <= tolerance) {
            REAL(out)[cell] = current;
            break;
          }
          previous = current;
          offset += order;
        }
      }
      R_CheckUserInterrupt();
    }
  } catch (const std::exception &error) {
    failure = error.what();
  }
  if (!failure.empty()) Rf_error("%s", failure.c_str());
  UNPROTECT(1);
  return out;
}


static SelNormConditionedQuadrature sampling_conditioned_quadrature(
    SEXP cluster_nodes, SEXP cluster_log_weights, SEXP cluster_orders,
    SEXP factor_nodes, SEXP factor_log_weights, SEXP factor_orders)
{
  const SEXP inputs[] = {cluster_nodes, cluster_log_weights, cluster_orders,
    factor_nodes, factor_log_weights, factor_orders};
  const char *labels[] = {"cluster_nodes", "cluster_log_weights", "cluster_orders",
    "factor_nodes", "factor_log_weights", "factor_orders"};
  for (int i = 0; i < 6; ++i) {
    check_real(inputs[i], labels[i]);
    for (R_xlen_t j = 0; j < XLENGTH(inputs[i]); ++j) {
      if (!std::isfinite(REAL(inputs[i])[j])) Rf_error("'%s' must contain finite values.", labels[i]);
    }
  }
  for (int offset : {0, 3}) {
    const int count = Rf_length(inputs[offset + 2]);
    if (count < 3) Rf_error("Conditional selection quadrature must contain at least three rules.");
    int total = 0, previous = 0;
    for (int i = 0; i < count; ++i) {
      const double order = REAL(inputs[offset + 2])[i];
      if (order <= previous || order != std::floor(order)) {
        Rf_error("Conditional selection quadrature orders must be increasing positive integers.");
      }
      previous = static_cast<int>(order);
      total += previous;
    }
    if (Rf_length(inputs[offset]) != total || Rf_length(inputs[offset + 1]) != total) {
      Rf_error("Conditional selection quadrature arrays do not match their orders.");
    }
  }
  return {REAL(cluster_nodes), REAL(cluster_log_weights), REAL(cluster_orders),
    Rf_length(cluster_orders), REAL(factor_nodes), REAL(factor_log_weights),
    REAL(factor_orders), Rf_length(factor_orders)};
}

static SEXP selnorm_conditioned_normalizer_or_projection_batch(
    SEXP means, SEXP diagonal, SEXP loading, SEXP rank, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP sign, SEXP kernel_mode,
    SEXP telescope_probabilities, SEXP vector_rule, SEXP group_index, SEXP qmc,
    SEXP initial_points, SEXP max_points, SEXP scrambles, SEXP relative_tolerance,
    SEXP cluster_nodes, SEXP cluster_log_weights, SEXP cluster_orders,
    SEXP factor_nodes, SEXP factor_log_weights, SEXP factor_orders,
    SEXP z, SEXP kind)
{
  const SelNormConditionedQuadrature quadrature = sampling_conditioned_quadrature(
    cluster_nodes, cluster_log_weights, cluster_orders, factor_nodes,
    factor_log_weights, factor_orders);
  check_real(means, "means");
  check_real(diagonal, "diagonal");
  check_real(loading, "loading");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_real(qmc, "qmc");
  check_real(relative_tolerance, "relative_tolerance");
  check_integer(rank, "rank");
  check_integer(sign, "sign");
  check_integer(kernel_mode, "kernel_mode");
  check_integer(group_index, "group_index");
  check_integer(initial_points, "initial_points");
  check_integer(max_points, "max_points");
  check_integer(scrambles, "scrambles");
  check_logical(telescope_probabilities, "telescope_probabilities");
  const int S = matrix_nrow(means, "means"), K = matrix_ncol(means, "means");
  const int B = matrix_ncol(omega, "omega");
  if (S < 1 || K < 1 || B < 1 || Rf_length(rank) != 1 || INTEGER(rank)[0] < 0 ||
      Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1) ||
      Rf_length(initial_points) != 1 || INTEGER(initial_points)[0] < 1 ||
      Rf_length(max_points) != 1 || INTEGER(max_points)[0] < INTEGER(initial_points)[0] ||
      Rf_length(scrambles) != 1 || INTEGER(scrambles)[0] < 2 ||
      Rf_length(relative_tolerance) != 1 || !(REAL(relative_tolerance)[0] > 0.0) ||
      !std::isfinite(REAL(relative_tolerance)[0]) ||
      LOGICAL(telescope_probabilities)[0] == NA_LOGICAL) {
    Rf_error("Conditional sampling selection controls are invalid.");
  }
  const int R = INTEGER(rank)[0], maximum = INTEGER(max_points)[0], shifts = INTEGER(scrambles)[0];
  const int mode_length = Rf_length(kernel_mode), rule_length = validate_vector_rule(vector_rule, S);
  const bool project = z != R_NilValue;
  int projection_kind = 0, size = 0;
  if (project) {
    check_real(z, "z");
    check_integer(kind, "kind");
    if (Rf_length(kind) != 1 || INTEGER(kind)[0] < 0 || INTEGER(kind)[0] > 2 || Rf_length(z) < 1) {
      Rf_error("Factor projection requires a nonempty grid and density, CDF, or two-sided-tail kind.");
    }
    projection_kind = INTEGER(kind)[0];
    size = Rf_length(z);
    for (int index = 0; index < size; ++index) {
      if (std::isnan(REAL(z)[index]) ||
          (projection_kind != 1 && !std::isfinite(REAL(z)[index])) ||
          (projection_kind == 2 && REAL(z)[index] < 0.0)) {
        Rf_error("Factor projection grid values are invalid for the requested kind.");
      }
    }
  }
  if (Rf_length(sei) != K || Rf_length(group_index) != K || Rf_length(z_lower) != B ||
      Rf_length(z_upper) != B || (mode_length != 1 && mode_length != S) ||
      (!(!project && R == 0 && XLENGTH(qmc) == 1) &&
       XLENGTH(qmc) != static_cast<R_xlen_t>(shifts) * maximum * std::max(2 * K, project ? K + R : R))) {
    Rf_error("Conditional sampling selection inputs have invalid dimensions.");
  }
  validate_matrix_dims_named(diagonal, "diagonal", S, K, "means");
  validate_matrix_dims_named(loading, "loading", S, K * R, "means");
  validate_matrix_dims_named(omega, "omega", S, B, "means");
  validate_omega_matrix(REAL(omega), S, B);
  if (!selnorm_is_descending_step_partition(REAL(z_lower), REAL(z_upper), B)) {
    Rf_error("Selection bounds do not define a descending step partition.");
  }
  for (int row = 0; row < K; ++row) {
    if (!(REAL(sei)[row] > 0.0) || !std::isfinite(REAL(sei)[row]) || INTEGER(group_index)[row] < 1) {
      Rf_error("Conditional sampling selection outcome inputs are invalid.");
    }
  }
  for (R_xlen_t i = 0; i < XLENGTH(means); ++i) {
    if (!std::isfinite(REAL(means)[i]) || !std::isfinite(REAL(diagonal)[i]) || REAL(diagonal)[i] < 0.0) {
      Rf_error("Conditional sampling selection Gaussian inputs must be finite, with nonnegative candidate diagonal variances.");
    }
  }
  for (R_xlen_t i = 0; i < XLENGTH(loading); ++i) {
    if (!std::isfinite(REAL(loading)[i])) Rf_error("'loading' must contain finite values.");
  }
  for (R_xlen_t i = 0; i < XLENGTH(qmc); ++i) {
    if (!(REAL(qmc)[i] > 0.0 && REAL(qmc)[i] < 1.0)) {
      Rf_error("'qmc' must contain finite values strictly between zero and one.");
    }
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  const char *labels[] = {project ? "log_density" : "log_normalizer", "relative_mcse",
    project ? "density" : "relative_change"};
  for (int i = 0; i < 3; ++i) {
    SET_VECTOR_ELT(out, i, project && i == 2 ? Rf_allocMatrix(REALSXP, S, size) : Rf_allocVector(REALSXP, S));
    SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  // Kernel modes are validated before the row loop so parallel workers never
  // raise errors from a worker thread.
  for (int s = 0; s < S; ++s) {
    const int mode = scalar_or_row_int(INTEGER(kernel_mode), mode_length, s, "kernel_mode");
    if (mode != SELKERNEL_NORMAL && mode != SELKERNEL_STEP) {
      Rf_error("Conditional sampling selection requires a normal or step kernel.");
    }
  }
  const int threads = robma_batch_threads(S);
  robma_parallel::for_rows(S, threads, [&](int s) {
    std::vector<double> row_m(static_cast<std::size_t>(K));
    std::vector<double> row_d(static_cast<std::size_t>(K));
    std::vector<double> row_l(static_cast<std::size_t>(K * R));
    std::vector<double> row_w(static_cast<std::size_t>(B));
    std::vector<double> row_projected(static_cast<std::size_t>(size));
    const int mode = scalar_or_row_int(INTEGER(kernel_mode), mode_length, s, "kernel_mode");
    for (int row = 0; row < K; ++row) {
      row_m[row] = REAL(means)[s + S * row];
      row_d[row] = REAL(diagonal)[s + S * row];
    }
    for (int i = 0; i < K * R; ++i) row_l[i] = REAL(loading)[s + S * i];
    for (int bin = 0; bin < B; ++bin) row_w[bin] = REAL(omega)[s + S * bin];
    if (project) {
      double normalizer = 0.0;
      cpp_selnorm_factor_projection(row_m.data(), row_d.data(), row_l.data(), K, R,
        REAL(sei), row_w.data(), B, REAL(z_lower), REAL(z_upper), INTEGER(sign)[0],
        mode, LOGICAL(telescope_probabilities)[0],
        scalar_or_row_int(INTEGER(vector_rule), rule_length, s, "vector_rule"),
        INTEGER(group_index), REAL(qmc), INTEGER(initial_points)[0], maximum,
        shifts, REAL(relative_tolerance)[0], quadrature, REAL(z), size,
        projection_kind, row_projected.data(), &normalizer, REAL(VECTOR_ELT(out, 1)) + s);
      REAL(VECTOR_ELT(out, 0))[s] = -normalizer;
      for (int index = 0; index < size; ++index) REAL(VECTOR_ELT(out, 2))[s + S * index] = row_projected[index];
    } else REAL(VECTOR_ELT(out, 0))[s] = cpp_selnorm_conditioned_log_normalizer(row_m.data(),
        row_d.data(), row_l.data(), K, R, REAL(sei), row_w.data(), B, REAL(z_lower), REAL(z_upper),
        INTEGER(sign)[0], mode, LOGICAL(telescope_probabilities)[0],
        scalar_or_row_int(INTEGER(vector_rule), rule_length, s, "vector_rule"),
        INTEGER(group_index), REAL(qmc), INTEGER(initial_points)[0], maximum,
        shifts, REAL(relative_tolerance)[0], REAL(VECTOR_ELT(out, 1)) + s,
        REAL(VECTOR_ELT(out, 2)) + s, quadrature);
  });
  UNPROTECT(2);
  return out;
}

extern "C" SEXP RoBMA_selnorm_conditioned_normalizer_batch(
    SEXP means, SEXP diagonal, SEXP loading, SEXP rank, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP sign, SEXP kernel_mode,
    SEXP telescope_probabilities, SEXP vector_rule, SEXP group_index, SEXP qmc,
    SEXP initial_points, SEXP max_points, SEXP scrambles, SEXP relative_tolerance,
    SEXP cluster_nodes, SEXP cluster_log_weights, SEXP cluster_orders,
    SEXP factor_nodes, SEXP factor_log_weights, SEXP factor_orders)
{
  return selnorm_conditioned_normalizer_or_projection_batch(means, diagonal,
    loading, rank, sei, omega, z_lower, z_upper, sign, kernel_mode,
    telescope_probabilities, vector_rule, group_index, qmc, initial_points,
    max_points, scrambles, relative_tolerance, cluster_nodes, cluster_log_weights,
    cluster_orders, factor_nodes, factor_log_weights, factor_orders,
    R_NilValue, R_NilValue);
}

extern "C" SEXP RoBMA_selnorm_factor_projection_batch(
    SEXP means, SEXP diagonal, SEXP loading, SEXP rank, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP sign, SEXP kernel_mode,
    SEXP telescope_probabilities, SEXP vector_rule, SEXP group_index, SEXP qmc,
    SEXP initial_points, SEXP max_points, SEXP scrambles, SEXP relative_tolerance,
    SEXP cluster_nodes, SEXP cluster_log_weights, SEXP cluster_orders,
    SEXP factor_nodes, SEXP factor_log_weights, SEXP factor_orders,
    SEXP z, SEXP kind)
{
  return selnorm_conditioned_normalizer_or_projection_batch(means, diagonal,
    loading, rank, sei, omega, z_lower, z_upper, sign, kernel_mode,
    telescope_probabilities, vector_rule, group_index, qmc, initial_points,
    max_points, scrambles, relative_tolerance, cluster_nodes, cluster_log_weights,
    cluster_orders, factor_nodes, factor_log_weights, factor_orders,
    z, kind);
}

extern "C" SEXP RoBMA_selnorm_sampling_conditioned_batch(
    SEXP yi, SEXP means, SEXP sampling_lower, SEXP diagonal, SEXP loading,
    SEXP rank, SEXP sampling_auxiliary, SEXP random_auxiliary, SEXP sei,
    SEXP omega, SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
    SEXP kernel_mode, SEXP telescope_probabilities, SEXP vector_rule,
    SEXP group_index, SEXP qmc, SEXP initial_points, SEXP max_points,
    SEXP scrambles, SEXP relative_tolerance,
    SEXP cluster_nodes, SEXP cluster_log_weights, SEXP cluster_orders,
    SEXP factor_nodes, SEXP factor_log_weights, SEXP factor_orders)
{
  const SelNormConditionedQuadrature quadrature = sampling_conditioned_quadrature(
    cluster_nodes, cluster_log_weights, cluster_orders, factor_nodes,
    factor_log_weights, factor_orders);
  check_real(yi, "yi");
  check_real(means, "means");
  check_real(sampling_lower, "sampling_lower");
  check_real(diagonal, "diagonal");
  check_real(loading, "loading");
  check_real(sampling_auxiliary, "sampling_auxiliary");
  check_real(random_auxiliary, "random_auxiliary");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_real(qmc, "qmc");
  check_real(relative_tolerance, "relative_tolerance");
  check_integer(rank, "rank");
  check_integer(obs_bin, "obs_bin");
  check_integer(sign, "sign");
  check_integer(kernel_mode, "kernel_mode");
  check_integer(group_index, "group_index");
  check_integer(initial_points, "initial_points");
  check_integer(max_points, "max_points");
  check_integer(scrambles, "scrambles");
  check_logical(telescope_probabilities, "telescope_probabilities");
  const int S = matrix_nrow(means, "means");
  const int K = matrix_ncol(means, "means");
  const int B = matrix_ncol(omega, "omega");
  const int lower_size = K * (K + 1) / 2;
  const bool fixed_covariance = Rf_getAttrib(sampling_lower, R_DimSymbol) == R_NilValue;
  if (S < 1 || K < 1 || B < 1 || Rf_length(rank) != 1 || INTEGER(rank)[0] < 0 ||
      Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1) ||
      Rf_length(initial_points) != 1 || INTEGER(initial_points)[0] < 1 ||
      Rf_length(max_points) != 1 || INTEGER(max_points)[0] < INTEGER(initial_points)[0] ||
      Rf_length(scrambles) != 1 || INTEGER(scrambles)[0] < 2 ||
      Rf_length(relative_tolerance) != 1 || !(REAL(relative_tolerance)[0] > 0.0) ||
      !std::isfinite(REAL(relative_tolerance)[0]) ||
      LOGICAL(telescope_probabilities)[0] == NA_LOGICAL) {
    Rf_error("Conditional sampling selection controls are invalid.");
  }
  const int R = INTEGER(rank)[0];
  const int points = INTEGER(initial_points)[0];
  const int maximum = INTEGER(max_points)[0];
  const int shifts = INTEGER(scrambles)[0];
  const int rule_length = validate_vector_rule(vector_rule, S);
  const int mode_length = Rf_length(kernel_mode);
  if (Rf_length(yi) != K || (fixed_covariance && Rf_length(sampling_lower) != lower_size) ||
      Rf_length(sei) != K || Rf_length(obs_bin) != K || Rf_length(group_index) != K ||
      Rf_length(z_lower) != B || Rf_length(z_upper) != B ||
      (mode_length != 1 && mode_length != S) ||
      (!(R == 0 && XLENGTH(qmc) == 1) &&
       XLENGTH(qmc) != static_cast<R_xlen_t>(shifts) * maximum * std::max(2 * K, R))) {
    Rf_error("Conditional sampling selection inputs have invalid dimensions.");
  }
  validate_matrix_dims_named(diagonal, "diagonal", S, K, "means");
  if (!fixed_covariance) {
    validate_matrix_dims_named(sampling_lower, "sampling_lower", S, lower_size, "means");
  }
  validate_matrix_dims_named(loading, "loading", S, K * R, "means");
  validate_matrix_dims_named(sampling_auxiliary, "sampling_auxiliary", S, K, "means");
  validate_matrix_dims_named(random_auxiliary, "random_auxiliary", S, K, "means");
  validate_matrix_dims_named(omega, "omega", S, B, "means");
  validate_omega_matrix(REAL(omega), S, B);
  if (!selnorm_is_descending_step_partition(REAL(z_lower), REAL(z_upper), B)) {
    Rf_error("Selection bounds do not define a descending step partition.");
  }
  for (int row = 0; row < K; ++row) {
    if (!std::isfinite(REAL(yi)[row]) || !std::isfinite(REAL(sei)[row]) ||
        !(REAL(sei)[row] > 0.0) || INTEGER(obs_bin)[row] < 1 ||
        INTEGER(obs_bin)[row] > B || INTEGER(group_index)[row] < 1) {
      Rf_error("Conditional sampling selection outcome inputs are invalid.");
    }
  }
  for (R_xlen_t i = 0; i < XLENGTH(sampling_lower); ++i) {
    if (!std::isfinite(REAL(sampling_lower)[i])) Rf_error("'sampling_lower' must contain finite values.");
  }
  for (R_xlen_t i = 0; i < XLENGTH(means); ++i) {
    if (!std::isfinite(REAL(means)[i]) || !std::isfinite(REAL(diagonal)[i]) ||
        REAL(diagonal)[i] < 0.0 || !std::isfinite(REAL(sampling_auxiliary)[i]) ||
        !std::isfinite(REAL(random_auxiliary)[i])) {
      Rf_error("Conditional sampling selection Gaussian inputs must be finite, with nonnegative candidate diagonal variances.");
    }
  }
  for (R_xlen_t i = 0; i < XLENGTH(loading); ++i) {
    if (!std::isfinite(REAL(loading)[i])) Rf_error("'loading' must contain finite values.");
  }
  for (R_xlen_t i = 0; i < XLENGTH(qmc); ++i) {
    if (!(REAL(qmc)[i] > 0.0 && REAL(qmc)[i] < 1.0)) {
      Rf_error("'qmc' must contain finite values strictly between zero and one.");
    }
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 6));
  const char *labels[] = {"log_lik", "e_star", "delta", "log_normalizer",
    "relative_mcse", "relative_change"};
  for (int i = 0; i < 6; ++i) {
    SET_VECTOR_ELT(out, i, (i == 1 || i == 2) ? Rf_allocMatrix(REALSXP, S, K) :
      Rf_allocVector(REALSXP, S));
    SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  std::vector<double> mean(K), d(K), l(K * R), e0(K), u0(K), weights(B), e(K), delta(K), retained(lower_size);
  std::string failure;
  try {
    for (int s = 0; s < S; ++s) {
      const int mode = scalar_or_row_int(INTEGER(kernel_mode), mode_length, s, "kernel_mode");
      if (mode != SELKERNEL_NORMAL && mode != SELKERNEL_STEP) {
        throw std::runtime_error("Conditional sampling selection requires a normal or step kernel.");
      }
      for (int row = 0; row < K; ++row) {
        mean[row] = REAL(means)[s + S * row];
        d[row] = REAL(diagonal)[s + S * row];
        e0[row] = REAL(sampling_auxiliary)[s + S * row];
        u0[row] = REAL(random_auxiliary)[s + S * row];
      }
      for (int i = 0; i < K * R; ++i) l[i] = REAL(loading)[s + S * i];
      for (int i = 0; i < lower_size; ++i) retained[i] = REAL(sampling_lower)[fixed_covariance ? i : s + S * i];
      for (int bin = 0; bin < B; ++bin) weights[bin] = REAL(omega)[s + S * bin];
      REAL(VECTOR_ELT(out, 0))[s] = cpp_selnorm_sampling_conditioned_lpdf(
        REAL(yi), mean.data(), retained.data(), d.data(), l.data(), K, R,
        e0.data(), u0.data(), REAL(sei), weights.data(), B, REAL(z_lower), REAL(z_upper),
        INTEGER(obs_bin), INTEGER(sign)[0], mode, LOGICAL(telescope_probabilities)[0],
        scalar_or_row_int(INTEGER(vector_rule), rule_length, s, "vector_rule"),
        INTEGER(group_index), REAL(qmc), points, maximum, shifts, REAL(relative_tolerance)[0],
        REAL(VECTOR_ELT(out, 4)) + s, REAL(VECTOR_ELT(out, 5)) + s,
        e.data(), delta.data(), REAL(VECTOR_ELT(out, 3)) + s, quadrature);
      for (int row = 0; row < K; ++row) {
        REAL(VECTOR_ELT(out, 1))[s + S * row] = e[row];
        REAL(VECTOR_ELT(out, 2))[s + S * row] = delta[row];
      }
      R_CheckUserInterrupt();
    }
  } catch (const std::exception &error) {
    failure = error.what();
  }
  if (!failure.empty()) Rf_error("%s", failure.c_str());
  UNPROTECT(2);
  return out;
}
