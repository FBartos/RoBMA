// Evaluation dispatchers try exact specialized representations from cheapest
// to most general. Returning false means "not applicable here" (or a recoverable
// factorization failure), never permission to approximate the covariance.

double plan_log_likelihood_values(
    CovariancePlan &plan,
    const double *mean_values,
    const std::vector<CovarianceFactor> &factors,
    const double *extra_values)
{
  for (int i = 0; i < plan.n; ++i) {
    if (!std::isfinite(extra_values[i]) || extra_values[i] < 0.0) {
      Rf_error("Known-V bridge extra variances must be finite and non-negative.");
    }
  }

  double log_likelihood = 0.0;
  for (BlockPlan &block : plan.blocks) {
    // Ordering affects runtime only. The final dense Cholesky call evaluates
    // the same block covariance when no specialization accepts the state.
    double block_value = 0.0;
    if (!markov_block_log_likelihood(
      plan,
      block,
      factors,
      mean_values,
      extra_values,
      &block_value
    ) && !fixed_known_group_block_log_likelihood(
      plan,
      block,
      factors,
      mean_values,
      extra_values,
      &block_value
    ) && !low_rank_block_log_likelihood(
      plan,
      block,
      factors,
      mean_values,
      extra_values,
      &block_value
    ) && !block_base_low_rank_log_likelihood(
      plan,
      block,
      factors,
      mean_values,
      extra_values,
      &block_value
    ) && !spectral_block_log_likelihood(
      plan,
      block,
      factors,
      mean_values,
      extra_values,
      &block_value
    ) && !root_dense_block_log_likelihood(
      plan,
      block,
      factors,
      mean_values,
      extra_values,
      &block_value
    )) {
      block_value = dense_block_log_likelihood(
        plan,
        block,
        factors,
        mean_values,
        extra_values
      );
    }
    log_likelihood += block_value;
  }
  return log_likelihood;
}

double plan_log_likelihood(CovariancePlan &plan, SEXP mean,
                           const std::vector<CovarianceFactor> &factors,
                           SEXP extra_variance)
{
  require_real_vector(mean, plan.n, "mean");
  require_real_vector(extra_variance, plan.n, "extra_variance");
  return plan_log_likelihood_values(
    plan,
    REAL(mean),
    factors,
    REAL(extra_variance)
  );
}

int require_batched_plan_inputs(
    const CovariancePlan &plan,
    SEXP means,
    SEXP random_covariance_states,
    SEXP extra_variances)
{
  if (TYPEOF(means) != REALSXP || !Rf_isMatrix(means) ||
      TYPEOF(extra_variances) != REALSXP || !Rf_isMatrix(extra_variances)) {
    Rf_error("Known-V batched means and extra variances must be numeric matrices.");
  }
  SEXP mean_dim = Rf_getAttrib(means, R_DimSymbol);
  SEXP extra_dim = Rf_getAttrib(extra_variances, R_DimSymbol);
  if (INTEGER(mean_dim)[0] != plan.n ||
      INTEGER(extra_dim)[0] != plan.n ||
      INTEGER(mean_dim)[1] != INTEGER(extra_dim)[1]) {
    Rf_error("Known-V batched matrices have inconsistent dimensions.");
  }
  if (TYPEOF(random_covariance_states) != VECSXP ||
      XLENGTH(random_covariance_states) != INTEGER(mean_dim)[1]) {
    Rf_error("Known-V batched random covariance states have inconsistent length.");
  }
  return INTEGER(mean_dim)[1];
}

void plan_location_quadratic_values(
    CovariancePlan &plan,
    const double *mean,
    const double *basis,
    const std::vector<CovarianceFactor> &factors,
    const double *extra,
    double *linear,
    double *quadratic)
{
  double linear_value = 0.0;
  double quadratic_value = 0.0;
  for (BlockPlan &block : plan.blocks) {
    const int size = static_cast<int>(block.index.size());
    double *covariance = plan.covariance_work.data();
    double *residual = plan.residual_work.data();
    fill_dense_block(
      plan,
      block,
      factors,
      mean,
      extra,
      covariance,
      residual
    );

    plan.conditional_rhs_work.resize(static_cast<size_t>(2 * size));
    double *right_hand_sides = plan.conditional_rhs_work.data();
    std::copy(residual, residual + size, right_hand_sides);
    for (int local = 0; local < size; ++local) {
      const int global = block.index[static_cast<size_t>(local)];
      if (!std::isfinite(basis[global])) {
        Rf_error("Known-V location bases must be finite.");
      }
      right_hand_sides[static_cast<size_t>(size + local)] = basis[global];
    }

    int info = 0;
    F77_CALL(dpotrf)("L", &size, covariance, &size, &info FCONE);
    if (info != 0) {
      Rf_error("Known-V integrated covariance is not positive definite.");
    }
    const int n_rhs = 2;
    F77_CALL(dpotrs)(
      "L",
      &size,
      &n_rhs,
      covariance,
      &size,
      right_hand_sides,
      &size,
      &info FCONE
    );
    if (info != 0) {
      Rf_error("Known-V location quadratic solve failed.");
    }
    for (int local = 0; local < size; ++local) {
      const int global = block.index[static_cast<size_t>(local)];
      linear_value += basis[global] *
        right_hand_sides[static_cast<size_t>(local)];
      quadratic_value += basis[global] *
        right_hand_sides[static_cast<size_t>(size + local)];
    }
  }

  *linear = linear_value;
  *quadratic = quadratic_value;
}

bool factor_dense_latent_system(
    CovariancePlan &plan,
    const std::vector<std::vector<int>> &active_columns,
    int size,
    int rank,
    const double *diagonal,
    const double *residual,
    const double *U)
{
  assemble_dense_latent_system(
    plan,
    active_columns,
    size,
    rank,
    diagonal,
    residual,
    U
  );
  int info = 0;
  F77_CALL(dpotrf)(
    "L",
    &rank,
    plan.small_matrix_work.data(),
    &rank,
    &info FCONE
  );
  return info == 0;
}

bool solve_latent_system(
    CovariancePlan &plan,
    BlockPlan &block,
    int rank,
    int n_rhs,
    bool sparse,
    double *rhs)
{
  if (!sparse) {
    int info = 0;
    F77_CALL(dpotrs)(
      "L",
      &rank,
      &n_rhs,
      plan.small_matrix_work.data(),
      &rank,
      rhs,
      &rank,
      &info FCONE
    );
    return info == 0;
  }

  cholmod_dense rhs_dense = {};
  rhs_dense.nrow = static_cast<size_t>(rank);
  rhs_dense.ncol = static_cast<size_t>(n_rhs);
  rhs_dense.nzmax = static_cast<size_t>(rank) *
    static_cast<size_t>(n_rhs);
  rhs_dense.d = static_cast<size_t>(rank);
  rhs_dense.x = rhs;
  rhs_dense.z = nullptr;
  rhs_dense.xtype = CHOLMOD_REAL;
  rhs_dense.dtype = CHOLMOD_DOUBLE;
  cholmod_dense *solution = M_cholmod_solve(
    CHOLMOD_A,
    block.sparse_factor,
    &rhs_dense,
    &plan.cholmod_common_state
  );
  if (solution == nullptr) {
    return false;
  }
  std::copy(
    static_cast<const double *>(solution->x),
    static_cast<const double *>(solution->x) + rhs_dense.nzmax,
    rhs
  );
  M_cholmod_free_dense(&solution, &plan.cholmod_common_state);
  return true;
}

bool finish_conditional_block(
    const BlockPlan &block,
    int size,
    int rank,
    const std::vector<std::vector<int>> &active_columns,
    const double *base_precision_diagonal,
    const double *base_precision_residual,
    const double *precision_factor,
    const double *solutions,
    const ConditionalOutput &output)
{
  // Q = B^-1 - C K^-1 C' and Qr = B^-1 r - C K^-1 h. Conditional log density,
  // residual, and variance are all derived from this same joint precision; none
  // is an independent marginal observation target.
  const double *latent_residual_solution = solutions;
  const bool residual_only = precision_residual_only(output);
  for (int observation = 0; observation < size; ++observation) {
    double residual_correction = 0.0;
    for (int column : active_columns[static_cast<size_t>(observation)]) {
      const double precision_value = precision_factor[static_cast<size_t>(
        observation + size * column
      )];
      residual_correction += precision_value *
        latent_residual_solution[static_cast<size_t>(column)];
    }
    const double precision_residual =
      base_precision_residual[static_cast<size_t>(observation)] -
      residual_correction;
    if (!std::isfinite(precision_residual)) {
      return false;
    }
    if (residual_only) {
      output.precision_residual[
        block.index[static_cast<size_t>(observation)]
      ] = precision_residual;
      continue;
    }

    const double *precision_solution = solutions +
      static_cast<size_t>((observation + 1) * rank);
    double diagonal_correction = 0.0;
    for (int column : active_columns[static_cast<size_t>(observation)]) {
      const double precision_value = precision_factor[static_cast<size_t>(
        observation + size * column
      )];
      diagonal_correction += precision_value *
        precision_solution[static_cast<size_t>(column)];
    }
    const double precision_diagonal =
      base_precision_diagonal[static_cast<size_t>(observation)] -
      diagonal_correction;
    const double cancellation_bound =
      static_cast<double>(rank + 2) *
      std::numeric_limits<double>::epsilon() *
      (std::abs(base_precision_diagonal[static_cast<size_t>(observation)]) +
       std::abs(diagonal_correction));
    if (!(precision_diagonal > cancellation_bound)) {
      // Do not clamp a cancellation-damaged diagonal. Decline this route so the
      // caller recomputes the same conditional density from the dense precision.
      return false;
    }
    if (!store_conditional_output(
        block.index[static_cast<size_t>(observation)],
        precision_diagonal,
        precision_residual,
        output
      )) {
      return false;
    }
  }
  return true;
}

bool diagonal_low_rank_conditional_values(
    CovariancePlan &plan,
    BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    const ConditionalOutput &output)
{
  if (!block.low_rank_eligible &&
      !(block.sparse_factor_eligible && !block.sparse_factor_block_base)) {
    return false;
  }
  const int size = static_cast<int>(block.index.size());
  const int rank = block.rank;
  double *diagonal = plan.diagonal_work.data();
  double *residual = plan.residual_work.data();
  for (int observation = 0; observation < size; ++observation) {
    const int global = block.index[static_cast<size_t>(observation)];
    const double variance = plan.sampling_covariance[static_cast<size_t>(
      global + plan.n * global
    )] + extra[global];
    if (!(variance > 0.0) || !std::isfinite(variance) ||
        !std::isfinite(mean[global])) {
      return false;
    }
    diagonal[static_cast<size_t>(observation)] = variance;
    residual[static_cast<size_t>(observation)] =
      plan.y[static_cast<size_t>(global)] - mean[global];
  }
  if (rank == 0) {
    const double two_pi = 6.283185307179586476925286766559;
    for (int observation = 0; observation < size; ++observation) {
      const double variance = diagonal[static_cast<size_t>(observation)];
      const double value = residual[static_cast<size_t>(observation)];
      const int global = block.index[static_cast<size_t>(observation)];
      if (output.log_density != nullptr) {
        output.log_density[global] = -0.5 * (
          std::log(two_pi * variance) + value * value / variance
        );
      }
      if (output.residual != nullptr) {
        output.residual[global] = value;
      }
      if (output.variance != nullptr) {
        output.variance[global] = variance;
      }
      if (output.precision_residual != nullptr) {
        output.precision_residual[global] = value / variance;
      }
    }
    return true;
  }

  double *U = fill_low_rank_factor(plan, block, factors);
  bool sparse = false;
  bool factored = false;
  if (block.sparse_factor_eligible && !block.sparse_factor_block_base) {
    factored = factor_sparse_latent_system(
      plan,
      block,
      size,
      rank,
      diagonal,
      residual,
      U
    );
    sparse = factored;
  }
  if (!factored && block.low_rank_eligible) {
    factored = factor_dense_latent_system(
      plan,
      block.active_columns,
      size,
      rank,
      diagonal,
      residual,
      U
    );
    sparse = false;
  }
  if (!factored) {
    return false;
  }

  const int n_rhs = precision_residual_only(output) ? 1 : size + 1;
  const size_t system_size = static_cast<size_t>(rank) *
    static_cast<size_t>(n_rhs);
  plan.conditional_rhs_work.resize(system_size);
  double *systems = plan.conditional_rhs_work.data();
  std::copy(
    plan.rhs_work.begin(),
    plan.rhs_work.begin() + rank,
    systems
  );
  if (n_rhs > 1) {
    std::fill(
      systems + rank,
      systems + system_size,
      0.0
    );
    for (int observation = 0; observation < size; ++observation) {
      const double inverse_variance = 1.0 /
        diagonal[static_cast<size_t>(observation)];
      for (int column : block.precision_active_columns[
          static_cast<size_t>(observation)
        ]) {
        systems[static_cast<size_t>(
          column + rank * (observation + 1)
        )] = U[static_cast<size_t>(observation + size * column)] *
          inverse_variance;
      }
    }
  }
  if (!solve_latent_system(
      plan,
      block,
      rank,
      n_rhs,
      sparse,
      systems
    )) {
    return false;
  }

  plan.conditional_base_work.resize(
    static_cast<size_t>(size) * static_cast<size_t>(rank)
  );
  double *precision_factor = plan.conditional_base_work.data();
  std::fill(
    precision_factor,
    precision_factor +
      static_cast<size_t>(size) * static_cast<size_t>(rank),
    0.0
  );
  for (int observation = 0; observation < size; ++observation) {
    const double inverse_variance = 1.0 /
      diagonal[static_cast<size_t>(observation)];
    for (int column : block.precision_active_columns[
        static_cast<size_t>(observation)
      ]) {
      precision_factor[static_cast<size_t>(
        observation + size * column
      )] = U[static_cast<size_t>(observation + size * column)] *
        inverse_variance;
    }
    residual[static_cast<size_t>(observation)] *= inverse_variance;
    diagonal[static_cast<size_t>(observation)] = inverse_variance;
  }
  return finish_conditional_block(
    block,
    size,
    rank,
    block.precision_active_columns,
    diagonal,
    residual,
    precision_factor,
    systems,
    output
  );
}

bool block_base_low_rank_conditional_values(
    CovariancePlan &plan,
    BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    const ConditionalOutput &output)
{
  if (!block.block_base_eligible &&
      !(block.sparse_factor_eligible && block.sparse_factor_block_base)) {
    return false;
  }
  const int size = static_cast<int>(block.index.size());
  const int rank = block.rank;
  double *base_precision_diagonal = plan.diagonal_work.data();
  double *residual = plan.residual_work.data();
  double *U = fill_low_rank_factor(plan, block, factors);
  plan.conditional_base_work.resize(
    static_cast<size_t>(size) * static_cast<size_t>(rank + 1)
  );
  double *base_precision = plan.conditional_base_work.data();
  for (int local = 0; local < size; ++local) {
    const int global = block.index[static_cast<size_t>(local)];
    if (!std::isfinite(mean[global])) {
      return false;
    }
    residual[static_cast<size_t>(local)] =
      plan.y[static_cast<size_t>(global)] - mean[global];
  }

  for (const SamplingBlock &sampling_block : block.sampling_blocks) {
    if (!prepare_block_base_sampling_block(
        plan,
        block,
        sampling_block,
        extra,
        size,
        rank,
        residual,
        U,
        nullptr,
        nullptr,
        base_precision_diagonal,
        base_precision
      )) {
      return false;
    }
  }

  plan.unit_diagonal_work.resize(static_cast<size_t>(size));
  std::fill(
    plan.unit_diagonal_work.begin(),
    plan.unit_diagonal_work.end(),
    1.0
  );
  const double *unit_diagonal = plan.unit_diagonal_work.data();
  bool sparse = false;
  bool factored = false;
  if (block.sparse_factor_eligible && block.sparse_factor_block_base) {
    factored = factor_sparse_latent_system(
      plan,
      block,
      size,
      rank,
      unit_diagonal,
      residual,
      U
    );
    sparse = factored;
  }
  if (!factored && block.block_base_eligible) {
    factored = factor_dense_latent_system(
      plan,
      block.sparse_factor_active_columns,
      size,
      rank,
      unit_diagonal,
      residual,
      U
    );
    sparse = false;
  }
  if (!factored) {
    return false;
  }

  const int n_rhs = precision_residual_only(output) ? 1 : size + 1;
  const size_t system_size = static_cast<size_t>(rank) *
    static_cast<size_t>(n_rhs);
  plan.conditional_rhs_work.resize(system_size);
  double *systems = plan.conditional_rhs_work.data();
  std::copy(
    plan.rhs_work.begin(),
    plan.rhs_work.begin() + rank,
    systems
  );
  if (n_rhs > 1) {
    std::fill(
      systems + rank,
      systems + system_size,
      0.0
    );
    for (int observation = 0; observation < size; ++observation) {
      for (int column : block.precision_active_columns[
          static_cast<size_t>(observation)
        ]) {
        systems[static_cast<size_t>(
          column + rank * (observation + 1)
        )] = base_precision[static_cast<size_t>(
          observation + size * (column + 1)
        )];
      }
    }
  }
  if (!solve_latent_system(
      plan,
      block,
      rank,
      n_rhs,
      sparse,
      systems
    )) {
    return false;
  }

  return finish_conditional_block(
    block,
    size,
    rank,
    block.precision_active_columns,
    base_precision_diagonal,
    base_precision,
    base_precision + size,
    systems,
    output
  );
}

void plan_conditional_values(
    CovariancePlan &plan,
    const double *mean_values,
    const std::vector<CovarianceFactor> &factors,
    const double *extra_values,
    const ConditionalOutput &output)
{
  for (int i = 0; i < plan.n; ++i) {
    if (!std::isfinite(extra_values[i]) || extra_values[i] < 0.0) {
      Rf_error("Known-V extra variances must be finite and non-negative.");
    }
  }

  for (BlockPlan &block : plan.blocks) {
    if (block.index.size() == 1) {
      double *covariance = plan.covariance_work.data();
      double *residual = plan.residual_work.data();
      fill_dense_block(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        covariance,
        residual
      );
      const double variance = covariance[0];
      const double value = residual[0];
      if (!(variance > 0.0) || !std::isfinite(variance) ||
          !std::isfinite(value)) {
        Rf_error("Known-V conditional covariance is not positive definite.");
      }
      const int global = block.index.front();
      if (output.log_density != nullptr) {
        const double two_pi = 6.283185307179586476925286766559;
        output.log_density[global] = -0.5 * (
          std::log(two_pi * variance) + value * value / variance
        );
      }
      if (output.residual != nullptr) {
        output.residual[global] = value;
      }
      if (output.variance != nullptr) {
        output.variance[global] = variance;
      }
      if (output.precision_residual != nullptr) {
        output.precision_residual[global] = value / variance;
      }
      continue;
    }

    // Specialized conditional routes are exact precision identities. Dense
    // Cholesky/inversion below is the target-preserving reference fallback.
    if (markov_block_conditional_values(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        output
      ) || fixed_known_group_block_conditional_values(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        output
      ) || spectral_block_conditional_values(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        output
      ) || diagonal_low_rank_conditional_values(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        output
      ) || block_base_low_rank_conditional_values(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        output
      )) {
      continue;
    }
    const int size = static_cast<int>(block.index.size());
    double *covariance = plan.covariance_work.data();
    double *precision_residual = plan.residual_work.data();
    if (!prepare_root_dense_block(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        covariance,
        precision_residual
      )) {
      fill_dense_block(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        covariance,
        precision_residual
      );
    }

    int info = 0;
    F77_CALL(dpotrf)("L", &size, covariance, &size, &info FCONE);
    if (info != 0) {
      Rf_error("Known-V conditional covariance is not positive definite.");
    }
    const int one = 1;
    F77_CALL(dpotrs)(
      "L",
      &size,
      &one,
      covariance,
      &size,
      precision_residual,
      &size,
      &info FCONE
    );
    if (info != 0) {
      Rf_error("Known-V conditional covariance solve failed.");
    }
    if (precision_residual_only(output)) {
      for (int local = 0; local < size; ++local) {
        const double weighted_residual = precision_residual[
          static_cast<size_t>(local)
        ];
        if (!std::isfinite(weighted_residual)) {
          Rf_error("Known-V precision residual is invalid.");
        }
        output.precision_residual[block.index[
          static_cast<size_t>(local)
        ]] = weighted_residual;
      }
      continue;
    }
    F77_CALL(dpotri)("L", &size, covariance, &size, &info FCONE);
    if (info != 0) {
      Rf_error("Known-V conditional precision could not be computed.");
    }
    for (int local = 0; local < size; ++local) {
      const double precision_diagonal = covariance[
        static_cast<size_t>(local + size * local)
      ];
      const double weighted_residual = precision_residual[
        static_cast<size_t>(local)
      ];
      if (!store_conditional_output(
          block.index[static_cast<size_t>(local)],
          precision_diagonal,
          weighted_residual,
          output
        )) {
        Rf_error("Known-V conditional precision is invalid.");
      }
    }
  }
}

SEXP plan_conditional_log_likelihood(
    CovariancePlan &plan, SEXP mean,
    const std::vector<CovarianceFactor> &factors,
    SEXP extra_variance)
{
  require_real_vector(mean, plan.n, "mean");
  require_real_vector(extra_variance, plan.n, "extra_variance");
  SEXP output = PROTECT(Rf_allocVector(REALSXP, plan.n));
  const ConditionalOutput destination = {
    REAL(output),
    nullptr,
    nullptr,
    nullptr
  };
  plan_conditional_values(
    plan,
    REAL(mean),
    factors,
    REAL(extra_variance),
    destination
  );
  UNPROTECT(1);
  return output;
}


bool affine_covariance_array_dimensions(
    SEXP covariance, int draws, int size, bool *draw_specific)
{
  if (TYPEOF(covariance) != REALSXP) {
    return false;
  }
  SEXP dimensions = Rf_getAttrib(covariance, R_DimSymbol);
  if (TYPEOF(dimensions) != INTSXP) {
    return false;
  }
  const int n_dimensions = Rf_length(dimensions);
  if (n_dimensions == 2 && INTEGER(dimensions)[0] == size &&
      INTEGER(dimensions)[1] == size) {
    *draw_specific = false;
    return true;
  }
  if (n_dimensions == 3 && INTEGER(dimensions)[0] == draws &&
      INTEGER(dimensions)[1] == size && INTEGER(dimensions)[2] == size) {
    *draw_specific = true;
    return true;
  }
  return false;
}


double affine_covariance_array_value(
    const double *values, bool draw_specific, int draws, int size,
    int draw, int row, int column)
{
  if (!draw_specific) {
    return values[static_cast<size_t>(row + size * column)];
  }
  return values[
    static_cast<size_t>(draw) + static_cast<size_t>(draws) *
      static_cast<size_t>(row + size * column)
  ];
}


std::vector<double> affine_rebase_candidates(
    const double *coefficients, int grids, double reference)
{
  std::vector<double> candidates;
  candidates.reserve(static_cast<size_t>(grids));
  for (int grid = 0; grid < grids; ++grid) {
    const double coefficient = coefficients[grid];
    if (coefficient == reference ||
        std::find(candidates.begin(), candidates.end(), coefficient) !=
          candidates.end()) {
      continue;
    }
    candidates.push_back(coefficient);
  }
  if (candidates.empty()) {
    return candidates;
  }

  std::vector<double> ordered = candidates;
  std::sort(ordered.begin(), ordered.end());
  const size_t middle = ordered.size() / 2;
  const double median = ordered.size() % 2 == 0 ?
    (ordered[middle - 1] + ordered[middle]) / 2.0 : ordered[middle];
  std::stable_sort(
    candidates.begin(),
    candidates.end(),
    [median](double left, double right) {
      return std::abs(left - median) < std::abs(right - median);
    }
  );
  return candidates;
}


bool affine_spectral_plan(
    int size, const double *base, const double *update,
    double reference, const std::vector<double> &rebase_candidates,
    std::vector<double> *root, std::vector<double> *eigenvectors,
    std::vector<double> *eigenvalues, std::vector<double> *work,
    double *plan_reference, double *log_determinant_base)
{
  std::vector<double> references;
  references.reserve(rebase_candidates.size() + 1L);
  references.push_back(reference);
  references.insert(
    references.end(), rebase_candidates.begin(), rebase_candidates.end()
  );
  const double one = 1.0;
  const int lwork = static_cast<int>(work->size());

  for (double candidate : references) {
    const double delta = candidate - reference;
    for (int index = 0; index < size * size; ++index) {
      (*root)[static_cast<size_t>(index)] =
        base[static_cast<size_t>(index)] +
        delta * update[static_cast<size_t>(index)];
      (*eigenvectors)[static_cast<size_t>(index)] =
        update[static_cast<size_t>(index)];
    }

    int info = 0;
    F77_CALL(dpotrf)("L", &size, root->data(), &size, &info FCONE);
    if (info != 0) {
      continue;
    }
    F77_CALL(dtrsm)(
      "L", "L", "N", "N", &size, &size, &one,
      root->data(), &size, eigenvectors->data(), &size
      FCONE FCONE FCONE FCONE
    );
    F77_CALL(dtrsm)(
      "R", "L", "T", "N", &size, &size, &one,
      root->data(), &size, eigenvectors->data(), &size
      FCONE FCONE FCONE FCONE
    );
    for (int column = 0; column < size; ++column) {
      for (int row = column + 1; row < size; ++row) {
        (*eigenvectors)[static_cast<size_t>(row + size * column)] =
          (*eigenvectors)[static_cast<size_t>(column + size * row)];
      }
    }
    F77_CALL(dsyev)(
      "V", "L", &size, eigenvectors->data(), &size,
      eigenvalues->data(), work->data(), &lwork, &info FCONE FCONE
    );
    if (info != 0) {
      continue;
    }

    double log_determinant = 0.0;
    bool finite = true;
    for (int index = 0; index < size; ++index) {
      const double diagonal = (*root)[static_cast<size_t>(index + size * index)];
      finite = finite && diagonal > 0.0 && std::isfinite(diagonal) &&
        std::isfinite((*eigenvalues)[static_cast<size_t>(index)]);
      log_determinant += 2.0 * std::log(diagonal);
    }
    if (!finite || !std::isfinite(log_determinant)) {
      continue;
    }
    *plan_reference = candidate;
    *log_determinant_base = log_determinant;
    return true;
  }
  return false;
}


SEXP plan_affine_grid_log_likelihood(
    CovariancePlan &plan, SEXP means, SEXP base_covariances,
    SEXP update_covariances, SEXP reference_coefficient, SEXP coefficients)
{
  if (!plan.factors.empty()) {
    Rf_error("Known-V affine grid plans cannot contain random factors.");
  }
  if (TYPEOF(means) != REALSXP || !Rf_isMatrix(means)) {
    Rf_error("Known-V affine grid means must be a numeric matrix.");
  }
  SEXP mean_dimensions = Rf_getAttrib(means, R_DimSymbol);
  const int draws = INTEGER(mean_dimensions)[1];
  if (INTEGER(mean_dimensions)[0] != plan.n || draws < 1) {
    Rf_error("Known-V affine grid means have inconsistent dimensions.");
  }
  if (TYPEOF(reference_coefficient) != REALSXP ||
      XLENGTH(reference_coefficient) != 1 ||
      !std::isfinite(REAL(reference_coefficient)[0]) ||
      TYPEOF(coefficients) != REALSXP || XLENGTH(coefficients) < 1 ||
      XLENGTH(coefficients) > static_cast<R_xlen_t>(2147483647)) {
    Rf_error("Known-V affine grid coefficients must be finite numeric values.");
  }
  const int grids = static_cast<int>(XLENGTH(coefficients));
  const double reference = REAL(reference_coefficient)[0];
  const double *coefficient_values = REAL(coefficients);
  for (int grid = 0; grid < grids; ++grid) {
    if (!std::isfinite(coefficient_values[grid])) {
      Rf_error("Known-V affine grid coefficients must be finite numeric values.");
    }
  }

  bool base_draw_specific = false;
  bool update_draw_specific = false;
  if (!affine_covariance_array_dimensions(
        base_covariances, draws, plan.n, &base_draw_specific
      ) || !affine_covariance_array_dimensions(
        update_covariances, draws, plan.n, &update_draw_specific
      ) || base_draw_specific != update_draw_specific) {
    Rf_error("Known-V affine covariance arrays have inconsistent dimensions.");
  }
  for (R_xlen_t index = 0; index < XLENGTH(base_covariances); ++index) {
    if (!std::isfinite(REAL(base_covariances)[index]) ||
        !std::isfinite(REAL(update_covariances)[index])) {
      Rf_error("Known-V affine covariance arrays must be finite.");
    }
  }
  for (R_xlen_t index = 0; index < XLENGTH(means); ++index) {
    if (!std::isfinite(REAL(means)[index])) {
      Rf_error("Known-V affine grid means must be finite.");
    }
  }

  size_t max_block_size = 0;
  for (const BlockPlan &block : plan.blocks) {
    max_block_size = std::max(max_block_size, block.index.size());
  }
  const int max_size = static_cast<int>(max_block_size);
  std::vector<double> base(max_block_size * max_block_size);
  std::vector<double> update(max_block_size * max_block_size);
  std::vector<double> root(max_block_size * max_block_size);
  std::vector<double> eigenvectors(max_block_size * max_block_size);
  std::vector<double> eigenvalues(max_block_size);
  std::vector<double> residual(max_block_size);
  std::vector<double> spectral(max_block_size);
  std::vector<double> work(static_cast<size_t>(std::max(1, 3 * max_size - 1)));
  const std::vector<double> rebase_candidates = affine_rebase_candidates(
    coefficient_values, grids, reference
  );
  const double *base_values = REAL(base_covariances);
  const double *update_values = REAL(update_covariances);
  const double *mean_values = REAL(means);
  const double log_two_pi = 1.837877066409345483560659472811;
  const int one_integer = 1;
  const double one = 1.0;
  const double zero = 0.0;

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, grids, draws));
  std::fill(
    REAL(out),
    REAL(out) + static_cast<size_t>(grids) * static_cast<size_t>(draws),
    0.0
  );
  for (const BlockPlan &block : plan.blocks) {
    const int size = static_cast<int>(block.index.size());
    for (int draw = 0; draw < draws; ++draw) {
      if ((draw & 63) == 0) {
        R_CheckUserInterrupt();
      }
      for (int column = 0; column < size; ++column) {
        const int global_column = block.index[static_cast<size_t>(column)];
        for (int row = 0; row < size; ++row) {
          const int global_row = block.index[static_cast<size_t>(row)];
          const size_t local = static_cast<size_t>(row + size * column);
          base[local] = affine_covariance_array_value(
            base_values, base_draw_specific, draws, plan.n,
            draw, global_row, global_column
          );
          update[local] = affine_covariance_array_value(
            update_values, update_draw_specific, draws, plan.n,
            draw, global_row, global_column
          );
        }
      }
      for (int column = 0; column < size; ++column) {
        for (int row = column + 1; row < size; ++row) {
          const size_t lower = static_cast<size_t>(row + size * column);
          const size_t upper = static_cast<size_t>(column + size * row);
          if (base[lower] != base[upper] || update[lower] != update[upper]) {
            UNPROTECT(1);
            Rf_error("Known-V affine covariance blocks must be exactly symmetric.");
          }
        }
      }

      double plan_reference = reference;
      double log_determinant_base = 0.0;
      if (!affine_spectral_plan(
            size, base.data(), update.data(), reference, rebase_candidates,
            &root, &eigenvectors, &eigenvalues, &work,
            &plan_reference, &log_determinant_base
          )) {
        UNPROTECT(1);
        return R_NilValue;
      }
      const double *draw_mean = mean_values +
        static_cast<size_t>(draw) * static_cast<size_t>(plan.n);
      for (int row = 0; row < size; ++row) {
        const int global = block.index[static_cast<size_t>(row)];
        residual[static_cast<size_t>(row)] =
          plan.y[static_cast<size_t>(global)] - draw_mean[global];
      }
      F77_CALL(dtrsv)(
        "L", "N", "N", &size, root.data(), &size,
        residual.data(), &one_integer FCONE FCONE FCONE
      );
      F77_CALL(dgemv)(
        "T", &size, &size, &one, eigenvectors.data(), &size,
        residual.data(), &one_integer, &zero, spectral.data(), &one_integer FCONE
      );

      for (int grid = 0; grid < grids; ++grid) {
        const double delta = coefficient_values[grid] - plan_reference;
        double log_determinant = log_determinant_base;
        double quadratic = 0.0;
        bool valid = true;
        for (int index = 0; index < size; ++index) {
          const double denominator = 1.0 + delta *
            eigenvalues[static_cast<size_t>(index)];
          if (!(denominator > 0.0) || !std::isfinite(denominator)) {
            valid = false;
            break;
          }
          const double value = spectral[static_cast<size_t>(index)];
          log_determinant += std::log(denominator);
          quadratic += value * value / denominator;
        }
        if (!valid || !std::isfinite(log_determinant) ||
            !std::isfinite(quadratic)) {
          UNPROTECT(1);
          return R_NilValue;
        }
        REAL(out)[static_cast<size_t>(grid) +
          static_cast<size_t>(grids) * static_cast<size_t>(draw)] +=
          -0.5 * (
            static_cast<double>(size) * log_two_pi +
            log_determinant + quadratic
          );
      }
    }
  }

  UNPROTECT(1);
  return out;
}


SEXP plan_factor_grid_log_likelihood(
    CovariancePlan &plan, SEXP means, SEXP random_covariance_states,
    SEXP extra_variances, SEXP update_grid)
{
  const int draws = require_batched_plan_inputs(
    plan,
    means,
    random_covariance_states,
    extra_variances
  );
  if (TYPEOF(update_grid) != VECSXP) {
    Rf_error("Known-V factor grid update must be a named list.");
  }
  SEXP family = list_element(update_grid, "family");
  SEXP factor_index_sexp = list_element(update_grid, "factor_index");
  if (TYPEOF(family) != STRSXP || XLENGTH(family) != 1 ||
      TYPEOF(factor_index_sexp) != INTSXP ||
      XLENGTH(factor_index_sexp) != 1) {
    Rf_error("Known-V factor grid metadata are invalid.");
  }
  const char *family_value = CHAR(STRING_ELT(family, 0));
  const bool factor_update = std::strcmp(family_value, "factor") == 0;
  const bool markov_update = std::strcmp(family_value, "markov") == 0;
  if (!factor_update && !markov_update) {
    Rf_error("Known-V factor grid family is unsupported.");
  }
  const int factor_index = INTEGER(factor_index_sexp)[0] - 1;
  if (factor_index < 0 ||
      factor_index >= static_cast<int>(plan.factors.size())) {
    Rf_error("Known-V factor grid identifies an invalid random factor.");
  }
  const PlanFactor &stored = plan.factors[static_cast<size_t>(factor_index)];
  if (stored.type == DENSE_FACTOR) {
    Rf_error("Known-V factor grids require a compiled coefficient factor.");
  }
  const int q = stored.n_columns;
  SEXP coefficient_scale = list_element(update_grid, "coefficient_scale");
  if (TYPEOF(coefficient_scale) != REALSXP ||
      !Rf_isMatrix(coefficient_scale)) {
    Rf_error("Known-V factor grid coefficient scales must be a numeric matrix.");
  }
  SEXP scale_dim = Rf_getAttrib(coefficient_scale, R_DimSymbol);
  if (INTEGER(scale_dim)[0] != draws || INTEGER(scale_dim)[1] != q) {
    Rf_error("Known-V factor grid coefficient scales have inconsistent dimensions.");
  }

  SEXP candidate_factor = R_NilValue;
  SEXP candidate_scale = R_NilValue;
  SEXP candidate_transition = R_NilValue;
  SEXP candidate_innovation = R_NilValue;
  int grids = 0;
  int component_index = -1;
  if (factor_update) {
    candidate_factor = list_element(update_grid, "coefficient_cholesky");
    candidate_scale = list_element(update_grid, "candidate_scale");
    SEXP component = list_element(update_grid, "component_index");
    if (TYPEOF(component) != INTSXP || XLENGTH(component) != 1 ||
        TYPEOF(candidate_scale) != REALSXP ||
        !Rf_isMatrix(candidate_scale)) {
      Rf_error("Known-V factor scale-grid metadata are invalid.");
    }
    component_index = INTEGER(component)[0] - 1;
    SEXP candidate_scale_dim = Rf_getAttrib(candidate_scale, R_DimSymbol);
    grids = INTEGER(candidate_scale_dim)[0];
    if (component_index < 0 || component_index >= q || grids < 1 ||
        INTEGER(candidate_scale_dim)[1] != draws) {
      Rf_error("Known-V factor scale-grid dimensions are inconsistent.");
    }
  } else {
    candidate_factor = list_element(update_grid, "candidate_cholesky");
    candidate_transition = list_element(update_grid, "candidate_transition");
    candidate_innovation = list_element(
      update_grid,
      "candidate_innovation_variance"
    );
    if (stored.coefficient_structure != MARKOV_COEFFICIENT ||
        TYPEOF(candidate_transition) != REALSXP ||
        !Rf_isMatrix(candidate_transition) ||
        TYPEOF(candidate_innovation) != REALSXP ||
        !Rf_isMatrix(candidate_innovation)) {
      Rf_error("Known-V Markov grid metadata are invalid.");
    }
    SEXP transition_dim = Rf_getAttrib(candidate_transition, R_DimSymbol);
    SEXP innovation_dim = Rf_getAttrib(candidate_innovation, R_DimSymbol);
    grids = INTEGER(transition_dim)[0];
    if (grids < 1 || INTEGER(transition_dim)[1] != q - 1 ||
        INTEGER(innovation_dim)[0] != grids ||
        INTEGER(innovation_dim)[1] != q - 1) {
      Rf_error("Known-V Markov grid dimensions are inconsistent.");
    }
  }
  SEXP factor_dim = Rf_getAttrib(candidate_factor, R_DimSymbol);
  if (TYPEOF(candidate_factor) != REALSXP ||
      TYPEOF(factor_dim) != INTSXP || Rf_length(factor_dim) != 3 ||
      INTEGER(factor_dim)[0] != (factor_update ? draws : grids) ||
      INTEGER(factor_dim)[1] != q || INTEGER(factor_dim)[2] != q) {
    Rf_error("Known-V factor grid Cholesky states have inconsistent dimensions.");
  }

  const double *scale_values = REAL(coefficient_scale);
  const double *candidate_factor_values = REAL(candidate_factor);
  const double *candidate_scale_values = factor_update ?
    REAL(candidate_scale) : nullptr;
  const double *transition_values = markov_update ?
    REAL(candidate_transition) : nullptr;
  const double *innovation_values = markov_update ?
    REAL(candidate_innovation) : nullptr;
  for (R_xlen_t index = 0; index < XLENGTH(coefficient_scale); ++index) {
    if (!std::isfinite(scale_values[index]) || scale_values[index] < 0.0) {
      Rf_error("Known-V factor grid coefficient scales must be finite and non-negative.");
    }
  }
  for (R_xlen_t index = 0; index < XLENGTH(candidate_factor); ++index) {
    if (!std::isfinite(candidate_factor_values[index])) {
      Rf_error("Known-V factor grid Cholesky states must be finite.");
    }
  }
  if (factor_update) {
    for (R_xlen_t index = 0; index < XLENGTH(candidate_scale); ++index) {
      if (!std::isfinite(candidate_scale_values[index]) ||
          candidate_scale_values[index] < 0.0) {
        Rf_error("Known-V candidate coefficient scales must be finite and non-negative.");
      }
    }
  } else {
    for (R_xlen_t index = 0; index < XLENGTH(candidate_transition); ++index) {
      if (!std::isfinite(transition_values[index]) ||
          !std::isfinite(innovation_values[index]) ||
          !(innovation_values[index] > 0.0)) {
        Rf_error("Known-V candidate Markov states are invalid.");
      }
    }
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, grids, draws));
  const double *mean_values = REAL(means);
  const double *extra_values = REAL(extra_variances);
  std::vector<double> coefficient_factor(static_cast<size_t>(q * q));
  std::vector<double> current_scale(static_cast<size_t>(q));
  std::vector<double> transition(static_cast<size_t>(std::max(0, q - 1)));
  std::vector<double> innovation(static_cast<size_t>(std::max(0, q - 1)));

  for (int draw = 0; draw < draws; ++draw) {
    if ((draw & 31) == 0) {
      R_CheckUserInterrupt();
    }
    std::vector<CovarianceFactor> factors = covariance_states(
      VECTOR_ELT(random_covariance_states, draw),
      plan
    );
    CovarianceFactor &focal = factors[static_cast<size_t>(factor_index)];
    for (int column = 0; column < q; ++column) {
      current_scale[static_cast<size_t>(column)] = scale_values[
        static_cast<size_t>(draw) +
        static_cast<size_t>(draws) * static_cast<size_t>(column)
      ];
    }
    for (int grid = 0; grid < grids; ++grid) {
      if (factor_update) {
        current_scale[static_cast<size_t>(component_index)] =
          candidate_scale_values[
            static_cast<size_t>(grid) +
            static_cast<size_t>(grids) * static_cast<size_t>(draw)
          ];
      }
      const int factor_states = factor_update ? draws : grids;
      const int factor_state = factor_update ? draw : grid;
      for (int column = 0; column < q; ++column) {
        for (int row = 0; row < q; ++row) {
          coefficient_factor[static_cast<size_t>(row + q * column)] =
            candidate_factor_values[
              static_cast<size_t>(factor_state) +
              static_cast<size_t>(factor_states) *
                static_cast<size_t>(row + q * column)
            ] * current_scale[static_cast<size_t>(row)];
        }
      }
      focal.coefficient_factor = coefficient_factor.data();
      if (focal.coefficient_structure == MARKOV_COEFFICIENT) {
        focal.coefficient_scale = current_scale.data();
      }
      if (markov_update) {
        for (int column = 0; column < q - 1; ++column) {
          transition[static_cast<size_t>(column)] = transition_values[
            static_cast<size_t>(grid) +
            static_cast<size_t>(grids) * static_cast<size_t>(column)
          ];
          innovation[static_cast<size_t>(column)] = innovation_values[
            static_cast<size_t>(grid) +
            static_cast<size_t>(grids) * static_cast<size_t>(column)
          ];
        }
        focal.markov_transition = transition.data();
        focal.markov_innovation_variance = innovation.data();
      }
      REAL(out)[static_cast<size_t>(grid) +
        static_cast<size_t>(grids) * static_cast<size_t>(draw)] =
        plan_log_likelihood_values(
          plan,
          mean_values + static_cast<size_t>(draw) *
            static_cast<size_t>(plan.n),
          factors,
          extra_values + static_cast<size_t>(draw) *
            static_cast<size_t>(plan.n)
        );
    }
  }

  UNPROTECT(1);
  return out;
}
