double plan_log_likelihood(CovariancePlan &plan, SEXP mean,
                           const std::vector<CovarianceFactor> &factors,
                           SEXP extra_variance)
{
  require_real_vector(mean, plan.n, "mean");
  require_real_vector(extra_variance, plan.n, "extra_variance");
  const double *mean_values = REAL(mean);
  const double *extra_values = REAL(extra_variance);
  for (int i = 0; i < plan.n; ++i) {
    if (!std::isfinite(extra_values[i]) || extra_values[i] < 0.0) {
      Rf_error("Known-V bridge extra variances must be finite and non-negative.");
    }
  }

  double log_likelihood = 0.0;
  for (BlockPlan &block : plan.blocks) {
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
    double *output)
{
  // Q = B^-1 - C K^-1 C' and Qr = B^-1 r - C K^-1 h.
  const double *latent_residual_solution = solutions;
  const double log_two_pi = 1.837877066409345483560659472811;
  for (int observation = 0; observation < size; ++observation) {
    const double *precision_solution = solutions +
      static_cast<size_t>((observation + 1) * rank);
    double residual_correction = 0.0;
    double diagonal_correction = 0.0;
    for (int column : active_columns[static_cast<size_t>(observation)]) {
      const double precision_value = precision_factor[static_cast<size_t>(
        observation + size * column
      )];
      residual_correction += precision_value *
        latent_residual_solution[static_cast<size_t>(column)];
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
    const double precision_residual =
      base_precision_residual[static_cast<size_t>(observation)] -
      residual_correction;
    if (!(precision_diagonal > cancellation_bound) ||
        !std::isfinite(precision_residual)) {
      return false;
    }
    output[block.index[static_cast<size_t>(observation)]] = 0.5 * (
      std::log(precision_diagonal) - log_two_pi -
      precision_residual * precision_residual / precision_diagonal
    );
  }
  return true;
}

bool diagonal_low_rank_conditional_log_likelihood(
    CovariancePlan &plan,
    BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *output)
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
    const double log_two_pi = 1.837877066409345483560659472811;
    for (int observation = 0; observation < size; ++observation) {
      const double variance = diagonal[static_cast<size_t>(observation)];
      const double value = residual[static_cast<size_t>(observation)];
      output[block.index[static_cast<size_t>(observation)]] = -0.5 * (
        log_two_pi + std::log(variance) + value * value / variance
      );
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

  const size_t system_size = static_cast<size_t>(rank) *
    static_cast<size_t>(size + 1);
  plan.conditional_rhs_work.resize(system_size);
  double *systems = plan.conditional_rhs_work.data();
  std::copy(
    plan.rhs_work.begin(),
    plan.rhs_work.begin() + rank,
    systems
  );
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
  if (!solve_latent_system(
      plan,
      block,
      rank,
      size + 1,
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

bool block_base_low_rank_conditional_log_likelihood(
    CovariancePlan &plan,
    BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *output)
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

  const size_t system_size = static_cast<size_t>(rank) *
    static_cast<size_t>(size + 1);
  plan.conditional_rhs_work.resize(system_size);
  double *systems = plan.conditional_rhs_work.data();
  std::copy(
    plan.rhs_work.begin(),
    plan.rhs_work.begin() + rank,
    systems
  );
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
  if (!solve_latent_system(
      plan,
      block,
      rank,
      size + 1,
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

SEXP plan_conditional_log_likelihood(
    CovariancePlan &plan, SEXP mean,
    const std::vector<CovarianceFactor> &factors,
    SEXP extra_variance)
{
  require_real_vector(mean, plan.n, "mean");
  require_real_vector(extra_variance, plan.n, "extra_variance");
  const double *mean_values = REAL(mean);
  const double *extra_values = REAL(extra_variance);
  for (int i = 0; i < plan.n; ++i) {
    if (!std::isfinite(extra_values[i]) || extra_values[i] < 0.0) {
      Rf_error("Known-V extra variances must be finite and non-negative.");
    }
  }

  SEXP output = PROTECT(Rf_allocVector(REALSXP, plan.n));
  double *output_values = REAL(output);
  const double log_two_pi = 1.837877066409345483560659472811;
  for (BlockPlan &block : plan.blocks) {
    if (markov_block_conditional_log_likelihood(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        output_values
      ) || fixed_known_group_block_conditional_log_likelihood(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        output_values
      ) || diagonal_low_rank_conditional_log_likelihood(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        output_values
      ) || block_base_low_rank_conditional_log_likelihood(
        plan,
        block,
        factors,
        mean_values,
        extra_values,
        output_values
      )) {
      continue;
    }
    const int size = static_cast<int>(block.index.size());
    double *covariance = plan.covariance_work.data();
    double *precision_residual = plan.residual_work.data();
    fill_dense_block(
      plan,
      block,
      factors,
      mean_values,
      extra_values,
      covariance,
      precision_residual
    );

    int info = 0;
    F77_CALL(dpotrf)("L", &size, covariance, &size, &info FCONE);
    if (info != 0) {
      UNPROTECT(1);
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
      UNPROTECT(1);
      Rf_error("Known-V conditional covariance solve failed.");
    }
    F77_CALL(dpotri)("L", &size, covariance, &size, &info FCONE);
    if (info != 0) {
      UNPROTECT(1);
      Rf_error("Known-V conditional precision could not be computed.");
    }
    for (int local = 0; local < size; ++local) {
      const double precision_diagonal = covariance[
        static_cast<size_t>(local + size * local)
      ];
      const double weighted_residual = precision_residual[
        static_cast<size_t>(local)
      ];
      if (!(precision_diagonal > 0.0) ||
          !std::isfinite(precision_diagonal) ||
          !std::isfinite(weighted_residual)) {
        UNPROTECT(1);
        Rf_error("Known-V conditional precision is invalid.");
      }
      output_values[block.index[static_cast<size_t>(local)]] = 0.5 * (
        std::log(precision_diagonal) - log_two_pi -
        weighted_residual * weighted_residual / precision_diagonal
      );
    }
  }
  UNPROTECT(1);
  return output;
}
