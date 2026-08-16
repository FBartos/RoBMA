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
