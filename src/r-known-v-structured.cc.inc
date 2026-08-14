bool prepare_block_base_sampling_block(
    CovariancePlan &plan,
    const BlockPlan &block,
    const SamplingBlock &sampling_block,
    const double *extra,
    int size,
    int rank,
    double *residual,
    double *U,
    double *log_determinant,
    double *base_quadratic,
    double *base_precision_diagonal,
    double *base_precision)
{
  const int block_size = static_cast<int>(sampling_block.local_rows.size());
  double *covariance = plan.covariance_work.data();
  double *base_rhs = plan.base_rhs_work.data();
  bool fixed_base = sampling_block.factor_eligible;
  for (int row = 0; row < block_size; ++row) {
    const int local = sampling_block.local_rows[static_cast<size_t>(row)];
    const int global = block.index[static_cast<size_t>(local)];
    fixed_base = fixed_base && extra[global] == 0.0;
  }
  if (fixed_base) {
    std::copy(
      sampling_block.factor.begin(),
      sampling_block.factor.end(),
      covariance
    );
    if (log_determinant != nullptr) {
      *log_determinant += sampling_block.log_determinant;
    }
  } else {
    std::copy(
      sampling_block.covariance.begin(),
      sampling_block.covariance.end(),
      covariance
    );
  }
  for (int row = 0; row < block_size; ++row) {
    const int local = sampling_block.local_rows[static_cast<size_t>(row)];
    const int global = block.index[static_cast<size_t>(local)];
    if (!fixed_base) {
      covariance[static_cast<size_t>(row + block_size * row)] +=
        extra[global];
    }
    base_rhs[static_cast<size_t>(row)] = residual[static_cast<size_t>(local)];
    for (int column = 0; column < rank; ++column) {
      base_rhs[static_cast<size_t>(row + block_size * (column + 1))] =
        U[static_cast<size_t>(local + size * column)];
    }
  }

  int info = 0;
  if (!fixed_base) {
    F77_CALL(dpotrf)(
      "L", &block_size, covariance, &block_size, &info FCONE
    );
    if (info != 0) {
      return false;
    }
    if (log_determinant != nullptr) {
      for (int row = 0; row < block_size; ++row) {
        *log_determinant += 2.0 * std::log(covariance[static_cast<size_t>(
          row + block_size * row
        )]);
      }
    }
  }

  const double one = 1.0;
  const int rhs_columns = rank + 1;
  F77_CALL(dtrsm)(
    "L", "L", "N", "N", &block_size, &rhs_columns, &one,
    covariance, &block_size, base_rhs, &block_size
    FCONE FCONE FCONE FCONE
  );
  for (int row = 0; row < block_size; ++row) {
    const int local = sampling_block.local_rows[static_cast<size_t>(row)];
    const double transformed_residual = base_rhs[static_cast<size_t>(row)];
    residual[static_cast<size_t>(local)] = transformed_residual;
    if (base_quadratic != nullptr) {
      *base_quadratic += transformed_residual * transformed_residual;
    }
    for (int column = 0; column < rank; ++column) {
      U[static_cast<size_t>(local + size * column)] = base_rhs[
        static_cast<size_t>(row + block_size * (column + 1))
      ];
    }
  }

  if (base_precision == nullptr) {
    return true;
  }
  F77_CALL(dtrsm)(
    "L", "L", "T", "N", &block_size, &rhs_columns, &one,
    covariance, &block_size, base_rhs, &block_size
    FCONE FCONE FCONE FCONE
  );
  for (int row = 0; row < block_size; ++row) {
    const int local = sampling_block.local_rows[static_cast<size_t>(row)];
    for (int column = 0; column <= rank; ++column) {
      base_precision[static_cast<size_t>(local + size * column)] =
        base_rhs[static_cast<size_t>(row + block_size * column)];
    }
  }
  if (fixed_base) {
    for (int row = 0; row < block_size; ++row) {
      const int local = sampling_block.local_rows[static_cast<size_t>(row)];
      base_precision_diagonal[static_cast<size_t>(local)] =
        sampling_block.precision_diagonal[static_cast<size_t>(row)];
    }
  } else {
    F77_CALL(dpotri)(
      "L", &block_size, covariance, &block_size, &info FCONE
    );
    if (info != 0) {
      return false;
    }
    for (int row = 0; row < block_size; ++row) {
      const int local = sampling_block.local_rows[static_cast<size_t>(row)];
      base_precision_diagonal[static_cast<size_t>(local)] = covariance[
        static_cast<size_t>(row + block_size * row)
      ];
    }
  }
  return true;
}

bool block_base_low_rank_log_likelihood(
    CovariancePlan &plan,
    BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *value)
{
  if (!block.block_base_eligible &&
      !(block.sparse_factor_eligible && block.sparse_factor_block_base)) {
    return false;
  }
  const int size = static_cast<int>(block.index.size());
  const int rank = block.rank;
  double *diagonal = plan.diagonal_work.data();
  double *residual = plan.residual_work.data();
  double *U = fill_low_rank_factor(plan, block, factors);
  std::fill(diagonal, diagonal + size, 1.0);
  for (int local = 0; local < size; ++local) {
    const int global = block.index[static_cast<size_t>(local)];
    if (!std::isfinite(mean[global])) {
      return false;
    }
    residual[static_cast<size_t>(local)] =
      plan.y[static_cast<size_t>(global)] - mean[global];
  }

  double log_determinant = 0.0;
  double base_quadratic = 0.0;
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
        &log_determinant,
        &base_quadratic,
        nullptr,
        nullptr
      )) {
      return false;
    }
  }

  if (block.sparse_factor_block_base && sparse_factor_log_likelihood(
      plan,
      block,
      size,
      rank,
      diagonal,
      residual,
      U,
      log_determinant,
      base_quadratic,
      value
    )) {
    return true;
  }
  if (!block.block_base_eligible) {
    return false;
  }
  return diagonal_low_rank_log_likelihood(
    plan,
    size,
    rank,
    diagonal,
    residual,
    U,
    log_determinant,
    base_quadratic,
    value
  );
}

bool spectral_block_log_likelihood(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *value)
{
  if (!block.spectral_eligible) {
    return false;
  }
  const int size = static_cast<int>(block.index.size());
  const int rank = block.rank;
  const double shift = extra[block.index.front()];
  if (!std::isfinite(shift) || shift < 0.0) {
    return false;
  }
  for (int local = 0; local < size; ++local) {
    const int global = block.index[static_cast<size_t>(local)];
    if (extra[global] != shift || !std::isfinite(mean[global])) {
      return false;
    }
  }

  double *diagonal = plan.diagonal_work.data();
  double *residual = plan.residual_work.data();
  double log_determinant = 0.0;
  double base_quadratic = 0.0;
  for (int eigen = 0; eigen < size; ++eigen) {
    const double variance = block.sampling_eigenvalues[
      static_cast<size_t>(eigen)
    ] + shift;
    if (!(variance > 0.0) || !std::isfinite(variance)) {
      return false;
    }
    double transformed_residual = 0.0;
    for (int local = 0; local < size; ++local) {
      const int global = block.index[static_cast<size_t>(local)];
      transformed_residual += block.sampling_eigenvectors[
        static_cast<size_t>(local + size * eigen)
      ] * (plan.y[static_cast<size_t>(global)] - mean[global]);
    }
    diagonal[static_cast<size_t>(eigen)] = variance;
    residual[static_cast<size_t>(eigen)] = transformed_residual;
    log_determinant += std::log(variance);
    base_quadratic += transformed_residual * transformed_residual / variance;
  }

  double *U = fill_low_rank_factor(plan, block, factors);
  double *transformed_U = plan.covariance_work.data();
  for (int column = 0; column < rank; ++column) {
    for (int eigen = 0; eigen < size; ++eigen) {
      double transformed = 0.0;
      for (int local = 0; local < size; ++local) {
        transformed += block.sampling_eigenvectors[
          static_cast<size_t>(local + size * eigen)
        ] * U[static_cast<size_t>(local + size * column)];
      }
      transformed_U[static_cast<size_t>(eigen + size * column)] = transformed;
    }
  }

  return diagonal_low_rank_log_likelihood(
    plan,
    size,
    rank,
    diagonal,
    residual,
    transformed_U,
    log_determinant,
    base_quadratic,
    value
  );
}

bool markov_block_log_likelihood(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *value)
{
  if (!block.markov_eligible) {
    return false;
  }
  const CovarianceFactor &factor = factors[static_cast<size_t>(
    block.markov_factor
  )];
  const int q = factor.n_columns;
  double state_mean = 0.0;
  double state_variance = 1.0;
  double log_likelihood = 0.0;
  const double log_two_pi = 1.837877066409345483560659472811;
  for (int level = 0; level < q; ++level) {
    if (level > 0) {
      const double transition = factor.markov_transition[level - 1];
      state_mean *= transition;
      state_variance = transition * transition * state_variance +
        factor.markov_innovation_variance[level - 1];
      if (!(state_variance > 0.0) || !std::isfinite(state_variance)) {
        return false;
      }
    }
    for (int local : block.markov_rows_by_level[static_cast<size_t>(level)]) {
      const int global = block.index[static_cast<size_t>(local)];
      const double observation_variance = plan.sampling_covariance[
        static_cast<size_t>(global + plan.n * global)
      ] + extra[global];
      if (!(observation_variance > 0.0) ||
          !std::isfinite(observation_variance) ||
          !std::isfinite(mean[global])) {
        return false;
      }
      double loading = block.markov_loading[static_cast<size_t>(local)] *
        factor.coefficient_scale[level];
      if (factor.type == ROW_GROUP_FACTOR) {
        loading *= factor.row_scale[global];
      }
      const double residual = plan.y[static_cast<size_t>(global)] -
        mean[global] - loading * state_mean;
      const double innovation_variance = observation_variance +
        loading * loading * state_variance;
      if (!(innovation_variance > 0.0) ||
          !std::isfinite(innovation_variance) ||
          !std::isfinite(residual)) {
        return false;
      }
      log_likelihood += -0.5 * (
        log_two_pi + std::log(innovation_variance) +
        residual * residual / innovation_variance
      );
      state_mean += state_variance * loading * residual /
        innovation_variance;
      state_variance *= observation_variance / innovation_variance;
      if (!(state_variance > 0.0) || !std::isfinite(state_mean)) {
        return false;
      }
    }
  }
  *value = log_likelihood;
  return true;
}

bool markov_block_conditional_log_likelihood(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *output)
{
  if (!block.markov_eligible) {
    return false;
  }
  const CovarianceFactor &factor = factors[static_cast<size_t>(
    block.markov_factor
  )];
  const int q = factor.n_columns;
  plan.markov_diagonal_work.assign(static_cast<size_t>(q), 0.0);
  plan.markov_off_diagonal_work.assign(
    static_cast<size_t>(q - 1),
    0.0
  );
  plan.markov_rhs_work.assign(static_cast<size_t>(q), 0.0);
  double *diagonal = plan.markov_diagonal_work.data();
  double *off_diagonal = plan.markov_off_diagonal_work.data();
  double *rhs = plan.markov_rhs_work.data();
  double *observation_variance = plan.diagonal_work.data();
  double *observation_loading = plan.low_rank_work.data();
  double *observation_residual = plan.residual_work.data();
  diagonal[0] = 1.0;
  for (int level = 1; level < q; ++level) {
    const double transition = factor.markov_transition[level - 1];
    const double innovation = factor.markov_innovation_variance[level - 1];
    const double inverse_innovation = 1.0 / innovation;
    diagonal[level - 1] += transition * transition * inverse_innovation;
    diagonal[level] += inverse_innovation;
    off_diagonal[level - 1] = -transition * inverse_innovation;
  }
  for (int level = 0; level < q; ++level) {
    for (int local : block.markov_rows_by_level[static_cast<size_t>(level)]) {
      const int global = block.index[static_cast<size_t>(local)];
      const double variance = plan.sampling_covariance[
        static_cast<size_t>(global + plan.n * global)
      ] + extra[global];
      if (!(variance > 0.0) || !std::isfinite(variance) ||
          !std::isfinite(mean[global])) {
        return false;
      }
      double loading = block.markov_loading[static_cast<size_t>(local)] *
        factor.coefficient_scale[level];
      if (factor.type == ROW_GROUP_FACTOR) {
        loading *= factor.row_scale[global];
      }
      const double residual = plan.y[static_cast<size_t>(global)] - mean[global];
      const double inverse_variance = 1.0 / variance;
      diagonal[level] += loading * loading * inverse_variance;
      rhs[level] += loading * residual * inverse_variance;
      observation_variance[local] = variance;
      observation_loading[local] = loading;
      observation_residual[local] = residual;
    }
  }

  plan.markov_factor_work.resize(static_cast<size_t>(q - 1));
  double *subdiagonal_factor = plan.markov_factor_work.data();
  if (!(diagonal[0] > 0.0) || !std::isfinite(diagonal[0])) {
    return false;
  }
  for (int level = 1; level < q; ++level) {
    subdiagonal_factor[level - 1] =
      off_diagonal[level - 1] / diagonal[level - 1];
    diagonal[level] -= subdiagonal_factor[level - 1] *
      off_diagonal[level - 1];
    if (!(diagonal[level] > 0.0) || !std::isfinite(diagonal[level])) {
      return false;
    }
  }

  plan.markov_solution_work.assign(rhs, rhs + q);
  double *solution = plan.markov_solution_work.data();
  for (int level = 1; level < q; ++level) {
    solution[level] -= subdiagonal_factor[level - 1] * solution[level - 1];
  }
  solution[q - 1] /= diagonal[q - 1];
  for (int level = q - 2; level >= 0; --level) {
    solution[level] = solution[level] / diagonal[level] -
      subdiagonal_factor[level] * solution[level + 1];
  }

  plan.markov_variance_work.resize(static_cast<size_t>(q));
  double *posterior_variance = plan.markov_variance_work.data();
  posterior_variance[q - 1] = 1.0 / diagonal[q - 1];
  for (int level = q - 2; level >= 0; --level) {
    posterior_variance[level] = 1.0 / diagonal[level] +
      subdiagonal_factor[level] * subdiagonal_factor[level] *
      posterior_variance[level + 1];
  }

  const double log_two_pi = 1.837877066409345483560659472811;
  for (int level = 0; level < q; ++level) {
    for (int local : block.markov_rows_by_level[static_cast<size_t>(level)]) {
      const double variance = observation_variance[local];
      const double loading = observation_loading[local];
      const double residual = observation_residual[local];
      const double precision_loading = loading / variance;
      const double precision_diagonal = 1.0 / variance -
        precision_loading * precision_loading * posterior_variance[level];
      const double precision_residual = residual / variance -
        precision_loading * solution[level];
      if (!(precision_diagonal > 0.0) ||
          !std::isfinite(precision_diagonal) ||
          !std::isfinite(precision_residual)) {
        return false;
      }
      output[block.index[static_cast<size_t>(local)]] = 0.5 * (
        std::log(precision_diagonal) - log_two_pi -
        precision_residual * precision_residual / precision_diagonal
      );
    }
  }
  return true;
}

bool prepare_fixed_known_group_block(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra)
{
  if (!block.fixed_known_group_eligible) {
    return false;
  }
  const int size = static_cast<int>(block.index.size());
  const int rank = block.rank;
  for (int local = 0; local < size; ++local) {
    const int global = block.index[static_cast<size_t>(local)];
    if (extra[global] != 0.0) {
      return false;
    }
  }

  const CovarianceFactor &factor = factors.front();
  double *scales = plan.scale_work.data();
  for (const LowRankGroup &group : block.low_rank_groups) {
    for (int coefficient = 0;
         coefficient < factor.n_columns;
         ++coefficient) {
      scales[static_cast<size_t>(group.column_offset + coefficient)] =
        factor.coefficient_factor[static_cast<size_t>(
          coefficient + factor.n_columns * coefficient
        )];
    }
  }

  double *residual = plan.residual_work.data();
  double *precision_residual = plan.diagonal_work.data();
  for (int local = 0; local < size; ++local) {
    const int global = block.index[static_cast<size_t>(local)];
    residual[static_cast<size_t>(local)] = plan.y[static_cast<size_t>(
      global
    )] - mean[global];
    if (!std::isfinite(residual[static_cast<size_t>(local)])) {
      Rf_error("Known-V bridge residuals must be finite.");
    }
    precision_residual[static_cast<size_t>(local)] =
      residual[static_cast<size_t>(local)] *
      block.fixed_known_group_precision[static_cast<size_t>(local)];
  }

  const double one = 1.0;
  const double zero = 0.0;
  const int increment = 1;
  double *rhs = plan.rhs_work.data();
  F77_CALL(dgemv)(
    "T", &size, &rank, &one, block.fixed_known_group_basis.data(),
    &size, precision_residual, &increment, &zero, rhs, &increment FCONE
  );
  for (int column = 0; column < rank; ++column) {
    rhs[static_cast<size_t>(column)] *= scales[static_cast<size_t>(column)];
  }

  double *small = plan.small_matrix_work.data();
  for (int column = 0; column < rank; ++column) {
    for (int row = column; row < rank; ++row) {
      small[static_cast<size_t>(row + rank * column)] =
        block.fixed_known_group_crossproduct[static_cast<size_t>(
          row + rank * column
        )] * scales[static_cast<size_t>(row)] *
        scales[static_cast<size_t>(column)];
    }
    small[static_cast<size_t>(column + rank * column)] += 1.0;
  }
  int info = 0;
  F77_CALL(dpotrf)("L", &rank, small, &rank, &info FCONE);
  if (info != 0) {
    return false;
  }
  F77_CALL(dpotrs)(
    "L", &rank, &increment, small, &rank, rhs, &rank, &info FCONE
  );
  if (info != 0) {
    return false;
  }

  for (int column = 0; column < rank; ++column) {
    scales[static_cast<size_t>(column)] *= rhs[static_cast<size_t>(column)];
  }
  const double minus_one = -1.0;
  F77_CALL(dgemv)(
    "N", &size, &rank, &minus_one, block.fixed_known_group_basis.data(),
    &size, scales, &increment, &one, residual, &increment FCONE
  );
  return true;
}

bool fixed_known_group_block_log_likelihood(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *value)
{
  if (!prepare_fixed_known_group_block(
      plan,
      block,
      factors,
      mean,
      extra
    )) {
    return false;
  }

  const int size = static_cast<int>(block.index.size());
  const int rank = block.rank;
  const int increment = 1;
  const double *small = plan.small_matrix_work.data();
  const double *rhs = plan.rhs_work.data();
  const double *residual = plan.residual_work.data();
  double log_determinant = block.fixed_known_group_log_determinant;
  for (int column = 0; column < rank; ++column) {
    log_determinant += 2.0 * std::log(small[static_cast<size_t>(
      column + rank * column
    )]);
  }
  double quadratic = F77_CALL(ddot)(
    &rank, rhs, &increment, rhs, &increment
  );
  for (int local = 0; local < size; ++local) {
    const double current = residual[static_cast<size_t>(local)];
    quadratic += current * current *
      block.fixed_known_group_precision[static_cast<size_t>(local)];
  }
  if (!std::isfinite(log_determinant) || !std::isfinite(quadratic)) {
    Rf_error("Known-V fixed-basis likelihood is not finite.");
  }
  const double log_two_pi = 1.837877066409345483560659472811;
  *value = -0.5 * (
    static_cast<double>(size) * log_two_pi + log_determinant + quadratic
  );
  return true;
}

bool fixed_known_group_block_conditional_log_likelihood(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *output)
{
  if (!prepare_fixed_known_group_block(
      plan,
      block,
      factors,
      mean,
      extra
    )) {
    return false;
  }

  const int size = static_cast<int>(block.index.size());
  const int rank = block.rank;
  const CovarianceFactor &factor = factors.front();
  double *transformed_basis = plan.low_rank_work.data();
  std::copy(
    block.fixed_known_group_basis.begin(),
    block.fixed_known_group_basis.end(),
    transformed_basis
  );
  for (const LowRankGroup &group : block.low_rank_groups) {
    for (int coefficient = 0;
         coefficient < factor.n_columns;
         ++coefficient) {
      const int column = group.column_offset + coefficient;
      const double scale = factor.coefficient_factor[static_cast<size_t>(
        coefficient + factor.n_columns * coefficient
      )];
      for (int local = 0; local < size; ++local) {
        transformed_basis[static_cast<size_t>(
          local + size * column
        )] *= scale;
      }
    }
  }

  const double one = 1.0;
  double *small = plan.small_matrix_work.data();
  F77_CALL(dtrsm)(
    "R", "L", "T", "N", &size, &rank, &one, small, &rank,
    transformed_basis, &size FCONE FCONE FCONE FCONE
  );
  const double log_two_pi = 1.837877066409345483560659472811;
  for (int local = 0; local < size; ++local) {
    const double precision = block.fixed_known_group_precision[
      static_cast<size_t>(local)
    ];
    double adjustment = 0.0;
    for (int column = 0; column < rank; ++column) {
      const double current = transformed_basis[static_cast<size_t>(
        local + size * column
      )];
      adjustment += current * current;
    }
    const double precision_diagonal = precision -
      precision * precision * adjustment;
    const double precision_residual = precision *
      plan.residual_work[static_cast<size_t>(local)];
    if (!(precision_diagonal > 0.0) ||
        !std::isfinite(precision_diagonal) ||
        !std::isfinite(precision_residual)) {
      return false;
    }
    output[block.index[static_cast<size_t>(local)]] = 0.5 * (
      std::log(precision_diagonal) - log_two_pi -
      precision_residual * precision_residual / precision_diagonal
    );
  }
  return true;
}
