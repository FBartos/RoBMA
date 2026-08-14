PlanFactor covariance_factor_plan(SEXP factor, int n)
{
  PlanFactor out = {};
  SEXP type = list_element(factor, "type");
  if (TYPEOF(type) != STRSXP || XLENGTH(type) != 1) {
    Rf_error("Random covariance factor type must be one string.");
  }
  const char *type_value = CHAR(STRING_ELT(type, 0));
  if (std::strcmp(type_value, "dense") == 0) {
    out.type = DENSE_FACTOR;
    return out;
  }

  SEXP model_matrix = list_element(factor, "model_matrix");
  if (TYPEOF(model_matrix) != REALSXP || !Rf_isMatrix(model_matrix)) {
    Rf_error("Random covariance factor model matrix must be double.");
  }
  SEXP model_dimensions = Rf_getAttrib(model_matrix, R_DimSymbol);
  out.n_columns = INTEGER(model_dimensions)[1];
  if (INTEGER(model_dimensions)[0] != n || out.n_columns < 1) {
    Rf_error("Random covariance factor model matrix has inconsistent dimensions.");
  }
  const double *model_values = REAL(model_matrix);
  for (int i = 0; i < n * out.n_columns; ++i) {
    if (!std::isfinite(model_values[i])) {
      Rf_error("Random covariance factor model matrix must be finite.");
    }
  }
  out.model_matrix.assign(
    model_values,
    model_values + static_cast<size_t>(n * out.n_columns)
  );

  SEXP group_map = list_element(factor, "group_map");
  if (TYPEOF(group_map) != INTSXP || XLENGTH(group_map) != n) {
    Rf_error("Random covariance factor group mapping is invalid.");
  }
  out.group_map.assign(INTEGER(group_map), INTEGER(group_map) + n);
  for (int group : out.group_map) {
    if (group < 1) {
      Rf_error("Random covariance factor group mapping is invalid.");
    }
  }

  out.coefficient_structure = DENSE_COEFFICIENT;
  SEXP coefficient_structure = optional_list_element(
    factor,
    "coefficient_structure"
  );
  if (coefficient_structure != R_NilValue) {
    if (TYPEOF(coefficient_structure) != STRSXP ||
        XLENGTH(coefficient_structure) != 1) {
      Rf_error("Random covariance coefficient structure must be one string.");
    }
    const char *structure_value = CHAR(STRING_ELT(coefficient_structure, 0));
    if (std::strcmp(structure_value, "diagonal") == 0) {
      out.coefficient_structure = DIAGONAL_COEFFICIENT;
    } else if (std::strcmp(structure_value, "markov") == 0) {
      out.coefficient_structure = MARKOV_COEFFICIENT;
    } else if (std::strcmp(structure_value, "dense") != 0) {
      Rf_error("Unknown random covariance coefficient structure.");
    }
  }
  SEXP initial_coefficient_factor = optional_list_element(
    factor,
    "coefficient_factor"
  );
  if (initial_coefficient_factor != R_NilValue) {
    require_real_square_matrix(
      initial_coefficient_factor,
      out.n_columns,
      "factor$coefficient_factor"
    );
    if (out.coefficient_structure == DIAGONAL_COEFFICIENT) {
      require_diagonal_coefficient_factor(
        REAL(initial_coefficient_factor),
        out.n_columns,
        "factor$coefficient_factor"
      );
    }
  }

  if (std::strcmp(type_value, "group") == 0) {
    out.type = GROUP_FACTOR;
    return out;
  }
  if (std::strcmp(type_value, "row_group") == 0) {
    out.type = ROW_GROUP_FACTOR;
    return out;
  }
  if (std::strcmp(type_value, "known_group") != 0) {
    Rf_error("Unknown random covariance factor type.");
  }

  SEXP group_covariance = list_element(factor, "group_covariance");
  if (TYPEOF(group_covariance) != REALSXP ||
      !Rf_isMatrix(group_covariance)) {
    Rf_error("Known group covariance factor must be a double matrix.");
  }
  SEXP group_dimensions = Rf_getAttrib(group_covariance, R_DimSymbol);
  out.n_groups = INTEGER(group_dimensions)[0];
  if (out.n_groups < 1 ||
      INTEGER(group_dimensions)[1] != out.n_groups) {
    Rf_error("Known group covariance factor must be square.");
  }
  for (int group : out.group_map) {
    if (group > out.n_groups) {
      Rf_error("Known group covariance factor mapping is out of bounds.");
    }
  }
  out.type = KNOWN_GROUP_FACTOR;
  const double *group_values = REAL(group_covariance);
  for (int i = 0; i < out.n_groups * out.n_groups; ++i) {
    if (!std::isfinite(group_values[i])) {
      Rf_error("Known group covariance factor must be finite.");
    }
  }
  for (int column = 0; column < out.n_groups; ++column) {
    for (int row = 0; row < out.n_groups; ++row) {
      if (group_values[row + out.n_groups * column] !=
          group_values[column + out.n_groups * row]) {
        Rf_error("Known group covariance factor must be exactly symmetric.");
      }
    }
  }
  {
    out.group_covariance.assign(
      group_values,
      group_values + static_cast<size_t>(out.n_groups * out.n_groups)
    );
    out.group_factor = out.group_covariance;
    int info = 0;
    F77_CALL(dpotrf)(
      "L",
      &out.n_groups,
      out.group_factor.data(),
      &out.n_groups,
      &info FCONE
    );
    if (info != 0) {
      Rf_error("Known group covariance factor must be positive definite.");
    }
    for (int column = 0; column < out.n_groups; ++column) {
      for (int row = 0; row < column; ++row) {
        out.group_factor[static_cast<size_t>(
          row + out.n_groups * column
        )] = 0.0;
      }
    }
  }
  return out;
}

void compile_fixed_known_group_basis(
    BlockPlan &block,
    const double *sampling,
    int n,
    const std::vector<PlanFactor> &factors)
{
  const int size = static_cast<int>(block.index.size());
  const int rank = block.rank;
  const PlanFactor &factor = factors.front();
  block.fixed_known_group_basis.assign(
    static_cast<size_t>(size * rank),
    0.0
  );
  for (const LowRankGroup &group : block.low_rank_groups) {
    for (int local : group.local_rows) {
      const int global = block.index[static_cast<size_t>(local)];
      const int group_index = factor.group_map[static_cast<size_t>(global)] - 1;
      const double group_loading = factor.group_factor[static_cast<size_t>(
        group_index + factor.n_groups * group.group_factor_column
      )];
      for (int coefficient = 0;
           coefficient < factor.n_columns;
           ++coefficient) {
        block.fixed_known_group_basis[static_cast<size_t>(
          local + size * (group.column_offset + coefficient)
        )] = factor.model_matrix[static_cast<size_t>(
          global + n * coefficient
        )] * group_loading;
      }
    }
  }

  block.fixed_known_group_precision.resize(static_cast<size_t>(size));
  block.fixed_known_group_log_determinant = 0.0;
  std::vector<double> weighted_basis = block.fixed_known_group_basis;
  for (int local = 0; local < size; ++local) {
    const int global = block.index[static_cast<size_t>(local)];
    const double variance = sampling[global + n * global];
    if (!(variance > 0.0) || !std::isfinite(variance)) {
      block.fixed_known_group_basis.clear();
      block.fixed_known_group_precision.clear();
      return;
    }
    block.fixed_known_group_precision[static_cast<size_t>(local)] =
      1.0 / variance;
    block.fixed_known_group_log_determinant += std::log(variance);
    const double inverse_root = 1.0 / std::sqrt(variance);
    for (int column = 0; column < rank; ++column) {
      weighted_basis[static_cast<size_t>(
        local + size * column
      )] *= inverse_root;
    }
  }

  block.fixed_known_group_crossproduct.resize(
    static_cast<size_t>(rank * rank)
  );
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dsyrk)(
    "L", "T", &rank, &size, &one, weighted_basis.data(), &size,
    &zero, block.fixed_known_group_crossproduct.data(), &rank FCONE FCONE
  );
  block.fixed_known_group_eligible = true;
}

BlockPlan make_block_plan(SEXP index, const double *sampling, int n,
                          const std::vector<PlanFactor> &factors)
{
  if (TYPEOF(index) != INTSXP || XLENGTH(index) == 0) {
    Rf_error("Each covariance block index must be a non-empty integer vector.");
  }
  BlockPlan out = {};
  const int size = static_cast<int>(XLENGTH(index));
  const int *values = INTEGER(index);
  out.index.resize(static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) {
    const int global = values[i] - 1;
    if (global < 0 || global >= n) {
      Rf_error("Covariance block indices are out of bounds.");
    }
    out.index[static_cast<size_t>(i)] = global;
  }

  bool sampling_diagonal = true;
  for (int col = 0; col < size; ++col) {
    for (int row = 0; row < size; ++row) {
      if (row != col && sampling[
          out.index[static_cast<size_t>(row)] +
          n * out.index[static_cast<size_t>(col)]
        ] != 0.0) {
        sampling_diagonal = false;
      }
    }
  }

  bool compatible = true;
  int rank = 0;
  for (size_t factor_i = 0; factor_i < factors.size(); ++factor_i) {
    const PlanFactor &factor = factors[factor_i];
    if (factor.type == DENSE_FACTOR) {
      compatible = false;
      continue;
    }
    if (factor.type == KNOWN_GROUP_FACTOR) {
      for (int root_column = 0;
           root_column < factor.n_groups;
           ++root_column) {
        LowRankGroup low_rank_group = {};
        low_rank_group.factor = static_cast<int>(factor_i);
        low_rank_group.column_offset = rank;
        low_rank_group.group_factor_column = root_column;
        for (int local = 0; local < size; ++local) {
          const int global = out.index[static_cast<size_t>(local)];
          const int group = factor.group_map[static_cast<size_t>(global)] - 1;
          if (group < 0 || group >= factor.n_groups) {
            Rf_error("Known group covariance factor mapping is out of bounds.");
          }
          if (factor.group_factor[static_cast<size_t>(
                group + factor.n_groups * root_column
              )] != 0.0) {
            low_rank_group.local_rows.push_back(local);
          }
        }
        if (!low_rank_group.local_rows.empty()) {
          out.low_rank_groups.push_back(low_rank_group);
          rank += factor.n_columns;
        }
      }
      continue;
    }
    std::vector<int> groups;
    for (int local = 0; local < size; ++local) {
      const int group = factor.group_map[
        static_cast<size_t>(out.index[static_cast<size_t>(local)])
      ];
      if (group < 1) {
        Rf_error("Random covariance factor group mapping is invalid.");
      }
      if (std::find(groups.begin(), groups.end(), group) == groups.end()) {
        groups.push_back(group);
      }
    }
    for (int group : groups) {
      LowRankGroup low_rank_group = {};
      low_rank_group.factor = static_cast<int>(factor_i);
      low_rank_group.column_offset = rank;
      low_rank_group.group_factor_column = -1;
      for (int local = 0; local < size; ++local) {
        const int global = out.index[static_cast<size_t>(local)];
        if (factor.group_map[static_cast<size_t>(global)] == group) {
          low_rank_group.local_rows.push_back(local);
        }
      }
      out.low_rank_groups.push_back(low_rank_group);
      rank += factor.n_columns;
    }
  }
  out.rank = rank;
  out.active_columns.resize(static_cast<size_t>(size));
  for (const LowRankGroup &group : out.low_rank_groups) {
    const int q = factors[static_cast<size_t>(group.factor)].n_columns;
    for (int local : group.local_rows) {
      std::vector<int> &columns = out.active_columns[
        static_cast<size_t>(local)
      ];
      for (int root_column = 0; root_column < q; ++root_column) {
        columns.push_back(group.column_offset + root_column);
      }
    }
  }
  for (std::vector<int> &columns : out.active_columns) {
    std::sort(columns.begin(), columns.end());
    columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
  }
  out.markov_eligible = false;
  out.fixed_known_group_eligible = false;
  out.markov_factor = -1;
  if (sampling_diagonal && factors.size() == 1 &&
      out.low_rank_groups.size() == 1) {
    const LowRankGroup &group = out.low_rank_groups.front();
    const PlanFactor &factor = factors[static_cast<size_t>(group.factor)];
    if ((factor.type == GROUP_FACTOR || factor.type == ROW_GROUP_FACTOR) &&
        factor.coefficient_structure == MARKOV_COEFFICIENT &&
        group.local_rows.size() == static_cast<size_t>(size)) {
      out.markov_rows_by_level.resize(
        static_cast<size_t>(factor.n_columns)
      );
      out.markov_loading.resize(static_cast<size_t>(size));
      bool one_state_per_row = true;
      for (int local = 0; local < size; ++local) {
        const int global = out.index[static_cast<size_t>(local)];
        int level = -1;
        double loading = 0.0;
        for (int column = 0; column < factor.n_columns; ++column) {
          const double value = factor.model_matrix[static_cast<size_t>(
            global + n * column
          )];
          if (value != 0.0) {
            if (level != -1) {
              one_state_per_row = false;
              break;
            }
            level = column;
            loading = value;
          }
        }
        if (!one_state_per_row || level == -1) {
          one_state_per_row = false;
          break;
        }
        out.markov_rows_by_level[static_cast<size_t>(level)].push_back(local);
        out.markov_loading[static_cast<size_t>(local)] = loading;
      }
      if (one_state_per_row) {
        out.markov_factor = group.factor;
        out.markov_eligible = true;
      } else {
        out.markov_rows_by_level.clear();
        out.markov_loading.clear();
      }
    }
  }
  out.root_eligible = compatible;
  out.low_rank_eligible = sampling_diagonal && compatible && rank > 0 &&
    rank < size;
  out.sparse_assembly_eligible = out.low_rank_eligible;
  out.sparse_factor_candidate = false;
  out.sparse_factor_eligible = false;
  out.sparse_factor_block_base = false;
  out.sparse_factor = nullptr;
  if (!sampling_diagonal && compatible && rank > 0) {
    std::vector<int> component(static_cast<size_t>(size), -1);
    std::vector<std::vector<int>> sampling_components;
    for (int seed = 0; seed < size; ++seed) {
      if (component[static_cast<size_t>(seed)] != -1) {
        continue;
      }
      const int component_id = static_cast<int>(sampling_components.size());
      std::vector<int> pending(1, seed);
      component[static_cast<size_t>(seed)] = component_id;
      std::vector<int> local_rows;
      while (!pending.empty()) {
        const int current_local = pending.back();
        pending.pop_back();
        local_rows.push_back(current_local);
        const int current_global = out.index[
          static_cast<size_t>(current_local)
        ];
        for (int candidate = 0; candidate < size; ++candidate) {
          if (component[static_cast<size_t>(candidate)] != -1) {
            continue;
          }
          const int candidate_global = out.index[
            static_cast<size_t>(candidate)
          ];
          if (sampling[current_global + n * candidate_global] != 0.0) {
            component[static_cast<size_t>(candidate)] = component_id;
            pending.push_back(candidate);
          }
        }
      }
      std::sort(local_rows.begin(), local_rows.end());
      sampling_components.push_back(local_rows);
    }
    if (sampling_components.size() > 1) {
      out.sampling_blocks.reserve(sampling_components.size());
      for (const std::vector<int> &local_rows : sampling_components) {
        SamplingBlock sampling_block = {};
        sampling_block.local_rows = local_rows;
        const int block_size = static_cast<int>(local_rows.size());
        sampling_block.covariance.resize(
          static_cast<size_t>(block_size * block_size)
        );
        for (int column = 0; column < block_size; ++column) {
          const int global_column = out.index[static_cast<size_t>(
            local_rows[static_cast<size_t>(column)]
          )];
          for (int row = 0; row < block_size; ++row) {
            const int global_row = out.index[static_cast<size_t>(
              local_rows[static_cast<size_t>(row)]
            )];
            sampling_block.covariance[static_cast<size_t>(
              row + block_size * column
            )] = sampling[global_row + n * global_column];
          }
        }
        sampling_block.factor = sampling_block.covariance;
        int info = 0;
        F77_CALL(dpotrf)(
          "L",
          &block_size,
          sampling_block.factor.data(),
          &block_size,
          &info FCONE
        );
        sampling_block.factor_eligible = info == 0;
        if (sampling_block.factor_eligible) {
          sampling_block.log_determinant = 0.0;
          for (int row = 0; row < block_size; ++row) {
            sampling_block.log_determinant += 2.0 * std::log(
              sampling_block.factor[static_cast<size_t>(
                row + block_size * row
              )]
            );
          }
          std::vector<double> precision = sampling_block.factor;
          F77_CALL(dpotri)(
            "L", &block_size, precision.data(), &block_size, &info FCONE
          );
          sampling_block.factor_eligible = info == 0;
          if (sampling_block.factor_eligible) {
            sampling_block.precision_diagonal.resize(
              static_cast<size_t>(block_size)
            );
            for (int row = 0; row < block_size; ++row) {
              sampling_block.precision_diagonal[static_cast<size_t>(row)] =
                precision[static_cast<size_t>(row + block_size * row)];
            }
          }
        }
        if (!sampling_block.factor_eligible) {
          sampling_block.factor.clear();
          sampling_block.precision_diagonal.clear();
        }
        out.sampling_blocks.push_back(sampling_block);
      }
    }
  }
  out.block_base_eligible = out.sampling_blocks.size() > 1 && rank < size;
  out.spectral_eligible = !sampling_diagonal && compatible && rank < size &&
    !out.block_base_eligible;
  if (out.spectral_eligible) {
    out.sampling_eigenvectors.resize(static_cast<size_t>(size * size));
    out.sampling_eigenvalues.resize(static_cast<size_t>(size));
    for (int column = 0; column < size; ++column) {
      const int global_column = out.index[static_cast<size_t>(column)];
      for (int row = 0; row < size; ++row) {
        const int global_row = out.index[static_cast<size_t>(row)];
        out.sampling_eigenvectors[static_cast<size_t>(row + size * column)] =
          sampling[global_row + n * global_column];
      }
    }
    const int lwork = std::max(1, 3 * size - 1);
    std::vector<double> work(static_cast<size_t>(lwork));
    int info = 0;
    F77_CALL(dsyev)(
      "V",
      "L",
      &size,
      out.sampling_eigenvectors.data(),
      &size,
      out.sampling_eigenvalues.data(),
      work.data(),
      &lwork,
      &info FCONE FCONE
    );
    if (info != 0) {
      Rf_error("Known-V sampling covariance eigendecomposition failed.");
    }
  }
  if (sampling_diagonal && factors.size() == 1 && rank > 0 && rank < size &&
      factors.front().type == KNOWN_GROUP_FACTOR &&
      factors.front().coefficient_structure == DIAGONAL_COEFFICIENT) {
    compile_fixed_known_group_basis(out, sampling, n, factors);
  }
  if (sampling_diagonal && compatible && rank > 0) {
    out.precision_active_columns = out.active_columns;
    compile_sparse_latent_pattern(out, out.active_columns);
  } else if (out.sampling_blocks.size() > 1) {
    std::vector<std::vector<int>> transformed_active_columns(
      static_cast<size_t>(size)
    );
    out.precision_active_columns.resize(static_cast<size_t>(size));
    for (const SamplingBlock &sampling_block : out.sampling_blocks) {
      std::vector<int> block_columns;
      for (int local : sampling_block.local_rows) {
        block_columns.insert(
          block_columns.end(),
          out.active_columns[static_cast<size_t>(local)].begin(),
          out.active_columns[static_cast<size_t>(local)].end()
        );
      }
      std::sort(block_columns.begin(), block_columns.end());
      block_columns.erase(
        std::unique(block_columns.begin(), block_columns.end()),
        block_columns.end()
      );
      std::vector<int> cumulative_columns;
      for (int local : sampling_block.local_rows) {
        cumulative_columns.insert(
          cumulative_columns.end(),
          out.active_columns[static_cast<size_t>(local)].begin(),
          out.active_columns[static_cast<size_t>(local)].end()
        );
        std::sort(cumulative_columns.begin(), cumulative_columns.end());
        cumulative_columns.erase(
          std::unique(cumulative_columns.begin(), cumulative_columns.end()),
          cumulative_columns.end()
        );
        transformed_active_columns[static_cast<size_t>(local)] =
          cumulative_columns;
        out.precision_active_columns[static_cast<size_t>(local)] =
          block_columns;
      }
    }
    out.sparse_factor_block_base = true;
    compile_sparse_latent_pattern(out, transformed_active_columns);
  }
  if (!out.sparse_values.empty()) {
    const size_t dense_entries = static_cast<size_t>(rank) *
      static_cast<size_t>(rank + 1) / 2;
    out.sparse_factor_candidate =
      out.sparse_row_indices.size() < dense_entries;
    if (!out.sparse_factor_candidate) {
      out.active_pair_positions.clear();
      out.sparse_column_pointers.clear();
      out.sparse_row_indices.clear();
      out.sparse_diagonal_positions.clear();
      out.sparse_values.clear();
    }
  }
  return out;
}

CovariancePlan *make_plan(SEXP y, SEXP sampling_covariance,
                          SEXP random_covariance_factors,
                          SEXP block_indices)
{
  const R_xlen_t n_xlen = XLENGTH(y);
  if (n_xlen > static_cast<R_xlen_t>(2147483647)) {
    Rf_error("Known-V covariance plan exceeds integer indexing limits.");
  }
  const int n = static_cast<int>(n_xlen);
  require_real_vector(y, n_xlen, "y");
  require_real_square_matrix(sampling_covariance, n, "sampling_covariance");
  if (TYPEOF(block_indices) != VECSXP || XLENGTH(block_indices) == 0) {
    Rf_error("Covariance block indices must be a non-empty list.");
  }
  if (TYPEOF(random_covariance_factors) != VECSXP) {
    Rf_error("Random covariance factor plans must be a list.");
  }

  CovariancePlan *plan = new CovariancePlan();
  plan->n = n;
  plan->y.assign(REAL(y), REAL(y) + n);
  plan->sampling_covariance.assign(
    REAL(sampling_covariance),
    REAL(sampling_covariance) + static_cast<size_t>(n * n)
  );
  for (double value : plan->y) {
    if (!std::isfinite(value)) {
      delete plan;
      Rf_error("Known-V bridge observations must be finite.");
    }
  }
  for (double value : plan->sampling_covariance) {
    if (!std::isfinite(value)) {
      delete plan;
      Rf_error("Known-V sampling covariance must be finite.");
    }
  }
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      if (plan->sampling_covariance[static_cast<size_t>(row + n * col)] !=
          plan->sampling_covariance[static_cast<size_t>(col + n * row)]) {
        delete plan;
        Rf_error("Known-V sampling covariance must be exactly symmetric.");
      }
    }
  }

  plan->factors.reserve(static_cast<size_t>(
    XLENGTH(random_covariance_factors)
  ));
  for (R_xlen_t factor_i = 0;
       factor_i < XLENGTH(random_covariance_factors);
       ++factor_i) {
    plan->factors.push_back(covariance_factor_plan(
      VECTOR_ELT(random_covariance_factors, factor_i),
      n
    ));
  }
  size_t max_block_size = 0;
  size_t max_rank = 0;
  size_t max_dense_rank = 0;
  size_t max_base_rhs_size = 0;
  std::vector<int> membership(static_cast<size_t>(n), -1);
  for (R_xlen_t block_i = 0; block_i < XLENGTH(block_indices); ++block_i) {
    BlockPlan block = make_block_plan(
      VECTOR_ELT(block_indices, block_i),
      plan->sampling_covariance.data(),
      n,
      plan->factors
    );
    for (int global : block.index) {
      if (membership[static_cast<size_t>(global)] != -1) {
        delete plan;
        Rf_error("Covariance block indices must not overlap.");
      }
      membership[static_cast<size_t>(global)] = static_cast<int>(block_i);
    }
    max_block_size = std::max(max_block_size, block.index.size());
    max_rank = std::max(max_rank, static_cast<size_t>(block.rank));
    if (block.low_rank_eligible || block.block_base_eligible ||
        block.spectral_eligible) {
      max_dense_rank = std::max(
        max_dense_rank,
        static_cast<size_t>(block.rank)
      );
    }
    for (const SamplingBlock &sampling_block : block.sampling_blocks) {
      max_base_rhs_size = std::max(
        max_base_rhs_size,
        sampling_block.local_rows.size() *
          static_cast<size_t>(block.rank + 1)
      );
    }
    plan->blocks.push_back(block);
  }
  if (std::find(membership.begin(), membership.end(), -1) != membership.end()) {
    delete plan;
    Rf_error("Covariance block indices must partition every observation.");
  }

  bool has_sparse_factor_candidate = false;
  for (const BlockPlan &block : plan->blocks) {
    has_sparse_factor_candidate = has_sparse_factor_candidate ||
      block.sparse_factor_candidate;
  }
  if (has_sparse_factor_candidate) {
    if (!M_cholmod_start(&plan->cholmod_common_state)) {
      delete plan;
      Rf_error("Known-V sparse latent solver could not be initialized.");
    }
    plan->cholmod_started = true;
    plan->cholmod_common_state.dbound = 0.0;
    plan->cholmod_common_state.final_ll = 1;
    for (BlockPlan &block : plan->blocks) {
      if (!block.sparse_factor_candidate) {
        continue;
      }
      cholmod_sparse matrix = sparse_matrix_view(block);
      block.sparse_factor = M_cholmod_analyze(
        &matrix,
        &plan->cholmod_common_state
      );
      block.sparse_factor_eligible = block.sparse_factor != nullptr;
    }
  }

  plan->covariance_work.resize(max_block_size * max_block_size);
  plan->residual_work.resize(max_block_size);
  plan->diagonal_work.resize(max_block_size);
  plan->low_rank_work.resize(max_block_size * max_rank);
  plan->small_matrix_work.resize(max_dense_rank * max_dense_rank);
  plan->rhs_work.resize(max_rank);
  plan->scale_work.resize(max_rank);
  plan->base_rhs_work.resize(max_base_rhs_size);
  return plan;
}

CovariancePlan *plan_pointer(SEXP pointer)
{
  if (TYPEOF(pointer) != EXTPTRSXP) {
    Rf_error("Known-V covariance plan must be an external pointer.");
  }
  CovariancePlan *plan = static_cast<CovariancePlan *>(
    R_ExternalPtrAddr(pointer)
  );
  if (plan == nullptr) {
    Rf_error("Known-V covariance plan is no longer valid.");
  }
  return plan;
}

std::vector<CovarianceFactor> covariance_states(
    SEXP states,
    const CovariancePlan &plan)
{
  if (TYPEOF(states) != VECSXP ||
      XLENGTH(states) != static_cast<R_xlen_t>(plan.factors.size())) {
    Rf_error("Known-V covariance plan received inconsistent factor states.");
  }

  std::vector<CovarianceFactor> out;
  out.reserve(plan.factors.size());
  for (size_t factor_i = 0; factor_i < plan.factors.size(); ++factor_i) {
    const PlanFactor &stored = plan.factors[factor_i];
    SEXP state = VECTOR_ELT(states, static_cast<R_xlen_t>(factor_i));
    CovarianceFactor value = {};
    value.type = stored.type;
    value.coefficient_structure = stored.coefficient_structure;
    value.n_columns = stored.n_columns;
    value.n_groups = stored.n_groups;

    if (stored.type == DENSE_FACTOR) {
      SEXP covariance = list_element(state, "covariance");
      require_real_square_matrix(
        covariance,
        plan.n,
        "factor_state$covariance"
      );
      value.covariance = REAL(covariance);
      out.push_back(value);
      continue;
    }

    SEXP coefficient_factor = list_element(state, "coefficient_factor");
    require_real_square_matrix(
      coefficient_factor,
      stored.n_columns,
      "factor_state$coefficient_factor"
    );
    const double *coefficient_values = REAL(coefficient_factor);
    for (int i = 0; i < stored.n_columns * stored.n_columns; ++i) {
      if (!std::isfinite(coefficient_values[i])) {
        Rf_error("Known-V coefficient factor state must be finite.");
      }
    }
    value.model_matrix = stored.model_matrix.data();
    value.group_map = stored.group_map.data();
    value.coefficient_factor = coefficient_values;

    if (stored.coefficient_structure == DIAGONAL_COEFFICIENT) {
      require_diagonal_coefficient_factor(
        coefficient_values,
        stored.n_columns,
        "factor_state$coefficient_factor"
      );
    }

    if (stored.coefficient_structure == MARKOV_COEFFICIENT) {
      SEXP coefficient_scale = list_element(state, "coefficient_scale");
      SEXP transition = list_element(state, "markov_transition");
      SEXP innovation = list_element(
        state,
        "markov_innovation_variance"
      );
      require_real_vector(
        coefficient_scale,
        stored.n_columns,
        "factor_state$coefficient_scale"
      );
      require_real_vector(
        transition,
        stored.n_columns - 1,
        "factor_state$markov_transition"
      );
      require_real_vector(
        innovation,
        stored.n_columns - 1,
        "factor_state$markov_innovation_variance"
      );
      value.coefficient_scale = REAL(coefficient_scale);
      value.markov_transition = REAL(transition);
      value.markov_innovation_variance = REAL(innovation);
      for (int column = 0; column < stored.n_columns; ++column) {
        if (!std::isfinite(value.coefficient_scale[column]) ||
            value.coefficient_scale[column] < 0.0) {
          Rf_error("Known-V Markov coefficient scales must be finite and non-negative.");
        }
      }
      for (int column = 0; column < stored.n_columns - 1; ++column) {
        if (!std::isfinite(value.markov_transition[column]) ||
            !std::isfinite(value.markov_innovation_variance[column]) ||
            !(value.markov_innovation_variance[column] > 0.0)) {
          Rf_error("Known-V Markov transition state is invalid.");
        }
      }
    }

    if (stored.type == ROW_GROUP_FACTOR) {
      SEXP row_scale = list_element(state, "row_scale");
      require_real_vector(row_scale, plan.n, "factor_state$row_scale");
      const double *row_scale_values = REAL(row_scale);
      for (int i = 0; i < plan.n; ++i) {
        if (!std::isfinite(row_scale_values[i]) || row_scale_values[i] < 0.0) {
          Rf_error("Known-V row-scale factor state must be finite and non-negative.");
        }
      }
      value.row_scale = row_scale_values;
    } else if (stored.type == KNOWN_GROUP_FACTOR) {
      value.group_covariance = stored.group_covariance.data();
    }
    out.push_back(value);
  }
  return out;
}

void finalize_plan(SEXP pointer)
{
  CovariancePlan *plan = static_cast<CovariancePlan *>(
    R_ExternalPtrAddr(pointer)
  );
  if (plan != nullptr) {
    delete plan;
    R_ClearExternalPtr(pointer);
  }
}
