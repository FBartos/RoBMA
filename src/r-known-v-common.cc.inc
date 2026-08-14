void require_real_vector(SEXP x, R_xlen_t length, const char *name)
{
  if (TYPEOF(x) != REALSXP || XLENGTH(x) != length) {
    Rf_error("'%s' must be a double vector of the expected length.", name);
  }
}

void require_real_square_matrix(SEXP x, int size, const char *name)
{
  if (TYPEOF(x) != REALSXP || !Rf_isMatrix(x)) {
    Rf_error("'%s' must be a double matrix.", name);
  }
  SEXP dimensions = Rf_getAttrib(x, R_DimSymbol);
  if (INTEGER(dimensions)[0] != size || INTEGER(dimensions)[1] != size) {
    Rf_error("'%s' has inconsistent dimensions.", name);
  }
}

SEXP list_element(SEXP x, const char *name)
{
  if (TYPEOF(x) != VECSXP) {
    Rf_error("Random covariance factors must be lists.");
  }
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  if (TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(x)) {
    Rf_error("Random covariance factors must be named lists.");
  }
  for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(x, i);
    }
  }
  Rf_error("Random covariance factor is missing '%s'.", name);
  return R_NilValue;
}

enum FactorType {
  DENSE_FACTOR,
  GROUP_FACTOR,
  ROW_GROUP_FACTOR,
  KNOWN_GROUP_FACTOR
};

enum CoefficientStructure {
  DENSE_COEFFICIENT,
  DIAGONAL_COEFFICIENT,
  MARKOV_COEFFICIENT
};

void require_diagonal_coefficient_factor(const double *factor, int size,
                                         const char *name)
{
  for (int column = 0; column < size; ++column) {
    for (int row = 0; row < size; ++row) {
      if (row != column && factor[row + size * column] != 0.0) {
        Rf_error("'%s' must be exactly diagonal.", name);
      }
    }
  }
}

struct CovarianceFactor {
  FactorType type;
  CoefficientStructure coefficient_structure;
  const double *covariance;
  const double *model_matrix;
  const int *group_map;
  const double *coefficient_factor;
  const double *group_covariance;
  const double *row_scale;
  const double *coefficient_scale;
  const double *markov_transition;
  const double *markov_innovation_variance;
  int n_columns;
  int n_groups;
};

struct PlanFactor {
  FactorType type;
  std::vector<double> model_matrix;
  std::vector<int> group_map;
  std::vector<double> group_covariance;
  std::vector<double> group_factor;
  CoefficientStructure coefficient_structure;
  int n_columns;
  int n_groups;
};

struct LowRankGroup {
  int factor;
  int column_offset;
  int group_factor_column;
  std::vector<int> local_rows;
};

struct SamplingBlock {
  std::vector<int> local_rows;
  std::vector<double> covariance;
  std::vector<double> factor;
  std::vector<double> precision_diagonal;
  double log_determinant;
  bool factor_eligible;
};

struct BlockPlan {
  std::vector<int> index;
  std::vector<LowRankGroup> low_rank_groups;
  std::vector<std::vector<int>> active_columns;
  std::vector<std::vector<int>> sparse_factor_active_columns;
  std::vector<std::vector<int>> precision_active_columns;
  std::vector<std::vector<int>> active_pair_positions;
  std::vector<int> sparse_column_pointers;
  std::vector<int> sparse_row_indices;
  std::vector<int> sparse_diagonal_positions;
  std::vector<double> sparse_values;
  std::vector<SamplingBlock> sampling_blocks;
  std::vector<double> sampling_eigenvalues;
  std::vector<double> sampling_eigenvectors;
  std::vector<std::vector<int>> markov_rows_by_level;
  std::vector<double> markov_loading;
  std::vector<double> fixed_known_group_basis;
  std::vector<double> fixed_known_group_crossproduct;
  std::vector<double> fixed_known_group_precision;
  double fixed_known_group_log_determinant;
  int markov_factor;
  int rank;
  bool fixed_known_group_eligible;
  bool markov_eligible;
  bool low_rank_eligible;
  bool sparse_assembly_eligible;
  bool sparse_factor_candidate;
  bool sparse_factor_eligible;
  bool sparse_factor_block_base;
  bool block_base_eligible;
  bool spectral_eligible;
  bool root_eligible;
  cholmod_factor *sparse_factor;
};

struct CovariancePlan {
  CovariancePlan() : n(0), cholmod_started(false) {}

  ~CovariancePlan()
  {
    if (cholmod_started) {
      for (BlockPlan &block : blocks) {
        if (block.sparse_factor != nullptr) {
          M_cholmod_free_factor(&block.sparse_factor, &cholmod_common_state);
        }
      }
      M_cholmod_finish(&cholmod_common_state);
    }
  }

  int n;
  bool cholmod_started;
  cholmod_common cholmod_common_state;
  std::vector<double> y;
  std::vector<double> sampling_covariance;
  std::vector<PlanFactor> factors;
  std::vector<BlockPlan> blocks;
  std::vector<double> covariance_work;
  std::vector<double> residual_work;
  std::vector<double> diagonal_work;
  std::vector<double> low_rank_work;
  std::vector<double> small_matrix_work;
  std::vector<double> rhs_work;
  std::vector<double> base_rhs_work;
  std::vector<double> conditional_rhs_work;
  std::vector<double> conditional_base_work;
  std::vector<double> unit_diagonal_work;
  std::vector<double> markov_diagonal_work;
  std::vector<double> markov_off_diagonal_work;
  std::vector<double> markov_factor_work;
  std::vector<double> markov_rhs_work;
  std::vector<double> markov_solution_work;
  std::vector<double> markov_variance_work;
  std::vector<double> scale_work;
};

cholmod_sparse sparse_matrix_view(BlockPlan &block)
{
  cholmod_sparse matrix = {};
  matrix.nrow = static_cast<size_t>(block.rank);
  matrix.ncol = static_cast<size_t>(block.rank);
  matrix.nzmax = block.sparse_values.size();
  matrix.p = block.sparse_column_pointers.data();
  matrix.i = block.sparse_row_indices.data();
  matrix.nz = nullptr;
  matrix.x = block.sparse_values.data();
  matrix.z = nullptr;
  matrix.stype = -1;
  matrix.itype = CHOLMOD_INT;
  matrix.xtype = CHOLMOD_REAL;
  matrix.dtype = CHOLMOD_DOUBLE;
  matrix.sorted = 1;
  matrix.packed = 1;
  return matrix;
}

SEXP optional_list_element(SEXP x, const char *name)
{
  if (TYPEOF(x) != VECSXP) {
    Rf_error("Random covariance factors must be lists.");
  }
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  if (TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(x)) {
    Rf_error("Random covariance factors must be named lists.");
  }
  for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(x, i);
    }
  }
  return R_NilValue;
}

void compile_sparse_latent_pattern(
    BlockPlan &block,
    const std::vector<std::vector<int>> &active_columns)
{
  const int rank = block.rank;
  block.sparse_factor_active_columns = active_columns;
  std::vector<std::vector<int>> rows_by_column(
    static_cast<size_t>(rank)
  );
  for (int column = 0; column < rank; ++column) {
    rows_by_column[static_cast<size_t>(column)].push_back(column);
  }
  for (const std::vector<int> &columns : active_columns) {
    for (size_t column_i = 0; column_i < columns.size(); ++column_i) {
      std::vector<int> &rows = rows_by_column[static_cast<size_t>(
        columns[column_i]
      )];
      for (size_t row_i = column_i; row_i < columns.size(); ++row_i) {
        rows.push_back(columns[row_i]);
      }
    }
  }

  block.sparse_column_pointers.resize(static_cast<size_t>(rank + 1));
  block.sparse_diagonal_positions.resize(static_cast<size_t>(rank));
  for (int column = 0; column < rank; ++column) {
    std::vector<int> &rows = rows_by_column[static_cast<size_t>(column)];
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    block.sparse_column_pointers[static_cast<size_t>(column)] =
      static_cast<int>(block.sparse_row_indices.size());
    block.sparse_row_indices.insert(
      block.sparse_row_indices.end(),
      rows.begin(),
      rows.end()
    );
    block.sparse_diagonal_positions[static_cast<size_t>(column)] =
      block.sparse_column_pointers[static_cast<size_t>(column)];
  }
  block.sparse_column_pointers[static_cast<size_t>(rank)] =
    static_cast<int>(block.sparse_row_indices.size());
  block.sparse_values.resize(block.sparse_row_indices.size());

  block.active_pair_positions.resize(active_columns.size());
  for (size_t observation = 0;
       observation < active_columns.size();
       ++observation) {
    const std::vector<int> &columns = active_columns[observation];
    std::vector<int> &positions = block.active_pair_positions[observation];
    positions.reserve(columns.size() * (columns.size() + 1) / 2);
    for (size_t column_i = 0; column_i < columns.size(); ++column_i) {
      const int column = columns[column_i];
      const int begin = block.sparse_column_pointers[
        static_cast<size_t>(column)
      ];
      const int end = block.sparse_column_pointers[
        static_cast<size_t>(column + 1)
      ];
      for (size_t row_i = column_i; row_i < columns.size(); ++row_i) {
        const int row = columns[row_i];
        const std::vector<int>::const_iterator found = std::lower_bound(
          block.sparse_row_indices.begin() + begin,
          block.sparse_row_indices.begin() + end,
          row
        );
        if (found == block.sparse_row_indices.begin() + end ||
            *found != row) {
          Rf_error("Known-V sparse latent pattern is inconsistent.");
        }
        positions.push_back(static_cast<int>(
          found - block.sparse_row_indices.begin()
        ));
      }
    }
  }
}

double design_covariance(const CovarianceFactor &factor, int row, int col,
                         int n)
{
  const int q = factor.n_columns;
  double value = 0.0;
  if (factor.coefficient_structure == DIAGONAL_COEFFICIENT) {
    for (int coefficient = 0; coefficient < q; ++coefficient) {
      const double scale = factor.coefficient_factor[
        coefficient + q * coefficient
      ];
      const double row_value =
        factor.model_matrix[row + n * coefficient] * scale;
      const double col_value =
        factor.model_matrix[col + n * coefficient] * scale;
      value += row_value * col_value;
    }
    if (factor.type == ROW_GROUP_FACTOR) {
      value *= factor.row_scale[row] * factor.row_scale[col];
    }
    return value;
  }
  for (int root = 0; root < q; ++root) {
    double row_value = 0.0;
    double col_value = 0.0;
    for (int coefficient = 0; coefficient < q; ++coefficient) {
      const double loading = factor.coefficient_factor[
        coefficient + q * root
      ];
      row_value += factor.model_matrix[row + n * coefficient] * loading;
      col_value += factor.model_matrix[col + n * coefficient] * loading;
    }
    value += row_value * col_value;
  }
  if (factor.type == ROW_GROUP_FACTOR) {
    value *= factor.row_scale[row] * factor.row_scale[col];
  }
  return value;
}

double covariance_contribution(const CovarianceFactor &factor,
                               int row, int col, int n)
{
  if (factor.type == DENSE_FACTOR) {
    return factor.covariance[row + n * col];
  }
  const int row_group = factor.group_map[row] - 1;
  const int col_group = factor.group_map[col] - 1;
  if (row_group < 0 || col_group < 0) {
    Rf_error("Random covariance factor group mapping is invalid.");
  }
  if (factor.type == KNOWN_GROUP_FACTOR) {
    if (row_group >= factor.n_groups || col_group >= factor.n_groups) {
      Rf_error("Known group covariance factor mapping is out of bounds.");
    }
    return factor.group_covariance[
      row_group + factor.n_groups * col_group
    ] * design_covariance(factor, row, col, n);
  }
  if (row_group != col_group) {
    return 0.0;
  }
  return design_covariance(factor, row, col, n);
}
