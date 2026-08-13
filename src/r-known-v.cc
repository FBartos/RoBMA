#include <R_ext/Error.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <Rinternals.h>
#include <Matrix/Matrix.h>
#include <Matrix/stubs.c>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

#ifndef FCONE
# define FCONE
#endif

namespace {

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
  int markov_factor;
  int rank;
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

std::vector<CovarianceFactor> covariance_factors(SEXP factors, int n)
{
  if (TYPEOF(factors) != VECSXP) {
    Rf_error("Random covariance factors must be a list.");
  }
  std::vector<CovarianceFactor> out;
  out.reserve(static_cast<size_t>(XLENGTH(factors)));
  for (R_xlen_t i = 0; i < XLENGTH(factors); ++i) {
    SEXP factor = VECTOR_ELT(factors, i);
    SEXP type = list_element(factor, "type");
    if (TYPEOF(type) != STRSXP || XLENGTH(type) != 1) {
      Rf_error("Random covariance factor type must be one string.");
    }
    const char *type_value = CHAR(STRING_ELT(type, 0));
    CovarianceFactor value = {};
    if (std::strcmp(type_value, "dense") == 0) {
      SEXP covariance = list_element(factor, "covariance");
      require_real_square_matrix(covariance, n, "factor$covariance");
      value.type = DENSE_FACTOR;
      value.covariance = REAL(covariance);
      out.push_back(value);
      continue;
    }

    SEXP model_matrix = list_element(factor, "model_matrix");
    if (TYPEOF(model_matrix) != REALSXP || !Rf_isMatrix(model_matrix)) {
      Rf_error("Random covariance factor model matrix must be double.");
    }
    SEXP model_dimensions = Rf_getAttrib(model_matrix, R_DimSymbol);
    const int n_columns = INTEGER(model_dimensions)[1];
    if (INTEGER(model_dimensions)[0] != n || n_columns < 1) {
      Rf_error("Random covariance factor model matrix has inconsistent dimensions.");
    }
    SEXP group_map = list_element(factor, "group_map");
    if (TYPEOF(group_map) != INTSXP || XLENGTH(group_map) != n) {
      Rf_error("Random covariance factor group mapping is invalid.");
    }
    SEXP coefficient_factor = list_element(factor, "coefficient_factor");
    require_real_square_matrix(
      coefficient_factor,
      n_columns,
      "factor$coefficient_factor"
    );
    value.model_matrix = REAL(model_matrix);
    value.group_map = INTEGER(group_map);
    value.coefficient_factor = REAL(coefficient_factor);
    value.n_columns = n_columns;
    value.coefficient_structure = DENSE_COEFFICIENT;
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
        value.coefficient_structure = DIAGONAL_COEFFICIENT;
      } else if (std::strcmp(structure_value, "markov") == 0) {
        value.coefficient_structure = MARKOV_COEFFICIENT;
      } else if (std::strcmp(structure_value, "dense") != 0) {
        Rf_error("Unknown random covariance coefficient structure.");
      }
    }
    if (value.coefficient_structure == MARKOV_COEFFICIENT) {
      SEXP coefficient_scale = list_element(factor, "coefficient_scale");
      SEXP transition = list_element(factor, "markov_transition");
      SEXP innovation = list_element(
        factor,
        "markov_innovation_variance"
      );
      require_real_vector(
        coefficient_scale,
        n_columns,
        "factor$coefficient_scale"
      );
      require_real_vector(
        transition,
        n_columns - 1,
        "factor$markov_transition"
      );
      require_real_vector(
        innovation,
        n_columns - 1,
        "factor$markov_innovation_variance"
      );
      value.coefficient_scale = REAL(coefficient_scale);
      value.markov_transition = REAL(transition);
      value.markov_innovation_variance = REAL(innovation);
      for (int column = 0; column < n_columns; ++column) {
        if (!std::isfinite(value.coefficient_scale[column]) ||
            value.coefficient_scale[column] < 0.0) {
          Rf_error("Random covariance Markov coefficient scales must be finite and non-negative.");
        }
      }
      for (int column = 0; column < n_columns - 1; ++column) {
        if (!std::isfinite(value.markov_transition[column]) ||
            !std::isfinite(value.markov_innovation_variance[column]) ||
            !(value.markov_innovation_variance[column] > 0.0)) {
          Rf_error("Random covariance Markov transition state is invalid.");
        }
      }
    } else if (value.coefficient_structure == DIAGONAL_COEFFICIENT) {
      require_diagonal_coefficient_factor(
        value.coefficient_factor,
        value.n_columns,
        "factor$coefficient_factor"
      );
    }

    if (std::strcmp(type_value, "group") == 0) {
      value.type = GROUP_FACTOR;
      out.push_back(value);
      continue;
    }
    if (std::strcmp(type_value, "row_group") == 0) {
      SEXP row_scale = list_element(factor, "row_scale");
      require_real_vector(row_scale, n, "factor$row_scale");
      value.type = ROW_GROUP_FACTOR;
      value.row_scale = REAL(row_scale);
      out.push_back(value);
      continue;
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
    const int n_groups = INTEGER(group_dimensions)[0];
    if (n_groups < 1 || INTEGER(group_dimensions)[1] != n_groups) {
      Rf_error("Known group covariance factor must be square.");
    }
    value.type = KNOWN_GROUP_FACTOR;
    value.group_covariance = REAL(group_covariance);
    value.n_groups = n_groups;
    out.push_back(value);
  }
  return out;
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

PlanFactor make_plan_factor(const CovarianceFactor &factor, int n)
{
  PlanFactor out = {};
  out.type = factor.type;
  out.coefficient_structure = factor.coefficient_structure;
  out.n_columns = factor.n_columns;
  out.n_groups = factor.n_groups;
  if (factor.type == DENSE_FACTOR) {
    return out;
  }
  out.model_matrix.assign(
    factor.model_matrix,
    factor.model_matrix + static_cast<size_t>(n * factor.n_columns)
  );
  out.group_map.assign(factor.group_map, factor.group_map + n);
  if (factor.type == KNOWN_GROUP_FACTOR) {
    out.group_covariance.assign(
      factor.group_covariance,
      factor.group_covariance +
        static_cast<size_t>(factor.n_groups * factor.n_groups)
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
  out.low_rank_eligible = sampling_diagonal && compatible && rank < size;
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
  std::vector<CovarianceFactor> current = covariance_factors(
    random_covariance_factors,
    n
  );

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

  plan->factors.reserve(current.size());
  for (const CovarianceFactor &factor : current) {
    plan->factors.push_back(make_plan_factor(factor, n));
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

void fill_dense_block(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *covariance,
    double *residual)
{
  const int size = static_cast<int>(block.index.size());
  for (int col = 0; col < size; ++col) {
    const int global_col = block.index[static_cast<size_t>(col)];
    for (int row = 0; row < size; ++row) {
      const int global_row = block.index[static_cast<size_t>(row)];
      double value = plan.sampling_covariance[
        static_cast<size_t>(global_row + plan.n * global_col)
      ];
      for (const CovarianceFactor &factor : factors) {
        value += covariance_contribution(
          factor,
          global_row,
          global_col,
          plan.n
        );
      }
      if (row == col) {
        value += extra[global_row];
      }
      if (!std::isfinite(value)) {
        Rf_error("Known-V bridge covariance must be finite.");
      }
      covariance[static_cast<size_t>(row + size * col)] = value;
    }
    if (!std::isfinite(mean[global_col])) {
      Rf_error("Known-V bridge means must be finite.");
    }
    residual[static_cast<size_t>(col)] =
      plan.y[static_cast<size_t>(global_col)] - mean[global_col];
  }
}

double cholesky_block_log_likelihood(int size, double *covariance,
                                     double *residual)
{
  int info = 0;
  F77_CALL(dpotrf)("L", &size, covariance, &size, &info FCONE);
  if (info != 0) {
    Rf_error("Known-V integrated bridge covariance is not positive definite.");
  }
  double log_determinant = 0.0;
  for (int row = 0; row < size; ++row) {
    const double diagonal = covariance[static_cast<size_t>(row + size * row)];
    log_determinant += 2.0 * std::log(diagonal);
    double value = residual[static_cast<size_t>(row)];
    for (int col = 0; col < row; ++col) {
      value -= covariance[static_cast<size_t>(row + size * col)] *
        residual[static_cast<size_t>(col)];
    }
    residual[static_cast<size_t>(row)] = value / diagonal;
  }
  double quadratic = 0.0;
  for (int row = 0; row < size; ++row) {
    const double value = residual[static_cast<size_t>(row)];
    quadratic += value * value;
  }
  const double log_two_pi = 1.837877066409345483560659472811;
  return -0.5 * (
    static_cast<double>(size) * log_two_pi +
    log_determinant +
    quadratic
  );
}

double dense_block_log_likelihood(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra)
{
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
  return cholesky_block_log_likelihood(size, covariance, residual);
}

double *fill_low_rank_factor(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors)
{
  const int size = static_cast<int>(block.index.size());
  const int rank = block.rank;
  double *U = plan.low_rank_work.data();
  if (rank == 0) {
    return U;
  }
  std::fill(U, U + static_cast<size_t>(size * rank), 0.0);
  for (const LowRankGroup &group : block.low_rank_groups) {
    const CovarianceFactor &factor = factors[
      static_cast<size_t>(group.factor)
    ];
    const int q = factor.n_columns;
    for (int local : group.local_rows) {
      const int global = block.index[static_cast<size_t>(local)];
      const double row_scale = factor.type == ROW_GROUP_FACTOR ?
        factor.row_scale[global] : 1.0;
      double group_scale = 1.0;
      if (factor.type == KNOWN_GROUP_FACTOR) {
        const PlanFactor &stored_factor = plan.factors[
          static_cast<size_t>(group.factor)
        ];
        const int group_index = factor.group_map[global] - 1;
        group_scale = stored_factor.group_factor[static_cast<size_t>(
          group_index + factor.n_groups * group.group_factor_column
        )];
      }
      for (int root_column = 0; root_column < q; ++root_column) {
        double transformed;
        if (factor.coefficient_structure == DIAGONAL_COEFFICIENT) {
          transformed = factor.model_matrix[
            global + plan.n * root_column
          ] * factor.coefficient_factor[
            root_column + q * root_column
          ];
        } else {
          transformed = 0.0;
          for (int coefficient = 0; coefficient < q; ++coefficient) {
            transformed += factor.model_matrix[
              global + plan.n * coefficient
            ] * factor.coefficient_factor[
              coefficient + q * root_column
            ];
          }
        }
        transformed *= row_scale * group_scale;
        if (!std::isfinite(transformed)) {
          Rf_error("Known-V bridge low-rank factor must be finite.");
        }
        U[static_cast<size_t>(
          local + size * (group.column_offset + root_column)
        )] = transformed;
      }
    }
  }
  return U;
}

bool root_dense_block_log_likelihood(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *value)
{
  if (!block.root_eligible) {
    return false;
  }
  const int size = static_cast<int>(block.index.size());
  double *covariance = plan.covariance_work.data();
  double *residual = plan.residual_work.data();
  for (int col = 0; col < size; ++col) {
    const int global_col = block.index[static_cast<size_t>(col)];
    for (int row = 0; row < size; ++row) {
      const int global_row = block.index[static_cast<size_t>(row)];
      double base = plan.sampling_covariance[static_cast<size_t>(
        global_row + plan.n * global_col
      )];
      if (row == col) {
        base += extra[global_row];
      }
      if (!std::isfinite(base)) {
        Rf_error("Known-V bridge covariance base must be finite.");
      }
      covariance[static_cast<size_t>(row + size * col)] = base;
    }
    if (!std::isfinite(mean[global_col])) {
      Rf_error("Known-V bridge means must be finite.");
    }
    residual[static_cast<size_t>(col)] =
      plan.y[static_cast<size_t>(global_col)] - mean[global_col];
  }

  double *U = fill_low_rank_factor(plan, block, factors);
  const int rank = block.rank;
  if (rank > 0) {
    const double one = 1.0;
    F77_CALL(dsyrk)(
      "L", "N", &size, &rank, &one, U, &size, &one,
      covariance, &size FCONE FCONE
    );
  }
  *value = cholesky_block_log_likelihood(size, covariance, residual);
  return true;
}

bool finish_low_rank_log_likelihood(
    CovariancePlan &plan,
    int size,
    int rank,
    const double *diagonal,
    double *residual,
    const double *U,
    double log_determinant,
    double base_quadratic,
    double *value)
{
  double *small = plan.small_matrix_work.data();
  double *rhs = plan.rhs_work.data();
  int info = 0;
  F77_CALL(dpotrf)("L", &rank, small, &rank, &info FCONE);
  if (info != 0) {
    Rf_error("Known-V bridge Woodbury system is not positive definite.");
  }
  for (int row = 0; row < rank; ++row) {
    const double root_diagonal = small[
      static_cast<size_t>(row + rank * row)
    ];
    log_determinant += 2.0 * std::log(root_diagonal);
    double transformed = rhs[static_cast<size_t>(row)];
    for (int column = 0; column < row; ++column) {
      transformed -= small[static_cast<size_t>(row + rank * column)] *
        rhs[static_cast<size_t>(column)];
    }
    rhs[static_cast<size_t>(row)] = transformed / root_diagonal;
  }
  double adjustment = 0.0;
  for (int row = 0; row < rank; ++row) {
    adjustment += rhs[static_cast<size_t>(row)] *
      rhs[static_cast<size_t>(row)];
  }
  double quadratic = base_quadratic - adjustment;
  const double cancellation_bound =
    static_cast<double>(size + rank + 1) *
    std::numeric_limits<double>::epsilon() *
    (std::abs(base_quadratic) + std::abs(adjustment));
  if (!(quadratic > cancellation_bound)) {
    const int one = 1;
    F77_CALL(dtrsv)(
      "L", "T", "N", &rank, small, &rank, rhs, &one FCONE FCONE FCONE
    );

    const double minus_one = -1.0;
    const double plus_one = 1.0;
    F77_CALL(dgemv)(
      "N", &size, &rank, &minus_one, U, &size, rhs, &one,
      &plus_one, residual, &one FCONE
    );
    quadratic = F77_CALL(ddot)(&rank, rhs, &one, rhs, &one);
    for (int observation = 0; observation < size; ++observation) {
      const double transformed = residual[static_cast<size_t>(observation)];
      quadratic += transformed * transformed /
        diagonal[static_cast<size_t>(observation)];
    }
  }
  if (!std::isfinite(log_determinant) || !std::isfinite(quadratic)) {
    Rf_error("Known-V bridge Woodbury likelihood is not finite.");
  }
  const double log_two_pi = 1.837877066409345483560659472811;
  *value = -0.5 * (
    static_cast<double>(size) * log_two_pi +
    log_determinant +
    quadratic
  );
  return true;
}

bool diagonal_low_rank_log_likelihood(
    CovariancePlan &plan,
    int size,
    int rank,
    const double *diagonal,
    double *residual,
    const double *U,
    double log_determinant,
    double base_quadratic,
    double *value)
{
  if (rank == 0) {
    const double log_two_pi = 1.837877066409345483560659472811;
    *value = -0.5 * (
      static_cast<double>(size) * log_two_pi +
      log_determinant +
      base_quadratic
    );
    return true;
  }

  double *small = plan.small_matrix_work.data();
  double *rhs = plan.rhs_work.data();
  std::fill(small, small + static_cast<size_t>(rank * rank), 0.0);
  std::fill(rhs, rhs + static_cast<size_t>(rank), 0.0);
  for (int column = 0; column < rank; ++column) {
    for (int row = 0; row < rank; ++row) {
      double entry = row == column ? 1.0 : 0.0;
      for (int observation = 0; observation < size; ++observation) {
        entry += U[static_cast<size_t>(observation + size * row)] *
          U[static_cast<size_t>(observation + size * column)] /
          diagonal[static_cast<size_t>(observation)];
      }
      small[static_cast<size_t>(row + rank * column)] = entry;
    }
    for (int observation = 0; observation < size; ++observation) {
      rhs[static_cast<size_t>(column)] +=
        U[static_cast<size_t>(observation + size * column)] *
        residual[static_cast<size_t>(observation)] /
        diagonal[static_cast<size_t>(observation)];
    }
  }

  return finish_low_rank_log_likelihood(
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

bool sparse_diagonal_low_rank_log_likelihood(
    CovariancePlan &plan,
    const BlockPlan &block,
    int size,
    int rank,
    const double *diagonal,
    double *residual,
    const double *U,
    double log_determinant,
    double base_quadratic,
    double *value)
{
  double *small = plan.small_matrix_work.data();
  double *rhs = plan.rhs_work.data();
  std::fill(small, small + static_cast<size_t>(rank * rank), 0.0);
  std::fill(rhs, rhs + static_cast<size_t>(rank), 0.0);
  for (int column = 0; column < rank; ++column) {
    small[static_cast<size_t>(column + rank * column)] = 1.0;
  }
  for (int observation = 0; observation < size; ++observation) {
    const double inverse_variance = 1.0 /
      diagonal[static_cast<size_t>(observation)];
    const double weighted_residual =
      residual[static_cast<size_t>(observation)] * inverse_variance;
    const std::vector<int> &columns = block.active_columns[
      static_cast<size_t>(observation)
    ];
    for (size_t column_i = 0; column_i < columns.size(); ++column_i) {
      const int column = columns[column_i];
      const double column_value = U[static_cast<size_t>(
        observation + size * column
      )];
      rhs[static_cast<size_t>(column)] +=
        column_value * weighted_residual;
      for (size_t row_i = column_i; row_i < columns.size(); ++row_i) {
        const int row = columns[row_i];
        small[static_cast<size_t>(row + rank * column)] +=
          U[static_cast<size_t>(observation + size * row)] *
          column_value * inverse_variance;
      }
    }
  }

  return finish_low_rank_log_likelihood(
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

void assemble_dense_latent_system(
    CovariancePlan &plan,
    const std::vector<std::vector<int>> &active_columns,
    int size,
    int rank,
    const double *diagonal,
    const double *residual,
    const double *U)
{
  double *small = plan.small_matrix_work.data();
  double *rhs = plan.rhs_work.data();
  std::fill(small, small + static_cast<size_t>(rank * rank), 0.0);
  std::fill(rhs, rhs + static_cast<size_t>(rank), 0.0);
  for (int column = 0; column < rank; ++column) {
    small[static_cast<size_t>(column + rank * column)] = 1.0;
  }
  for (int observation = 0; observation < size; ++observation) {
    const double inverse_variance = 1.0 /
      diagonal[static_cast<size_t>(observation)];
    const double weighted_residual =
      residual[static_cast<size_t>(observation)] * inverse_variance;
    const std::vector<int> &columns = active_columns[
      static_cast<size_t>(observation)
    ];
    for (size_t column_i = 0; column_i < columns.size(); ++column_i) {
      const int column = columns[column_i];
      const double column_value = U[static_cast<size_t>(
        observation + size * column
      )];
      rhs[static_cast<size_t>(column)] +=
        column_value * weighted_residual;
      for (size_t row_i = column_i; row_i < columns.size(); ++row_i) {
        const int row = columns[row_i];
        small[static_cast<size_t>(row + rank * column)] +=
          U[static_cast<size_t>(observation + size * row)] *
          column_value * inverse_variance;
      }
    }
  }
}

bool factor_sparse_latent_system(
    CovariancePlan &plan,
    BlockPlan &block,
    int size,
    int rank,
    const double *diagonal,
    const double *residual,
    const double *U)
{
  std::fill(
    block.sparse_values.begin(),
    block.sparse_values.end(),
    0.0
  );
  for (int column = 0; column < rank; ++column) {
    block.sparse_values[static_cast<size_t>(
      block.sparse_diagonal_positions[static_cast<size_t>(column)]
    )] = 1.0;
  }

  double *rhs = plan.rhs_work.data();
  std::fill(rhs, rhs + static_cast<size_t>(rank), 0.0);
  for (int observation = 0; observation < size; ++observation) {
    const double inverse_variance = 1.0 /
      diagonal[static_cast<size_t>(observation)];
    const double weighted_residual =
      residual[static_cast<size_t>(observation)] * inverse_variance;
    const std::vector<int> &columns = block.sparse_factor_active_columns[
      static_cast<size_t>(observation)
    ];
    const std::vector<int> &positions = block.active_pair_positions[
      static_cast<size_t>(observation)
    ];
    size_t position_i = 0;
    for (size_t column_i = 0; column_i < columns.size(); ++column_i) {
      const int column = columns[column_i];
      const double column_value = U[static_cast<size_t>(
        observation + size * column
      )];
      rhs[static_cast<size_t>(column)] +=
        column_value * weighted_residual;
      for (size_t row_i = column_i; row_i < columns.size(); ++row_i) {
        const int row = columns[row_i];
        block.sparse_values[static_cast<size_t>(positions[position_i])] +=
          U[static_cast<size_t>(observation + size * row)] *
          column_value * inverse_variance;
        ++position_i;
      }
    }
  }

  cholmod_sparse matrix = sparse_matrix_view(block);
  if (!M_cholmod_factorize(
      &matrix,
      block.sparse_factor,
      &plan.cholmod_common_state
    ) || block.sparse_factor->minor < static_cast<size_t>(rank)) {
    return false;
  }
  return true;
}

bool sparse_factor_log_likelihood(
    CovariancePlan &plan,
    BlockPlan &block,
    int size,
    int rank,
    const double *diagonal,
    double *residual,
    const double *U,
    double log_determinant,
    double base_quadratic,
    double *value)
{
  if (!block.sparse_factor_eligible || !factor_sparse_latent_system(
      plan,
      block,
      size,
      rank,
      diagonal,
      residual,
      U
    )) {
    return false;
  }
  log_determinant += M_cholmod_factor_ldetA(block.sparse_factor);

  double *rhs = plan.rhs_work.data();
  cholmod_dense rhs_dense = {};
  rhs_dense.nrow = static_cast<size_t>(rank);
  rhs_dense.ncol = 1;
  rhs_dense.nzmax = static_cast<size_t>(rank);
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
  const double *solution_values = static_cast<const double *>(solution->x);
  double adjustment = 0.0;
  for (int row = 0; row < rank; ++row) {
    adjustment += rhs[static_cast<size_t>(row)] *
      solution_values[static_cast<size_t>(row)];
  }
  std::copy(
    solution_values,
    solution_values + static_cast<size_t>(rank),
    rhs
  );
  M_cholmod_free_dense(&solution, &plan.cholmod_common_state);

  double quadratic = base_quadratic - adjustment;
  const double cancellation_bound =
    static_cast<double>(size + rank + 1) *
    std::numeric_limits<double>::epsilon() *
    (std::abs(base_quadratic) + std::abs(adjustment));
  if (!(quadratic > cancellation_bound)) {
    const int one = 1;
    const double minus_one = -1.0;
    const double plus_one = 1.0;
    F77_CALL(dgemv)(
      "N", &size, &rank, &minus_one, U, &size, rhs, &one,
      &plus_one, residual, &one FCONE
    );
    quadratic = F77_CALL(ddot)(&rank, rhs, &one, rhs, &one);
    for (int observation = 0; observation < size; ++observation) {
      const double transformed = residual[static_cast<size_t>(observation)];
      quadratic += transformed * transformed /
        diagonal[static_cast<size_t>(observation)];
    }
  }
  if (!std::isfinite(log_determinant) || !std::isfinite(quadratic)) {
    Rf_error("Known-V sparse Woodbury likelihood is not finite.");
  }
  const double log_two_pi = 1.837877066409345483560659472811;
  *value = -0.5 * (
    static_cast<double>(size) * log_two_pi +
    log_determinant +
    quadratic
  );
  return true;
}

bool low_rank_block_log_likelihood(
    CovariancePlan &plan,
    BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *value)
{
  if (!block.low_rank_eligible &&
      !(block.sparse_factor_eligible && !block.sparse_factor_block_base)) {
    return false;
  }
  const int size = static_cast<int>(block.index.size());
  const int rank = block.rank;
  double *diagonal = plan.diagonal_work.data();
  double *residual = plan.residual_work.data();
  double log_determinant = 0.0;
  double base_quadratic = 0.0;
  for (int row = 0; row < size; ++row) {
    const int global = block.index[static_cast<size_t>(row)];
    const double variance = plan.sampling_covariance[
      static_cast<size_t>(global + plan.n * global)
    ] + extra[global];
    if (!(variance > 0.0) || !std::isfinite(variance) ||
        !std::isfinite(mean[global])) {
      return false;
    }
    diagonal[static_cast<size_t>(row)] = variance;
    residual[static_cast<size_t>(row)] =
      plan.y[static_cast<size_t>(global)] - mean[global];
    log_determinant += std::log(variance);
    base_quadratic += residual[static_cast<size_t>(row)] *
      residual[static_cast<size_t>(row)] / variance;
  }

  double *U = fill_low_rank_factor(plan, block, factors);
  if (!block.sparse_factor_block_base && sparse_factor_log_likelihood(
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
  if (!block.low_rank_eligible) {
    return false;
  }
  if (block.sparse_assembly_eligible) {
    return sparse_diagonal_low_rank_log_likelihood(
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
    );
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
  double *covariance = plan.covariance_work.data();
  double *base_rhs = plan.base_rhs_work.data();
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
  const double one = 1.0;
  const int rhs_columns = rank + 1;
  for (const SamplingBlock &sampling_block : block.sampling_blocks) {
    const int block_size = static_cast<int>(sampling_block.local_rows.size());
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
      log_determinant += sampling_block.log_determinant;
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
      base_rhs[static_cast<size_t>(row)] = residual[
        static_cast<size_t>(local)
      ];
      for (int column = 0; column < rank; ++column) {
        base_rhs[static_cast<size_t>(
          row + block_size * (column + 1)
        )] = U[static_cast<size_t>(local + size * column)];
      }
    }

    if (!fixed_base) {
      int info = 0;
      F77_CALL(dpotrf)(
        "L", &block_size, covariance, &block_size, &info FCONE
      );
      if (info != 0) {
        return false;
      }
      for (int row = 0; row < block_size; ++row) {
        log_determinant += 2.0 * std::log(covariance[static_cast<size_t>(
          row + block_size * row
        )]);
      }
    }
    F77_CALL(dtrsm)(
      "L", "L", "N", "N", &block_size, &rhs_columns, &one,
      covariance, &block_size, base_rhs, &block_size
      FCONE FCONE FCONE FCONE
    );
    for (int row = 0; row < block_size; ++row) {
      const int local = sampling_block.local_rows[static_cast<size_t>(row)];
      const double transformed_residual = base_rhs[static_cast<size_t>(row)];
      residual[static_cast<size_t>(local)] = transformed_residual;
      base_quadratic += transformed_residual * transformed_residual;
      for (int column = 0; column < rank; ++column) {
        U[static_cast<size_t>(local + size * column)] = base_rhs[
          static_cast<size_t>(row + block_size * (column + 1))
        ];
      }
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
  double *covariance = plan.covariance_work.data();
  double *base_rhs = plan.base_rhs_work.data();
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

  const double one = 1.0;
  const int rhs_columns = rank + 1;
  for (const SamplingBlock &sampling_block : block.sampling_blocks) {
    const int block_size = static_cast<int>(sampling_block.local_rows.size());
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
      base_rhs[static_cast<size_t>(row)] = residual[
        static_cast<size_t>(local)
      ];
      for (int column = 0; column < rank; ++column) {
        base_rhs[static_cast<size_t>(
          row + block_size * (column + 1)
        )] = U[static_cast<size_t>(local + size * column)];
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
    }
    F77_CALL(dtrsm)(
      "L", "L", "N", "N", &block_size, &rhs_columns, &one,
      covariance, &block_size, base_rhs, &block_size
      FCONE FCONE FCONE FCONE
    );
    // The first solve gives L^-1[r,U]; the second gives B^-1[r,U].
    for (int row = 0; row < block_size; ++row) {
      const int local = sampling_block.local_rows[static_cast<size_t>(row)];
      residual[static_cast<size_t>(local)] = base_rhs[static_cast<size_t>(row)];
      for (int column = 0; column < rank; ++column) {
        U[static_cast<size_t>(local + size * column)] = base_rhs[
          static_cast<size_t>(row + block_size * (column + 1))
        ];
      }
    }
    F77_CALL(dtrsm)(
      "L", "L", "T", "N", &block_size, &rhs_columns, &one,
      covariance, &block_size, base_rhs, &block_size
      FCONE FCONE FCONE FCONE
    );
    for (int row = 0; row < block_size; ++row) {
      const int local = sampling_block.local_rows[static_cast<size_t>(row)];
      for (int column = 0; column <= rank; ++column) {
        base_precision[static_cast<size_t>(
          local + size * column
        )] = base_rhs[static_cast<size_t>(
          row + block_size * column
        )];
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

} // namespace

extern "C" SEXP RoBMA_known_v_covariance_plan_create(
    SEXP y,
    SEXP sampling_covariance,
    SEXP random_covariance_factors,
    SEXP block_indices)
{
  CovariancePlan *plan = make_plan(
    y,
    sampling_covariance,
    random_covariance_factors,
    block_indices
  );
  SEXP pointer = PROTECT(R_MakeExternalPtr(plan, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(pointer, finalize_plan, TRUE);
  int low_rank_blocks = 0;
  int markov_blocks = 0;
  int sparse_assembly_blocks = 0;
  int sparse_factor_blocks = 0;
  int block_base_blocks = 0;
  int spectral_blocks = 0;
  int root_dense_blocks = 0;
  int dense_blocks = 0;
  for (const BlockPlan &block : plan->blocks) {
    if (block.markov_eligible) {
      ++markov_blocks;
    } else if (block.sparse_factor_eligible && block.sparse_factor_block_base) {
      ++block_base_blocks;
      ++sparse_factor_blocks;
    } else if (block.sparse_factor_eligible) {
      ++low_rank_blocks;
      sparse_assembly_blocks += block.sparse_assembly_eligible ? 1 : 0;
      ++sparse_factor_blocks;
    } else if (block.low_rank_eligible) {
      ++low_rank_blocks;
      sparse_assembly_blocks += block.sparse_assembly_eligible ? 1 : 0;
    } else if (block.block_base_eligible) {
      ++block_base_blocks;
    } else if (block.spectral_eligible) {
      ++spectral_blocks;
    } else {
      ++dense_blocks;
      root_dense_blocks += block.root_eligible ? 1 : 0;
    }
  }
  SEXP low_rank_blocks_sexp = PROTECT(Rf_ScalarInteger(low_rank_blocks));
  SEXP markov_blocks_sexp = PROTECT(Rf_ScalarInteger(markov_blocks));
  SEXP sparse_assembly_blocks_sexp = PROTECT(Rf_ScalarInteger(
    sparse_assembly_blocks
  ));
  SEXP sparse_factor_blocks_sexp = PROTECT(Rf_ScalarInteger(
    sparse_factor_blocks
  ));
  SEXP block_base_blocks_sexp = PROTECT(Rf_ScalarInteger(block_base_blocks));
  SEXP spectral_blocks_sexp = PROTECT(Rf_ScalarInteger(spectral_blocks));
  SEXP root_dense_blocks_sexp = PROTECT(Rf_ScalarInteger(root_dense_blocks));
  SEXP dense_blocks_sexp = PROTECT(Rf_ScalarInteger(dense_blocks));
  Rf_setAttrib(
    pointer,
    Rf_install("markov_blocks"),
    markov_blocks_sexp
  );
  Rf_setAttrib(
    pointer,
    Rf_install("root_dense_blocks"),
    root_dense_blocks_sexp
  );
  Rf_setAttrib(
    pointer,
    Rf_install("low_rank_blocks"),
    low_rank_blocks_sexp
  );
  Rf_setAttrib(
    pointer,
    Rf_install("sparse_assembly_blocks"),
    sparse_assembly_blocks_sexp
  );
  Rf_setAttrib(
    pointer,
    Rf_install("sparse_factor_blocks"),
    sparse_factor_blocks_sexp
  );
  Rf_setAttrib(
    pointer,
    Rf_install("block_base_blocks"),
    block_base_blocks_sexp
  );
  Rf_setAttrib(
    pointer,
    Rf_install("spectral_blocks"),
    spectral_blocks_sexp
  );
  Rf_setAttrib(
    pointer,
    Rf_install("dense_blocks"),
    dense_blocks_sexp
  );
  UNPROTECT(9);
  return pointer;
}

extern "C" SEXP RoBMA_known_v_covariance_plan_loglik(
    SEXP pointer,
    SEXP mean,
    SEXP random_covariance_states,
    SEXP extra_variance)
{
  CovariancePlan *plan = plan_pointer(pointer);
  std::vector<CovarianceFactor> states = covariance_states(
    random_covariance_states,
    *plan
  );
  return Rf_ScalarReal(plan_log_likelihood(
    *plan,
    mean,
    states,
    extra_variance
  ));
}

extern "C" SEXP RoBMA_known_v_covariance_plan_conditional_loglik(
    SEXP pointer,
    SEXP mean,
    SEXP random_covariance_states,
    SEXP extra_variance)
{
  CovariancePlan *plan = plan_pointer(pointer);
  std::vector<CovarianceFactor> states = covariance_states(
    random_covariance_states,
    *plan
  );
  return plan_conditional_log_likelihood(
    *plan,
    mean,
    states,
    extra_variance
  );
}

extern "C" SEXP RoBMA_known_v_block_mvn_loglik(
    SEXP y,
    SEXP mean,
    SEXP sampling_covariance,
    SEXP random_covariance_factors,
    SEXP block_indices,
    SEXP extra_variance)
{
  CovariancePlan *plan = make_plan(
    y,
    sampling_covariance,
    random_covariance_factors,
    block_indices
  );
  std::vector<CovarianceFactor> factors = covariance_factors(
    random_covariance_factors,
    plan->n
  );
  const double value = plan_log_likelihood(
    *plan,
    mean,
    factors,
    extra_variance
  );
  delete plan;
  return Rf_ScalarReal(value);
}
