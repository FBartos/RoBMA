#include <R_ext/Error.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <Rinternals.h>

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

struct CovarianceFactor {
  FactorType type;
  const double *covariance;
  const double *model_matrix;
  const int *group_map;
  const double *coefficient_factor;
  const double *group_covariance;
  const double *row_scale;
  int n_columns;
  int n_groups;
};

struct PlanFactor {
  FactorType type;
  std::vector<double> model_matrix;
  std::vector<int> group_map;
  std::vector<double> group_covariance;
  std::vector<double> group_factor;
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
};

struct BlockPlan {
  std::vector<int> index;
  std::vector<LowRankGroup> low_rank_groups;
  std::vector<SamplingBlock> sampling_blocks;
  std::vector<double> sampling_eigenvalues;
  std::vector<double> sampling_eigenvectors;
  int rank;
  bool low_rank_eligible;
  bool block_base_eligible;
  bool spectral_eligible;
  bool root_eligible;
};

struct CovariancePlan {
  int n;
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
};

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
  out.root_eligible = compatible;
  out.low_rank_eligible = sampling_diagonal && compatible && rank < size;
  if (!sampling_diagonal && compatible && rank < size) {
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
        SamplingBlock sampling_block;
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
        out.sampling_blocks.push_back(sampling_block);
      }
    }
  }
  out.block_base_eligible = out.sampling_blocks.size() > 1;
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
    plan->blocks.push_back(block);
  }
  if (std::find(membership.begin(), membership.end(), -1) != membership.end()) {
    delete plan;
    Rf_error("Covariance block indices must partition every observation.");
  }

  plan->covariance_work.resize(max_block_size * max_block_size);
  plan->residual_work.resize(max_block_size);
  plan->diagonal_work.resize(max_block_size);
  plan->low_rank_work.resize(max_block_size * max_rank);
  plan->small_matrix_work.resize(max_rank * max_rank);
  plan->rhs_work.resize(max_rank);
  plan->base_rhs_work.resize(max_block_size * (max_rank + 1));
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
        double transformed = 0.0;
        for (int coefficient = 0; coefficient < q; ++coefficient) {
          transformed += factor.model_matrix[
            global + plan.n * coefficient
          ] * factor.coefficient_factor[
            coefficient + q * root_column
          ];
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

bool low_rank_block_log_likelihood(
    CovariancePlan &plan,
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *value)
{
  if (!block.low_rank_eligible) {
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
    const BlockPlan &block,
    const std::vector<CovarianceFactor> &factors,
    const double *mean,
    const double *extra,
    double *value)
{
  if (!block.block_base_eligible) {
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
    std::copy(
      sampling_block.covariance.begin(),
      sampling_block.covariance.end(),
      covariance
    );
    for (int row = 0; row < block_size; ++row) {
      const int local = sampling_block.local_rows[static_cast<size_t>(row)];
      const int global = block.index[static_cast<size_t>(local)];
      covariance[static_cast<size_t>(row + block_size * row)] += extra[global];
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
  for (const BlockPlan &block : plan.blocks) {
    double block_value = 0.0;
    if (!low_rank_block_log_likelihood(
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
  for (const BlockPlan &block : plan.blocks) {
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
  int block_base_blocks = 0;
  int spectral_blocks = 0;
  int root_dense_blocks = 0;
  int dense_blocks = 0;
  for (const BlockPlan &block : plan->blocks) {
    if (block.low_rank_eligible) {
      ++low_rank_blocks;
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
  SEXP block_base_blocks_sexp = PROTECT(Rf_ScalarInteger(block_base_blocks));
  SEXP spectral_blocks_sexp = PROTECT(Rf_ScalarInteger(spectral_blocks));
  SEXP root_dense_blocks_sexp = PROTECT(Rf_ScalarInteger(root_dense_blocks));
  SEXP dense_blocks_sexp = PROTECT(Rf_ScalarInteger(dense_blocks));
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
  UNPROTECT(6);
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
