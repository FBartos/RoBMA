#include <R_ext/Error.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
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

// Private implementation fragments remain in one translation unit so the
// exact solver dispatch and compiler-visible numerical kernels are unchanged.
#include "r-known-v-common.cc.inc"
#include "r-known-v-plan.cc.inc"
#include "r-known-v-low-rank.cc.inc"
#include "r-known-v-structured.cc.inc"
#include "r-known-v-evaluate.cc.inc"

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
  int fixed_known_group_blocks = 0;
  int sparse_assembly_blocks = 0;
  int sparse_factor_blocks = 0;
  int block_base_blocks = 0;
  int spectral_blocks = 0;
  int root_dense_blocks = 0;
  int dense_blocks = 0;
  for (const BlockPlan &block : plan->blocks) {
    if (block.markov_eligible) {
      ++markov_blocks;
    } else if (block.fixed_known_group_eligible) {
      ++fixed_known_group_blocks;
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
  SEXP fixed_known_group_blocks_sexp = PROTECT(Rf_ScalarInteger(
    fixed_known_group_blocks
  ));
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
    Rf_install("fixed_known_group_blocks"),
    fixed_known_group_blocks_sexp
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
  UNPROTECT(10);
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

extern "C" SEXP RoBMA_known_v_covariance_plan_loglik_batch(
    SEXP pointer,
    SEXP means,
    SEXP random_covariance_states,
    SEXP extra_variances)
{
  CovariancePlan *plan = plan_pointer(pointer);
  const int draws = require_batched_plan_inputs(
    *plan,
    means,
    random_covariance_states,
    extra_variances
  );
  SEXP out = PROTECT(Rf_allocVector(REALSXP, draws));
  const double *mean_values = REAL(means);
  const double *extra_values = REAL(extra_variances);
  for (int draw = 0; draw < draws; ++draw) {
    std::vector<CovarianceFactor> states = covariance_states(
      VECTOR_ELT(random_covariance_states, draw),
      *plan
    );
    REAL(out)[draw] = plan_log_likelihood_values(
      *plan,
      mean_values + static_cast<size_t>(draw) * plan->n,
      states,
      extra_values + static_cast<size_t>(draw) * plan->n
    );
  }
  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_known_v_covariance_plan_group_iid_variance_grid_loglik(
    SEXP pointer,
    SEXP means,
    SEXP group_variances,
    SEXP diagonal_variances)
{
  CovariancePlan *plan = plan_pointer(pointer);
  if (TYPEOF(means) != REALSXP || !Rf_isMatrix(means) ||
      TYPEOF(group_variances) != REALSXP ||
      !Rf_isMatrix(group_variances) ||
      TYPEOF(diagonal_variances) != REALSXP ||
      !Rf_isMatrix(diagonal_variances)) {
    Rf_error(
      "Known-V group-IID grid inputs must be numeric matrices."
    );
  }
  if (!plan->factors.empty()) {
    Rf_error("Known-V group-IID grid plans cannot contain random factors.");
  }
  SEXP mean_dim = Rf_getAttrib(means, R_DimSymbol);
  SEXP group_dim = Rf_getAttrib(group_variances, R_DimSymbol);
  SEXP diagonal_dim = Rf_getAttrib(diagonal_variances, R_DimSymbol);
  const int draws = INTEGER(mean_dim)[1];
  const int grids = INTEGER(group_dim)[0];
  if (INTEGER(mean_dim)[0] != plan->n ||
      grids < 1 ||
      INTEGER(group_dim)[1] != draws ||
      INTEGER(diagonal_dim)[0] != grids ||
      INTEGER(diagonal_dim)[1] != draws) {
    Rf_error("Known-V group-IID grid inputs have inconsistent dimensions.");
  }
  const double *mean_values = REAL(means);
  const double *group_values = REAL(group_variances);
  const double *diagonal_values = REAL(diagonal_variances);
  for (int index = 0; index < grids * draws; ++index) {
    if (!std::isfinite(group_values[index]) || group_values[index] < 0.0 ||
        !std::isfinite(diagonal_values[index]) ||
        diagonal_values[index] < 0.0) {
      Rf_error(
        "Known-V group-IID variances must be finite and non-negative."
      );
    }
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, grids, draws));
  std::fill(
    REAL(out),
    REAL(out) + static_cast<size_t>(grids) * static_cast<size_t>(draws),
    0.0
  );
  size_t max_block_size = 0;
  for (const BlockPlan &block : plan->blocks) {
    max_block_size = std::max(max_block_size, block.index.size());
  }
  std::vector<double> sampling(max_block_size * max_block_size);
  std::vector<double> eigenvectors(max_block_size * max_block_size);
  std::vector<double> eigenvalues(max_block_size);
  std::vector<double> group_loading(max_block_size);
  std::vector<double> residual(max_block_size);
  std::vector<double> transformed_residual(max_block_size);
  std::vector<double> dense_residual(max_block_size);
  std::vector<double> covariance(max_block_size * max_block_size);
  const int lwork = std::max(1, 3 * static_cast<int>(max_block_size) - 1);
  std::vector<double> eigen_work(static_cast<size_t>(lwork));
  const double log_two_pi = 1.837877066409345483560659472811;

  for (const BlockPlan &block : plan->blocks) {
    const int size = static_cast<int>(block.index.size());
    for (int column = 0; column < size; ++column) {
      const int global_column = block.index[static_cast<size_t>(column)];
      for (int row = 0; row < size; ++row) {
        const int global_row = block.index[static_cast<size_t>(row)];
        sampling[static_cast<size_t>(row + size * column)] =
          plan->sampling_covariance[static_cast<size_t>(
            global_row + plan->n * global_column
          )];
      }
    }
    std::copy(
      sampling.begin(),
      sampling.begin() + static_cast<size_t>(size * size),
      eigenvectors.begin()
    );
    int info = 0;
    F77_CALL(dsyev)(
      "V", "L", &size, eigenvectors.data(), &size,
      eigenvalues.data(), eigen_work.data(), &lwork, &info FCONE FCONE
    );
    if (info != 0) {
      Rf_error("Known-V sampling covariance eigendecomposition failed.");
    }
    for (int eigen = 0; eigen < size; ++eigen) {
      double loading = 0.0;
      for (int row = 0; row < size; ++row) {
        loading += eigenvectors[static_cast<size_t>(row + size * eigen)];
      }
      group_loading[static_cast<size_t>(eigen)] = loading;
    }

    for (int draw = 0; draw < draws; ++draw) {
      const double *draw_mean = mean_values +
        static_cast<size_t>(draw) * static_cast<size_t>(plan->n);
      for (int row = 0; row < size; ++row) {
        const int global = block.index[static_cast<size_t>(row)];
        const double mean = draw_mean[global];
        if (!std::isfinite(mean)) {
          Rf_error("Known-V bridge means must be finite.");
        }
        residual[static_cast<size_t>(row)] =
          plan->y[static_cast<size_t>(global)] - mean;
      }
      for (int eigen = 0; eigen < size; ++eigen) {
        double value = 0.0;
        for (int row = 0; row < size; ++row) {
          value += eigenvectors[static_cast<size_t>(row + size * eigen)] *
            residual[static_cast<size_t>(row)];
        }
        transformed_residual[static_cast<size_t>(eigen)] = value;
      }

      for (int grid = 0; grid < grids; ++grid) {
        const size_t variance_index = static_cast<size_t>(grid) +
          static_cast<size_t>(grids) * static_cast<size_t>(draw);
        const double diagonal_variance = diagonal_values[variance_index];
        const double group_variance = group_values[variance_index];

        double log_determinant = 0.0;
        double quadratic = 0.0;
        double group_precision = 0.0;
        double group_residual = 0.0;
        bool use_dense_reference = false;
        for (int eigen = 0; eigen < size; ++eigen) {
          const double variance =
            eigenvalues[static_cast<size_t>(eigen)] + diagonal_variance;
          if (!(variance > 0.0) || !std::isfinite(variance)) {
            use_dense_reference = true;
            break;
          }
          const double loading = group_loading[static_cast<size_t>(eigen)];
          const double transformed =
            transformed_residual[static_cast<size_t>(eigen)];
          log_determinant += std::log(variance);
          quadratic += transformed * transformed / variance;
          group_precision += loading * loading / variance;
          group_residual += loading * transformed / variance;
        }
        const double update = 1.0 + group_variance * group_precision;
        if (!(update > 0.0) || !std::isfinite(update)) {
          use_dense_reference = true;
        }

        double block_log_likelihood = 0.0;
        if (use_dense_reference) {
          for (int column = 0; column < size; ++column) {
            for (int row = 0; row < size; ++row) {
              covariance[static_cast<size_t>(row + size * column)] =
                sampling[static_cast<size_t>(row + size * column)] +
                group_variance +
                (row == column ? diagonal_variance : 0.0);
            }
          }
          std::copy(
            residual.begin(),
            residual.begin() + static_cast<size_t>(size),
            dense_residual.begin()
          );
          block_log_likelihood = cholesky_block_log_likelihood(
            size, covariance.data(), dense_residual.data()
          );
        } else {
          log_determinant += std::log(update);
          quadratic -= group_variance * group_residual * group_residual /
            update;
          block_log_likelihood = -0.5 * (
            static_cast<double>(size) * log_two_pi +
            log_determinant + quadratic
          );
        }
        REAL(out)[static_cast<size_t>(grid) +
          static_cast<size_t>(grids) * static_cast<size_t>(draw)] +=
          block_log_likelihood;
      }
    }
  }

  UNPROTECT(1);
  return out;
}


extern "C" SEXP RoBMA_known_v_covariance_plan_affine_grid_loglik(
    SEXP pointer,
    SEXP means,
    SEXP base_covariances,
    SEXP update_covariances,
    SEXP reference_coefficient,
    SEXP coefficients)
{
  CovariancePlan *plan = plan_pointer(pointer);
  return plan_affine_grid_log_likelihood(
    *plan,
    means,
    base_covariances,
    update_covariances,
    reference_coefficient,
    coefficients
  );
}


extern "C" SEXP RoBMA_known_v_covariance_plan_factor_grid_loglik(
    SEXP pointer,
    SEXP means,
    SEXP random_covariance_states,
    SEXP extra_variances,
    SEXP update_grid)
{
  CovariancePlan *plan = plan_pointer(pointer);
  return plan_factor_grid_log_likelihood(
    *plan,
    means,
    random_covariance_states,
    extra_variances,
    update_grid
  );
}



extern "C" SEXP RoBMA_known_v_covariance_plan_location_quadratic_batch(
    SEXP pointer,
    SEXP means,
    SEXP bases,
    SEXP random_covariance_states,
    SEXP extra_variances)
{
  CovariancePlan *plan = plan_pointer(pointer);
  const int draws = require_batched_plan_inputs(
    *plan,
    means,
    random_covariance_states,
    extra_variances
  );
  if (TYPEOF(bases) != REALSXP || !Rf_isMatrix(bases)) {
    Rf_error("Known-V location bases must be a numeric matrix.");
  }
  SEXP basis_dim = Rf_getAttrib(bases, R_DimSymbol);
  if (INTEGER(basis_dim)[0] != plan->n ||
      INTEGER(basis_dim)[1] != draws) {
    Rf_error("Known-V location bases have inconsistent dimensions.");
  }

  SEXP output = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP linear = PROTECT(Rf_allocVector(REALSXP, draws));
  SEXP quadratic = PROTECT(Rf_allocVector(REALSXP, draws));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("linear"));
  SET_STRING_ELT(names, 1, Rf_mkChar("quadratic"));
  SET_VECTOR_ELT(output, 0, linear);
  SET_VECTOR_ELT(output, 1, quadratic);
  Rf_setAttrib(output, R_NamesSymbol, names);

  const double *mean_values = REAL(means);
  const double *basis_values = REAL(bases);
  const double *extra_values = REAL(extra_variances);
  for (int draw = 0; draw < draws; ++draw) {
    std::vector<CovarianceFactor> states = covariance_states(
      VECTOR_ELT(random_covariance_states, draw),
      *plan
    );
    plan_location_quadratic_values(
      *plan,
      mean_values + static_cast<size_t>(draw) * plan->n,
      basis_values + static_cast<size_t>(draw) * plan->n,
      states,
      extra_values + static_cast<size_t>(draw) * plan->n,
      REAL(linear) + draw,
      REAL(quadratic) + draw
    );
  }
  UNPROTECT(4);
  return output;
}

extern "C" SEXP RoBMA_known_v_covariance_plan_conditional_loglik_batch(
    SEXP pointer,
    SEXP means,
    SEXP random_covariance_states,
    SEXP extra_variances)
{
  CovariancePlan *plan = plan_pointer(pointer);
  const int draws = require_batched_plan_inputs(
    *plan,
    means,
    random_covariance_states,
    extra_variances
  );
  SEXP output = PROTECT(Rf_allocMatrix(REALSXP, plan->n, draws));
  const double *mean_values = REAL(means);
  const double *extra_values = REAL(extra_variances);
  for (int draw = 0; draw < draws; ++draw) {
    std::vector<CovarianceFactor> states = covariance_states(
      VECTOR_ELT(random_covariance_states, draw),
      *plan
    );
    const ConditionalOutput destination = {
      REAL(output) + static_cast<size_t>(draw) * plan->n,
      nullptr,
      nullptr
    };
    plan_conditional_values(
      *plan,
      mean_values + static_cast<size_t>(draw) * plan->n,
      states,
      extra_values + static_cast<size_t>(draw) * plan->n,
      destination
    );
  }
  UNPROTECT(1);
  return output;
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

extern "C" SEXP RoBMA_known_v_covariance_plan_conditional_summary_batch(
    SEXP pointer,
    SEXP means,
    SEXP random_covariance_states,
    SEXP extra_variances)
{
  CovariancePlan *plan = plan_pointer(pointer);
  const int draws = require_batched_plan_inputs(
    *plan,
    means,
    random_covariance_states,
    extra_variances
  );
  SEXP output = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP residual = PROTECT(Rf_allocMatrix(REALSXP, plan->n, draws));
  SEXP variance = PROTECT(Rf_allocMatrix(REALSXP, plan->n, draws));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("residual"));
  SET_STRING_ELT(names, 1, Rf_mkChar("variance"));
  SET_VECTOR_ELT(output, 0, residual);
  SET_VECTOR_ELT(output, 1, variance);
  Rf_setAttrib(output, R_NamesSymbol, names);

  const double *mean_values = REAL(means);
  const double *extra_values = REAL(extra_variances);
  for (int draw = 0; draw < draws; ++draw) {
    std::vector<CovarianceFactor> states = covariance_states(
      VECTOR_ELT(random_covariance_states, draw),
      *plan
    );
    const ConditionalOutput destination = {
      nullptr,
      REAL(residual) + static_cast<size_t>(draw) * plan->n,
      REAL(variance) + static_cast<size_t>(draw) * plan->n
    };
    plan_conditional_values(
      *plan,
      mean_values + static_cast<size_t>(draw) * plan->n,
      states,
      extra_values + static_cast<size_t>(draw) * plan->n,
      destination
    );
  }
  UNPROTECT(4);
  return output;
}
