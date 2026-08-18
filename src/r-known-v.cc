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
