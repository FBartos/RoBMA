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
