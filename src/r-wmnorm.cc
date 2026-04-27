#include <Rinternals.h>
#include <R_ext/Error.h>

#include <cmath>
#include <vector>

#include "source/wmnorm.h"

static int matrix_nrow(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[0];
}

static int matrix_ncol(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[1];
}

static void check_real(SEXP x, const char *name)
{
  if (TYPEOF(x) != REALSXP) {
    Rf_error("'%s' must be numeric.", name);
  }
}

static void check_integer(SEXP x, const char *name)
{
  if (TYPEOF(x) != INTSXP) {
    Rf_error("'%s' must be integer.", name);
  }
}

static void validate_mix_mapping(const double *mapping, const int *mapping_max,
                                 int n_branches, int n_map, int n_crit,
                                 int n_omega)
{
  if (n_omega < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  for (int branch = 0; branch < n_branches; ++branch) {
    const int active_cuts = mapping_max[branch];
    if (active_cuts == NA_INTEGER || active_cuts < 0 ||
        active_cuts > n_map || active_cuts > n_crit) {
      Rf_error("'crit_yi_mapping_max' contains an invalid active cutoff count.");
    }

    for (int j = 0; j < active_cuts; ++j) {
      const double raw_index = mapping[j + n_map * branch];

      if (!std::isfinite(raw_index)) {
        Rf_error("'crit_yi_mapping' contains an invalid cutoff index.");
      }

      const int index = static_cast<int>(raw_index);
      if (raw_index != static_cast<double>(index) ||
          index < 1 || index > n_crit) {
        Rf_error("'crit_yi_mapping' contains an invalid cutoff index.");
      }
      if (index >= n_omega) {
        Rf_error("'omega' must have one more column than the largest mapped cutoff index.");
      }
    }
  }
}

static void validate_mix_inputs(SEXP value, SEXP mean, SEXP sd, SEXP omega,
                                SEXP crit_yi, SEXP bias_indicator,
                                SEXP crit_yi_mapping, SEXP crit_yi_mapping_max,
                                int *S, int *K, int *n_crit, int *n_map,
                                int *n_branches, int *n_omega)
{
  check_real(value, "value");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(omega, "omega");
  check_real(crit_yi, "crit_yi");
  check_integer(bias_indicator, "bias_indicator");
  check_real(crit_yi_mapping, "crit_yi_mapping");
  check_integer(crit_yi_mapping_max, "crit_yi_mapping_max");

  *S      = matrix_nrow(mean, "mean");
  *K      = matrix_ncol(mean, "mean");
  *n_crit = matrix_nrow(crit_yi, "crit_yi");
  *n_map  = matrix_nrow(crit_yi_mapping, "crit_yi_mapping");
  *n_branches = Rf_length(crit_yi_mapping_max);
  *n_omega    = matrix_ncol(omega, "omega");

  if (matrix_nrow(sd, "sd") != *S || matrix_ncol(sd, "sd") != *K) {
    Rf_error("'sd' dimensions must match 'mean'.");
  }
  if (matrix_nrow(omega, "omega") != *S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (matrix_ncol(crit_yi, "crit_yi") != *K) {
    Rf_error("'crit_yi' must have one column per observation.");
  }
  if (matrix_ncol(crit_yi_mapping, "crit_yi_mapping") != *n_branches) {
    Rf_error("'crit_yi_mapping' columns must match 'crit_yi_mapping_max'.");
  }
  if (Rf_length(value) != *K) {
    Rf_error("'value' must have one value per observation.");
  }
  if (Rf_length(bias_indicator) != *S) {
    Rf_error("'bias_indicator' must have one value per posterior sample.");
  }

  validate_mix_mapping(
    REAL(crit_yi_mapping), INTEGER(crit_yi_mapping_max),
    *n_branches, *n_map, *n_crit, *n_omega
  );
}

static int branch_mapping_max(int branch, const int *mapping_max, int n_branches, int n_map, int n_crit)
{
  if (branch == NA_INTEGER || branch < 1 || branch > n_branches) {
    Rf_error("'bias_indicator' contains an invalid branch index.");
  }

  const int out = mapping_max[branch - 1];
  if (out == NA_INTEGER || out < 0 || out > n_map || out > n_crit) {
    Rf_error("'crit_yi_mapping_max' contains an invalid active cutoff count.");
  }
  return out;
}

static void prepare_mix_row_metadata(const int *bias_indicator,
                                     const int *mapping_max,
                                     const double *mapping,
                                     int S, int n_branches, int n_map,
                                     int n_crit,
                                     std::vector<int> &active_cuts,
                                     std::vector<const double *> &mapping_ptr)
{
  active_cuts.resize(S);
  mapping_ptr.resize(S);

  for (int s = 0; s < S; ++s) {
    const int branch = bias_indicator[s];
    active_cuts[s] = branch_mapping_max(
      branch, mapping_max, n_branches, n_map, n_crit
    );
    mapping_ptr[s] = mapping + n_map * (branch - 1);
  }
}

extern "C" SEXP RoBMA_wnorm_mix_logpdf_matrix(SEXP yi, SEXP mean, SEXP sd,
                                               SEXP omega, SEXP crit_yi,
                                               SEXP bias_indicator,
                                               SEXP crit_yi_mapping,
                                               SEXP crit_yi_mapping_max)
{
  int S, K, n_crit, n_map, n_branches, n_omega;
  validate_mix_inputs(
    yi, mean, sd, omega, crit_yi, bias_indicator, crit_yi_mapping,
    crit_yi_mapping_max, &S, &K, &n_crit, &n_map, &n_branches, &n_omega
  );
  (void)n_omega;

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const double *yi_p       = REAL(yi);
  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *omega_p    = REAL(omega);
  const double *crit_p     = REAL(crit_yi);
  const int *bias_p        = INTEGER(bias_indicator);
  const double *mapping_p  = REAL(crit_yi_mapping);
  const int *mapping_max_p = INTEGER(crit_yi_mapping_max);
  double *out_p            = REAL(out);
  std::vector<int> active_cuts(S);
  std::vector<const double *> mapping_ptr(S);
  prepare_mix_row_metadata(
    bias_p, mapping_max_p, mapping_p, S, n_branches, n_map, n_crit,
    active_cuts, mapping_ptr
  );

  for (int k = 0; k < K; ++k) {
    const double *crit_k = crit_p + n_crit * k;

    for (int s = 0; s < S; ++s) {
      out_p[s + S * k] = cpp_wnorm_mix_lpdf(
        yi_p[k],
        mean_p[s + S * k],
        sd_p[s + S * k],
        crit_k,
        omega_p + s,
        mapping_ptr[s],
        active_cuts[s],
        S
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_wnorm_mix_cdf_matrix(SEXP q, SEXP mean, SEXP sd,
                                            SEXP omega, SEXP crit_yi,
                                            SEXP bias_indicator,
                                            SEXP crit_yi_mapping,
                                            SEXP crit_yi_mapping_max,
                                            SEXP lower_tail, SEXP log_p)
{
  int S, K, n_crit, n_map, n_branches, n_omega;
  validate_mix_inputs(
    q, mean, sd, omega, crit_yi, bias_indicator, crit_yi_mapping,
    crit_yi_mapping_max, &S, &K, &n_crit, &n_map, &n_branches, &n_omega
  );
  (void)n_omega;

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const bool lower = Rf_asLogical(lower_tail) == TRUE;
  const bool log   = Rf_asLogical(log_p) == TRUE;

  const double *q_p         = REAL(q);
  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *omega_p     = REAL(omega);
  const double *crit_p      = REAL(crit_yi);
  const int *bias_p         = INTEGER(bias_indicator);
  const double *mapping_p   = REAL(crit_yi_mapping);
  const int *mapping_max_p  = INTEGER(crit_yi_mapping_max);
  double *out_p             = REAL(out);
  std::vector<int> active_cuts(S);
  std::vector<const double *> mapping_ptr(S);
  prepare_mix_row_metadata(
    bias_p, mapping_max_p, mapping_p, S, n_branches, n_map, n_crit,
    active_cuts, mapping_ptr
  );

  for (int k = 0; k < K; ++k) {
    const double *crit_k = crit_p + n_crit * k;

    for (int s = 0; s < S; ++s) {
      out_p[s + S * k] = cpp_wnorm_mix_cdf(
        q_p[k],
        mean_p[s + S * k],
        sd_p[s + S * k],
        crit_k,
        omega_p + s,
        mapping_ptr[s],
        active_cuts[s],
        lower,
        log,
        S
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_wnorm_mix_logpdf_precomp_matrix(SEXP yi, SEXP mean,
                                                       SEXP sd, SEXP omega,
                                                       SEXP crit_yi,
                                                       SEXP bias_indicator,
                                                       SEXP crit_yi_mapping,
                                                       SEXP crit_yi_mapping_max,
                                                       SEXP log_norm)
{
  int S, K, n_crit, n_map, n_branches, n_omega;
  validate_mix_inputs(
    yi, mean, sd, omega, crit_yi, bias_indicator, crit_yi_mapping,
    crit_yi_mapping_max, &S, &K, &n_crit, &n_map, &n_branches, &n_omega
  );
  (void)n_omega;
  check_real(log_norm, "log_norm");
  if (matrix_nrow(log_norm, "log_norm") != S ||
      matrix_ncol(log_norm, "log_norm") != K) {
    Rf_error("'log_norm' dimensions must match 'mean'.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const double *yi_p        = REAL(yi);
  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *omega_p     = REAL(omega);
  const double *crit_p      = REAL(crit_yi);
  const int *bias_p         = INTEGER(bias_indicator);
  const double *mapping_p   = REAL(crit_yi_mapping);
  const int *mapping_max_p  = INTEGER(crit_yi_mapping_max);
  const double *log_norm_p  = REAL(log_norm);
  double *out_p             = REAL(out);
  std::vector<int> active_cuts(S);
  std::vector<const double *> mapping_ptr(S);
  prepare_mix_row_metadata(
    bias_p, mapping_max_p, mapping_p, S, n_branches, n_map, n_crit,
    active_cuts, mapping_ptr
  );

  for (int k = 0; k < K; ++k) {
    const double *crit_k = crit_p + n_crit * k;

    for (int s = 0; s < S; ++s) {
      out_p[s + S * k] = cpp_wnorm_mix_lpdf_precomputed(
        yi_p[k],
        mean_p[s + S * k],
        sd_p[s + S * k],
        crit_k,
        omega_p + s,
        mapping_ptr[s],
        active_cuts[s],
        log_norm_p[s + S * k],
        S
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_wnorm_mix_cdf_precomp_matrix(SEXP q, SEXP mean,
                                                    SEXP sd, SEXP omega,
                                                    SEXP crit_yi,
                                                    SEXP bias_indicator,
                                                    SEXP crit_yi_mapping,
                                                    SEXP crit_yi_mapping_max,
                                                    SEXP log_norm,
                                                    SEXP lower_tail,
                                                    SEXP log_p)
{
  int S, K, n_crit, n_map, n_branches, n_omega;
  validate_mix_inputs(
    q, mean, sd, omega, crit_yi, bias_indicator, crit_yi_mapping,
    crit_yi_mapping_max, &S, &K, &n_crit, &n_map, &n_branches, &n_omega
  );
  (void)n_omega;
  check_real(log_norm, "log_norm");
  if (matrix_nrow(log_norm, "log_norm") != S ||
      matrix_ncol(log_norm, "log_norm") != K) {
    Rf_error("'log_norm' dimensions must match 'mean'.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const bool lower = Rf_asLogical(lower_tail) == TRUE;
  const bool log   = Rf_asLogical(log_p) == TRUE;

  const double *q_p         = REAL(q);
  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *omega_p     = REAL(omega);
  const double *crit_p      = REAL(crit_yi);
  const int *bias_p         = INTEGER(bias_indicator);
  const double *mapping_p   = REAL(crit_yi_mapping);
  const int *mapping_max_p  = INTEGER(crit_yi_mapping_max);
  const double *log_norm_p  = REAL(log_norm);
  double *out_p             = REAL(out);
  std::vector<int> active_cuts(S);
  std::vector<const double *> mapping_ptr(S);
  prepare_mix_row_metadata(
    bias_p, mapping_max_p, mapping_p, S, n_branches, n_map, n_crit,
    active_cuts, mapping_ptr
  );

  for (int k = 0; k < K; ++k) {
    const double *crit_k = crit_p + n_crit * k;

    for (int s = 0; s < S; ++s) {
      out_p[s + S * k] = cpp_wnorm_mix_cdf_precomputed(
        q_p[k],
        mean_p[s + S * k],
        sd_p[s + S * k],
        crit_k,
        omega_p + s,
        mapping_ptr[s],
        active_cuts[s],
        log_norm_p[s + S * k],
        lower,
        log,
        S
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_wnorm_mix_moments_matrix(SEXP mean, SEXP sd,
                                                SEXP omega, SEXP crit_yi,
                                                SEXP bias_indicator,
                                                SEXP crit_yi_mapping,
                                                SEXP crit_yi_mapping_max)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(omega, "omega");
  check_real(crit_yi, "crit_yi");
  check_integer(bias_indicator, "bias_indicator");
  check_real(crit_yi_mapping, "crit_yi_mapping");
  check_integer(crit_yi_mapping_max, "crit_yi_mapping_max");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int n_crit = matrix_nrow(crit_yi, "crit_yi");
  const int n_map = matrix_nrow(crit_yi_mapping, "crit_yi_mapping");
  const int n_branches = Rf_length(crit_yi_mapping_max);
  const int n_omega = matrix_ncol(omega, "omega");

  if (matrix_nrow(sd, "sd") != S || matrix_ncol(sd, "sd") != K) {
    Rf_error("'sd' dimensions must match 'mean'.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (matrix_ncol(crit_yi, "crit_yi") != K) {
    Rf_error("'crit_yi' must have one column per observation.");
  }
  if (matrix_ncol(crit_yi_mapping, "crit_yi_mapping") != n_branches) {
    Rf_error("'crit_yi_mapping' columns must match 'crit_yi_mapping_max'.");
  }
  if (Rf_length(bias_indicator) != S) {
    Rf_error("'bias_indicator' must have one value per posterior sample.");
  }
  validate_mix_mapping(
    REAL(crit_yi_mapping), INTEGER(crit_yi_mapping_max),
    n_branches, n_map, n_crit, n_omega
  );

  SEXP out_mean   = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  SEXP out_second = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  SEXP out        = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names      = PROTECT(Rf_allocVector(STRSXP, 2));

  SET_STRING_ELT(names, 0, Rf_mkChar("mean"));
  SET_STRING_ELT(names, 1, Rf_mkChar("second"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *omega_p     = REAL(omega);
  const double *crit_p      = REAL(crit_yi);
  const int *bias_p         = INTEGER(bias_indicator);
  const double *mapping_p   = REAL(crit_yi_mapping);
  const int *mapping_max_p  = INTEGER(crit_yi_mapping_max);
  double *out_mean_p        = REAL(out_mean);
  double *out_second_p      = REAL(out_second);
  std::vector<int> active_cuts(S);
  std::vector<const double *> mapping_ptr(S);
  prepare_mix_row_metadata(
    bias_p, mapping_max_p, mapping_p, S, n_branches, n_map, n_crit,
    active_cuts, mapping_ptr
  );

  for (int k = 0; k < K; ++k) {
    const double *crit_k = crit_p + n_crit * k;

    for (int s = 0; s < S; ++s) {
      cpp_wnorm_mix_moments(
        mean_p[s + S * k],
        sd_p[s + S * k],
        crit_k,
        omega_p + s,
        mapping_ptr[s],
        active_cuts[s],
        out_mean_p + s + S * k,
        out_second_p + s + S * k,
        S
      );
    }
  }

  SET_VECTOR_ELT(out, 0, out_mean);
  SET_VECTOR_ELT(out, 1, out_second);

  UNPROTECT(4);
  return out;
}

extern "C" SEXP RoBMA_wnorm_mix_log_norm_matrix(SEXP mean, SEXP sd,
                                                 SEXP omega, SEXP crit_yi,
                                                 SEXP bias_indicator,
                                                 SEXP crit_yi_mapping,
                                                 SEXP crit_yi_mapping_max)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(omega, "omega");
  check_real(crit_yi, "crit_yi");
  check_integer(bias_indicator, "bias_indicator");
  check_real(crit_yi_mapping, "crit_yi_mapping");
  check_integer(crit_yi_mapping_max, "crit_yi_mapping_max");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int n_crit = matrix_nrow(crit_yi, "crit_yi");
  const int n_map = matrix_nrow(crit_yi_mapping, "crit_yi_mapping");
  const int n_branches = Rf_length(crit_yi_mapping_max);
  const int n_omega = matrix_ncol(omega, "omega");

  if (matrix_nrow(sd, "sd") != S || matrix_ncol(sd, "sd") != K) {
    Rf_error("'sd' dimensions must match 'mean'.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (matrix_ncol(crit_yi, "crit_yi") != K) {
    Rf_error("'crit_yi' must have one column per observation.");
  }
  if (matrix_ncol(crit_yi_mapping, "crit_yi_mapping") != n_branches) {
    Rf_error("'crit_yi_mapping' columns must match 'crit_yi_mapping_max'.");
  }
  if (Rf_length(bias_indicator) != S) {
    Rf_error("'bias_indicator' must have one value per posterior sample.");
  }
  validate_mix_mapping(
    REAL(crit_yi_mapping), INTEGER(crit_yi_mapping_max),
    n_branches, n_map, n_crit, n_omega
  );

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *omega_p     = REAL(omega);
  const double *crit_p      = REAL(crit_yi);
  const int *bias_p         = INTEGER(bias_indicator);
  const double *mapping_p   = REAL(crit_yi_mapping);
  const int *mapping_max_p  = INTEGER(crit_yi_mapping_max);
  double *out_p             = REAL(out);
  std::vector<int> active_cuts(S);
  std::vector<const double *> mapping_ptr(S);
  prepare_mix_row_metadata(
    bias_p, mapping_max_p, mapping_p, S, n_branches, n_map, n_crit,
    active_cuts, mapping_ptr
  );

  for (int k = 0; k < K; ++k) {
    const double *crit_k = crit_p + n_crit * k;

    for (int s = 0; s < S; ++s) {
      out_p[s + S * k] = cpp_wnorm_mix_log_norm(
        mean_p[s + S * k],
        sd_p[s + S * k],
        crit_k,
        omega_p + s,
        mapping_ptr[s],
        active_cuts[s],
        S
      );
    }
  }

  UNPROTECT(1);
  return out;
}
