#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Random.h>

#include <cmath>
#include <vector>

#include "source/selnorm.h"

extern "C" double Rf_dnorm4(double, double, double, int);

namespace {

int matrix_nrow(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[0];
}

int matrix_ncol(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[1];
}

void check_real(SEXP x, const char *name)
{
  if (TYPEOF(x) != REALSXP) {
    Rf_error("'%s' must be numeric.", name);
  }
}

void check_integer(SEXP x, const char *name)
{
  if (TYPEOF(x) != INTSXP) {
    Rf_error("'%s' must be integer.", name);
  }
}

int scalar_or_row_int(const int *x, int n, int row, const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (row >= n) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
  return x[row];
}

double scalar_or_row_real(const double *x, int n, int row, const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (row >= n) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
  return x[row];
}

void validate_matrix_dims(SEXP x, const char *name, int S, int K)
{
  if (matrix_nrow(x, name) != S || matrix_ncol(x, name) != K) {
    Rf_error("'%s' dimensions must match 'mu_num'.", name);
  }
}

void validate_matrix_dims_named(SEXP x, const char *name, int S, int K,
                                const char *ref)
{
  if (matrix_nrow(x, name) != S || matrix_ncol(x, name) != K) {
    Rf_error("'%s' dimensions must match '%s'.", name, ref);
  }
}

void check_logical(SEXP x, const char *name)
{
  if (TYPEOF(x) != LGLSXP || Rf_length(x) != 1) {
    Rf_error("'%s' must be a single logical value.", name);
  }
}

double scalar_or_col_real(const double *x, int n, int col, const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (col >= n) {
    Rf_error("'%s' must have length 1 or one value per observation.", name);
  }
  return x[col];
}

void validate_omega_matrix(const double *omega, int S, int B)
{
  for (int b = 0; b < B; ++b) {
    for (int s = 0; s < S; ++s) {
      const double w = omega[s + S * b];
      if (!std::isfinite(w) || w < 0) {
        Rf_error("'omega' must contain finite non-negative values.");
      }
    }
  }
}

void set_kernel_data(SelNormKernelData *data, int B, int n_segments,
                     SEXP sign, SEXP q, SEXP z_lower, SEXP z_upper,
                     SEXP phack_z_source, SEXP phack_z_dest,
                     SEXP segment_bounds, SEXP segment_step_bin,
                     SEXP segment_phack_region)
{
  data->n_bins               = B;
  data->n_segments           = n_segments;
  data->effect_sign          = INTEGER(sign)[0];
  data->q                    = INTEGER(q)[0];
  data->z_lower              = REAL(z_lower);
  data->z_upper              = REAL(z_upper);
  data->phack_z_source       = REAL(phack_z_source);
  data->phack_z_dest         = REAL(phack_z_dest);
  data->segment_bounds       = REAL(segment_bounds);
  data->segment_step_bin     = INTEGER(segment_step_bin);
  data->segment_phack_region = INTEGER(segment_phack_region);
  data->segment_step_bin_real = 0;
  data->segment_phack_region_real = 0;
}

int z_step_bin(double z, const SelNormKernelData &data)
{
  for (int b = 0; b < data.n_bins; ++b) {
    if (z >= data.z_lower[b] - 1e-12 && z <= data.z_upper[b] + 1e-12) {
      return b + 1;
    }
  }

  return z > data.z_upper[0] ? 1 : data.n_bins;
}

bool row_has_active_phack(int mode, double alpha, int phack_kind)
{
  return (mode == SELKERNEL_PHACK_POWER ||
    mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
}

}

extern "C" SEXP RoBMA_selnorm_kernel_loglik_matrix(
  SEXP yi, SEXP mu_num, SEXP sigma_num, SEXP mu_norm, SEXP sigma_norm,
  SEXP sei, SEXP weights, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP q, SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region)
{
  check_real(yi, "yi");
  check_real(mu_num, "mu_num");
  check_real(sigma_num, "sigma_num");
  check_real(mu_norm, "mu_norm");
  check_real(sigma_norm, "sigma_norm");
  check_real(sei, "sei");
  check_real(weights, "weights");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(obs_bin, "obs_bin");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");

  const int S = matrix_nrow(mu_num, "mu_num");
  const int K = matrix_ncol(mu_num, "mu_num");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims(sigma_num, "sigma_num", S, K);
  validate_matrix_dims(mu_norm, "mu_norm", S, K);
  validate_matrix_dims(sigma_norm, "sigma_norm", S, K);

  if (Rf_length(yi) != K || Rf_length(sei) != K || Rf_length(weights) != K ||
      Rf_length(obs_bin) != K) {
    Rf_error("'yi', 'sei', 'weights', and 'obs_bin' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region
  );

  const double *yi_p         = REAL(yi);
  const double *mu_num_p     = REAL(mu_num);
  const double *sigma_num_p  = REAL(sigma_num);
  const double *mu_norm_p    = REAL(mu_norm);
  const double *sigma_norm_p = REAL(sigma_norm);
  const double *sei_p        = REAL(sei);
  const double *weights_p    = REAL(weights);
  const double *omega_p      = REAL(omega);
  const double *alpha_p      = REAL(alpha);
  const int *phack_kind_p    = INTEGER(phack_kind);
  const int *kernel_mode_p   = INTEGER(kernel_mode);
  const int *obs_bin_p       = INTEGER(obs_bin);
  double *out_p              = REAL(out);

  validate_omega_matrix(omega_p, S, B);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_lpdf(
        yi_p[k],
        mu_num_p[cell],
        sigma_num_p[cell],
        mu_norm_p[cell],
        sigma_norm_p[cell],
        sei_p[k],
        weights_p[k],
        omega_p + s,
        obs_bin_p[k],
        scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha"),
        scalar_or_row_int(phack_kind_p, Rf_length(phack_kind), s, "phack_kind"),
        scalar_or_row_int(kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"),
        data,
        S,
        false
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_log_norm_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region
  );

  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *sei_p       = REAL(sei);
  const double *omega_p     = REAL(omega);
  const double *alpha_p     = REAL(alpha);
  const int *phack_kind_p   = INTEGER(phack_kind);
  const int *kernel_mode_p  = INTEGER(kernel_mode);
  double *out_p             = REAL(out);

  validate_omega_matrix(omega_p, S, B);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_log_norm(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha"),
        scalar_or_row_int(phack_kind_p, Rf_length(phack_kind), s, "phack_kind"),
        scalar_or_row_int(kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"),
        data,
        S,
        false
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_cdf_matrix(
  SEXP x, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP lower_tail)
{
  check_real(x, "q");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(lower_tail, "lower.tail");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if ((Rf_length(x) != 1 && Rf_length(x) != K) || Rf_length(sei) != K) {
    Rf_error("'q' must have length 1 or K and 'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region
  );

  const double *x_p          = REAL(x);
  const double *mean_p       = REAL(mean);
  const double *sd_p         = REAL(sd);
  const double *sei_p        = REAL(sei);
  const double *omega_p      = REAL(omega);
  const double *alpha_p      = REAL(alpha);
  const int *phack_kind_p    = INTEGER(phack_kind);
  const int *kernel_mode_p   = INTEGER(kernel_mode);
  const bool lower_tail_bool = LOGICAL(lower_tail)[0] == TRUE;
  double *out_p              = REAL(out);

  validate_omega_matrix(omega_p, S, B);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_cdf(
        scalar_or_col_real(x_p, Rf_length(x), k, "q"),
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha"),
        scalar_or_row_int(phack_kind_p, Rf_length(phack_kind), s, "phack_kind"),
        scalar_or_row_int(kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"),
        data,
        S,
        lower_tail_bool,
        false
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_moments_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out_mean = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  SEXP out_second = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));

  SET_STRING_ELT(names, 0, Rf_mkChar("mean"));
  SET_STRING_ELT(names, 1, Rf_mkChar("second"));
  SET_VECTOR_ELT(out, 0, out_mean);
  SET_VECTOR_ELT(out, 1, out_second);
  Rf_setAttrib(out, R_NamesSymbol, names);

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region
  );

  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  const double *omega_p    = REAL(omega);
  const double *alpha_p    = REAL(alpha);
  const int *phack_kind_p  = INTEGER(phack_kind);
  const int *kernel_mode_p = INTEGER(kernel_mode);
  double *mean_out_p       = REAL(out_mean);
  double *second_out_p     = REAL(out_second);

  validate_omega_matrix(omega_p, S, B);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      cpp_selnorm_kernel_moments(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha"),
        scalar_or_row_int(phack_kind_p, Rf_length(phack_kind), s, "phack_kind"),
        scalar_or_row_int(kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"),
        data,
        mean_out_p + cell,
        second_out_p + cell,
        S,
        false
      );
    }
  }

  UNPROTECT(4);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_rng_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }
  if (Rf_length(alpha) != 1 && Rf_length(alpha) != S) {
    Rf_error("'alpha' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(phack_kind) != 1 && Rf_length(phack_kind) != S) {
    Rf_error("'phack_kind' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(kernel_mode) != 1 && Rf_length(kernel_mode) != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region
  );

  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  const double *omega_p    = REAL(omega);
  const double *alpha_p    = REAL(alpha);
  const int *phack_kind_p  = INTEGER(phack_kind);
  const int *kernel_mode_p = INTEGER(kernel_mode);
  double *out_p            = REAL(out);

  validate_omega_matrix(omega_p, S, B);

  GetRNGstate();
  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_rng(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        unif_rand(),
        unif_rand(),
        scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha"),
        scalar_or_row_int(phack_kind_p, Rf_length(phack_kind), s, "phack_kind"),
        scalar_or_row_int(kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"),
        data,
        S,
        false
      );
    }
  }
  PutRNGstate();

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_weighted_summary(
  SEXP yi, SEXP mean, SEXP sd, SEXP sei, SEXP psis_weights, SEXP omega,
  SEXP alpha, SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region)
{
  check_real(yi, "yi");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(psis_weights, "psis_weights");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  validate_matrix_dims_named(psis_weights, "psis_weights", S, K, "mean");

  if (Rf_length(yi) != K || Rf_length(sei) != K) {
    Rf_error("'yi' and 'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }
  if (Rf_length(alpha) != 1 && Rf_length(alpha) != S) {
    Rf_error("'alpha' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(phack_kind) != 1 && Rf_length(phack_kind) != S) {
    Rf_error("'phack_kind' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(kernel_mode) != 1 && Rf_length(kernel_mode) != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out_cdf = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP out_mean = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP out_second = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));

  SET_STRING_ELT(names, 0, Rf_mkChar("cdf"));
  SET_STRING_ELT(names, 1, Rf_mkChar("mean"));
  SET_STRING_ELT(names, 2, Rf_mkChar("second"));
  SET_VECTOR_ELT(out, 0, out_cdf);
  SET_VECTOR_ELT(out, 1, out_mean);
  SET_VECTOR_ELT(out, 2, out_second);
  Rf_setAttrib(out, R_NamesSymbol, names);

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region
  );

  const double *yi_p        = REAL(yi);
  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *sei_p       = REAL(sei);
  const double *weights_p   = REAL(psis_weights);
  const double *omega_p     = REAL(omega);
  const double *alpha_p     = REAL(alpha);
  const int *phack_kind_p   = INTEGER(phack_kind);
  const int *kernel_mode_p  = INTEGER(kernel_mode);
  double *cdf_out_p         = REAL(out_cdf);
  double *mean_out_p        = REAL(out_mean);
  double *second_out_p      = REAL(out_second);

  validate_omega_matrix(omega_p, S, B);

  for (int k = 0; k < K; ++k) {
    double cdf_sum = 0;
    double mean_sum = 0;
    double second_sum = 0;

    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      const double w = weights_p[cell];
      const double alpha_s = scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha");
      const int phack_s = scalar_or_row_int(
        phack_kind_p, Rf_length(phack_kind), s, "phack_kind"
      );
      const int mode_s = scalar_or_row_int(
        kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"
      );
      double moment_mean = NA_REAL;
      double moment_second = NA_REAL;

      cpp_selnorm_kernel_moments(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        alpha_s,
        phack_s,
        mode_s,
        data,
        &moment_mean,
        &moment_second,
        S,
        false
      );

      cdf_sum += w * cpp_selnorm_kernel_cdf(
        yi_p[k],
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        alpha_s,
        phack_s,
        mode_s,
        data,
        S,
        true,
        false
      );
      mean_sum += w * moment_mean;
      second_sum += w * moment_second;
    }

    cdf_out_p[k] = cdf_sum;
    mean_out_p[k] = mean_sum;
    second_out_p[k] = second_sum;
  }

  UNPROTECT(5);
  return out;
}

extern "C" SEXP RoBMA_zcurve_normal_density_matrix(
  SEXP z_sequence, SEXP mean, SEXP sd, SEXP sei)
{
  check_real(z_sequence, "z_sequence");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int L = Rf_length(z_sequence);

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, L));

  const double *z_p    = REAL(z_sequence);
  const double *mean_p = REAL(mean);
  const double *sd_p   = REAL(sd);
  const double *sei_p  = REAL(sei);
  double *out_p        = REAL(out);

  for (int ell = 0; ell < L; ++ell) {
    for (int s = 0; s < S; ++s) {
      double dens_sum = 0;
      for (int k = 0; k < K; ++k) {
        const int cell = s + S * k;
        dens_sum += Rf_dnorm4(
          z_p[ell] * sei_p[k], mean_p[cell], sd_p[cell], false
        ) * sei_p[k];
      }
      out_p[s + S * ell] = dens_sum / K;
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_zcurve_density_matrix(
  SEXP z_sequence, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP extrapolate)
{
  check_real(z_sequence, "z_sequence");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(extrapolate, "extrapolate");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int L = Rf_length(z_sequence);
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }
  if (Rf_length(alpha) != 1 && Rf_length(alpha) != S) {
    Rf_error("'alpha' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(phack_kind) != 1 && Rf_length(phack_kind) != S) {
    Rf_error("'phack_kind' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(kernel_mode) != 1 && Rf_length(kernel_mode) != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, L));

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region
  );

  const double *z_p        = REAL(z_sequence);
  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  const double *omega_p    = REAL(omega);
  const double *alpha_p    = REAL(alpha);
  const int *phack_kind_p  = INTEGER(phack_kind);
  const int *kernel_mode_p = INTEGER(kernel_mode);
  const bool extrapolate_b = LOGICAL(extrapolate)[0] == TRUE;
  double *out_p            = REAL(out);

  validate_omega_matrix(omega_p, S, B);

  std::vector<double> log_norm(S * K);
  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      const double alpha_s = scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha");
      const int phack_s = scalar_or_row_int(
        phack_kind_p, Rf_length(phack_kind), s, "phack_kind"
      );
      const int mode_s = scalar_or_row_int(
        kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"
      );

      if (row_has_active_phack(mode_s, alpha_s, phack_s)) {
        log_norm[cell] = std::numeric_limits<double>::quiet_NaN();
      } else {
        log_norm[cell] = cpp_selnorm_kernel_log_norm(
          mean_p[cell],
          sd_p[cell],
          sei_p[k],
          omega_p + s,
          alpha_s,
          phack_s,
          mode_s,
          data,
          S,
          false
        );
      }
    }
  }

  for (int ell = 0; ell < L; ++ell) {
    for (int s = 0; s < S; ++s) {
      const double alpha_s = scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha");
      const int phack_s = scalar_or_row_int(
        phack_kind_p, Rf_length(phack_kind), s, "phack_kind"
      );
      const int mode_s = scalar_or_row_int(
        kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"
      );
      const bool use_step = mode_s == SELKERNEL_STEP ||
        mode_s == SELKERNEL_STEP_PHACK_POWER;

      double dens_sum = 0;
      for (int k = 0; k < K; ++k) {
        const int cell = s + S * k;
        double log_density = Rf_dnorm4(
          z_p[ell] * sei_p[k], mean_p[cell], sd_p[cell], true
        );

        if (use_step || extrapolate_b) {
          log_density -= log_norm[cell];
        }

        if (use_step && !extrapolate_b) {
          const int bin = z_step_bin(data.effect_sign * z_p[ell], data);
          const double w = omega_p[s + S * (bin - 1)];
          if (w > 0) {
            log_density += std::log(w);
          } else {
            log_density = -INFINITY;
          }
        }

        if (row_has_active_phack(mode_s, alpha_s, phack_s)) {
          log_density = std::numeric_limits<double>::quiet_NaN();
        }

        dens_sum += std::exp(log_density) * sei_p[k];
      }
      out_p[s + S * ell] = dens_sum / K;
    }
  }

  UNPROTECT(1);
  return out;
}
