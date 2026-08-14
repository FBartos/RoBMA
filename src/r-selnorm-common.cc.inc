// Shared selected-normal native argument contract. R matrices are column-major
// S posterior rows by K observations; selection controls may be scalar or vary
// by posterior row. These structures are non-owning call-lifetime views, and
// validation preserves the R-level kernel semantics.

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
                     SEXP segment_phack_region,
                     bool telescope_probabilities)
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
  data->trusted_step_partition = selnorm_is_descending_step_partition(
    data->z_lower, data->z_upper, B
  );
  data->telescope_probabilities = data->trusted_step_partition &&
    telescope_probabilities;
}

struct SelNormNativeArgs {
  int S;
  int B;
  int n_segments;
  int alpha_len;
  int phack_kind_len;
  int kernel_mode_len;
  const double *omega;
  const double *alpha;
  const int *phack_kind;
  const int *kernel_mode;
  SelNormKernelData data;

  const double *omega_row(int row) const
  {
    return omega + row;
  }

  double alpha_at(int row) const
  {
    return scalar_or_row_real(alpha, alpha_len, row, "alpha");
  }

  int phack_kind_at(int row) const
  {
    return scalar_or_row_int(phack_kind, phack_kind_len, row, "phack_kind");
  }

  int kernel_mode_at(int row) const
  {
    return scalar_or_row_int(kernel_mode, kernel_mode_len, row, "kernel_mode");
  }
};

void validate_selnorm_row_arg_length(SEXP x, int S, const char *name)
{
  if (Rf_length(x) != 1 && Rf_length(x) != S) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
}

SelNormNativeArgs validate_selnorm_native_args(
    SEXP omega, SEXP alpha, SEXP phack_kind, SEXP kernel_mode,
    SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
    SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
    SEXP segment_step_bin, SEXP segment_phack_region,
    SEXP telescope_probabilities, int S)
{
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
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int B = matrix_ncol(omega, "omega");
  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
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

  validate_selnorm_row_arg_length(alpha, S, "alpha");
  validate_selnorm_row_arg_length(phack_kind, S, "phack_kind");
  validate_selnorm_row_arg_length(kernel_mode, S, "kernel_mode");

  SelNormNativeArgs args;
  args.S               = S;
  args.B               = B;
  args.n_segments      = n_segments;
  args.alpha_len       = Rf_length(alpha);
  args.phack_kind_len  = Rf_length(phack_kind);
  args.kernel_mode_len = Rf_length(kernel_mode);
  args.omega           = REAL(omega);
  args.alpha           = REAL(alpha);
  args.phack_kind      = INTEGER(phack_kind);
  args.kernel_mode     = INTEGER(kernel_mode);

  validate_omega_matrix(args.omega, S, B);
  set_kernel_data(
    &args.data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  return args;
}
