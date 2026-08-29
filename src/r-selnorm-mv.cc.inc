extern "C" SEXP RoBMA_selnorm_mnorm_step_loglik_batch(
  SEXP yi, SEXP means, SEXP covariance_lower, SEXP sei, SEXP omega,
  SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP telescope_probabilities, SEXP kernel_mode, SEXP qmc,
  SEXP points, SEXP scrambles, SEXP relative_tolerance)
{
  check_real(yi, "yi");
  check_real(means, "means");
  check_real(covariance_lower, "covariance_lower");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(obs_bin, "obs_bin");
  check_integer(sign, "sign");
  check_logical(telescope_probabilities, "telescope_probabilities");
  check_integer(kernel_mode, "kernel_mode");
  check_real(qmc, "qmc");
  check_integer(points, "points");
  check_integer(scrambles, "scrambles");
  check_real(relative_tolerance, "relative_tolerance");

  const int S = matrix_nrow(means, "means");
  const int K = matrix_ncol(means, "means");
  const int L = K * (K + 1) / 2;
  const int B = matrix_ncol(omega, "omega");
  if (S < 1 || K < 1) {
    Rf_error("'means' must have positive dimensions.");
  }
  validate_matrix_dims_named(
    covariance_lower, "covariance_lower", S, L, "means"
  );
  if (matrix_nrow(omega, "omega") != S || B < 1) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(yi) != K || Rf_length(sei) != K ||
      Rf_length(obs_bin) != K) {
    Rf_error("'yi', 'sei', and 'obs_bin' must have one value per observation.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("Selection bounds must have one value per selection bin.");
  }
  if (Rf_length(sign) != 1 ||
      (INTEGER(sign)[0] != -1 && INTEGER(sign)[0] != 1)) {
    Rf_error("'sign' must be either -1 or 1.");
  }
  const int kernel_length = Rf_length(kernel_mode);
  if (kernel_length != 1 && kernel_length != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(points) != 1 || INTEGER(points)[0] < 1 ||
      Rf_length(scrambles) != 1 || INTEGER(scrambles)[0] < 2) {
    Rf_error("'points' and 'scrambles' contain invalid integration settings.");
  }
  const int n_points    = INTEGER(points)[0];
  const int n_scrambles = INTEGER(scrambles)[0];
  const R_xlen_t expected_qmc = static_cast<R_xlen_t>(n_scrambles) *
    static_cast<R_xlen_t>(n_points) * static_cast<R_xlen_t>(2 * K);
  if (XLENGTH(qmc) != expected_qmc) {
    Rf_error("'qmc' dimensions do not match the integration settings.");
  }
  if (Rf_length(relative_tolerance) != 1 ||
      !std::isfinite(REAL(relative_tolerance)[0]) ||
      REAL(relative_tolerance)[0] <= 0.0) {
    Rf_error("'relative_tolerance' must be one finite positive value.");
  }

  validate_omega_matrix(REAL(omega), S, B);
  for (int k = 0; k < K; ++k) {
    if (!std::isfinite(REAL(yi)[k]) || !std::isfinite(REAL(sei)[k]) ||
        REAL(sei)[k] <= 0.0 || INTEGER(obs_bin)[k] < 1 ||
        INTEGER(obs_bin)[k] > B) {
      Rf_error("Observed selection inputs are invalid.");
    }
  }

  SEXP log_density = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP relative_mcse = PROTECT(Rf_allocVector(REALSXP, S));
  std::vector<double> mean(static_cast<std::size_t>(K));
  std::vector<double> lower(static_cast<std::size_t>(L));
  std::vector<double> omega_row(static_cast<std::size_t>(B));

  for (int s = 0; s < S; ++s) {
    for (int k = 0; k < K; ++k) {
      mean[static_cast<std::size_t>(k)] = REAL(means)[s + S * k];
    }
    for (int l = 0; l < L; ++l) {
      lower[static_cast<std::size_t>(l)] =
        REAL(covariance_lower)[s + S * l];
    }
    for (int b = 0; b < B; ++b) {
      omega_row[static_cast<std::size_t>(b)] = REAL(omega)[s + S * b];
    }

    double error = 0.0;
    REAL(log_density)[s] = cpp_selnorm_mnorm_step_lpdf(
      REAL(yi), mean.data(), lower.data(), K, REAL(sei), omega_row.data(), B,
      REAL(z_lower), REAL(z_upper), INTEGER(obs_bin), INTEGER(sign)[0],
      LOGICAL(telescope_probabilities)[0],
      scalar_or_row_int(INTEGER(kernel_mode), kernel_length, s, "kernel_mode"),
      REAL(qmc), n_points, n_scrambles, &error
    );
    REAL(relative_mcse)[s] = error;
  }

  SEXP output = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names  = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_VECTOR_ELT(output, 0, log_density);
  SET_VECTOR_ELT(output, 1, relative_mcse);
  SET_STRING_ELT(names, 0, Rf_mkChar("log_density"));
  SET_STRING_ELT(names, 1, Rf_mkChar("relative_mcse"));
  Rf_setAttrib(output, R_NamesSymbol, names);

  UNPROTECT(4);
  return output;
}

namespace {

bool selnorm_mvn_exact_rank_one_factor(const std::vector<double> &covariance,
                                       int n, double *factor,
                                       double *rank_one)
{
  int pivot = 0;
  double largest = -1.0;
  for (int i = 0; i < n; ++i) {
    const double diagonal = covariance[static_cast<std::size_t>(i + n * i)];
    if (!(diagonal > 0.0)) return false;
    const double magnitude = std::sqrt(diagonal);
    if (magnitude > largest) {
      largest = magnitude;
      pivot = i;
    }
  }

  for (int i = 0; i < n; ++i) {
    const double pivot_value = covariance[static_cast<std::size_t>(pivot + n * i)];
    if (pivot_value == 0.0) return false;
    const double direction = i == pivot ? 1.0 : (pivot_value > 0.0 ? 1.0 : -1.0);
    rank_one[i] =
      direction * std::sqrt(covariance[static_cast<std::size_t>(i + n * i)]);
  }
  for (int column = 0; column < n; ++column) {
    for (int row = 0; row < n; ++row) {
      if (rank_one[row] * rank_one[column] !=
          covariance[static_cast<std::size_t>(row + n * column)]) {
        return false;
      }
    }
  }

  std::fill(factor, factor + static_cast<std::size_t>(n * n), 0.0);
  for (int column = 0; column < n; ++column) {
    factor[n * column] = rank_one[column];
  }
  return true;
}

bool selnorm_mvn_sampling_factor(const std::vector<double> &covariance,
                                 int n, double *factor,
                                 double *eigenvalues, double *work,
                                 int work_length, double *rank_one,
                                 double *cholesky, double *eigenvectors)
{
  if (selnorm_mvn_exact_rank_one_factor(
        covariance, n, factor, rank_one)) {
    return true;
  }

  std::copy(covariance.begin(), covariance.end(), eigenvectors);
  int info = 0;
  F77_CALL(dsyev)("V", "U", &n, eigenvectors, &n, eigenvalues,
                  work, &work_length, &info FCONE FCONE);
  if (info != 0) return false;

  double scale = 0.0;
  for (int i = 0; i < n; ++i) {
    scale = std::max(scale, std::abs(eigenvalues[i]));
  }
  const double operations = 4.0 * static_cast<double>(std::max(1, n));
  const double scaled_epsilon = operations *
    std::numeric_limits<double>::epsilon();
  const double tolerance = scale == 0.0 ? 0.0 :
    scaled_epsilon / (1.0 - scaled_epsilon) * scale;
  if (eigenvalues[0] < -tolerance) return false;

  std::copy(covariance.begin(), covariance.end(), cholesky);
  F77_CALL(dpotrf)("U", &n, cholesky, &n, &info FCONE);
  if (info == 0) {
    std::fill(factor, factor + static_cast<std::size_t>(n * n), 0.0);
    for (int column = 0; column < n; ++column) {
      for (int row = 0; row <= column; ++row) {
        factor[row + n * column] =
          cholesky[row + n * column];
      }
    }
    return true;
  }

  for (int i = 0; i < n; ++i) {
    if (std::abs(eigenvalues[i]) <= tolerance) eigenvalues[i] = 0.0;
    if (eigenvalues[i] < 0.0) return false;
  }
  for (int column = 0; column < n; ++column) {
    for (int row = 0; row < n; ++row) {
      double value = 0.0;
      for (int component = 0; component < n; ++component) {
        value += eigenvectors[row + n * component] *
          std::sqrt(eigenvalues[component]) *
          eigenvectors[column + n * component];
      }
      factor[row + n * column] = value;
    }
  }
  return true;
}

int selnorm_mvn_step_bin(double value, double sei, int sign,
                         const double *z_lower, int n_bins)
{
  const double score = static_cast<double>(sign) * value / sei;
  for (int bin = 0; bin < n_bins; ++bin) {
    if (score >= z_lower[bin]) return bin;
  }
  return -1;
}

}

extern "C" SEXP RoBMA_selnorm_mnorm_step_rng_batch(
  SEXP means, SEXP covariance, SEXP sei, SEXP omega,
  SEXP z_lower, SEXP z_upper, SEXP sign, SEXP kernel_mode,
  SEXP dependency_blocks, SEXP max_attempts)
{
  check_real(means, "means");
  check_real(covariance, "covariance");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(kernel_mode, "kernel_mode");
  check_integer(max_attempts, "max_attempts");

  const int S = matrix_nrow(means, "means");
  const int K = matrix_ncol(means, "means");
  const int B = matrix_ncol(omega, "omega");
  if (S < 1 || K < 1 || B < 1 || matrix_nrow(omega, "omega") != S) {
    Rf_error("Selected response inputs have invalid dimensions.");
  }
  SEXP covariance_dim = Rf_getAttrib(covariance, R_DimSymbol);
  if (TYPEOF(covariance_dim) != INTSXP || Rf_length(covariance_dim) != 3 ||
      INTEGER(covariance_dim)[0] != S ||
      INTEGER(covariance_dim)[1] != K ||
      INTEGER(covariance_dim)[2] != K) {
    Rf_error("'covariance' dimensions must match 'means'.");
  }
  if (Rf_length(sei) != K || Rf_length(z_lower) != B ||
      Rf_length(z_upper) != B) {
    Rf_error("Selected response vectors have invalid dimensions.");
  }
  if (Rf_length(sign) != 1 ||
      (INTEGER(sign)[0] != -1 && INTEGER(sign)[0] != 1)) {
    Rf_error("'sign' must be either -1 or 1.");
  }
  const int mode_length = Rf_length(kernel_mode);
  if (mode_length != 1 && mode_length != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(max_attempts) != 1 || INTEGER(max_attempts)[0] < 1) {
    Rf_error("'max_attempts' must be one positive integer.");
  }
  if (TYPEOF(dependency_blocks) != VECSXP || Rf_length(dependency_blocks) < 1) {
    Rf_error("'dependency_blocks' must be a non-empty list.");
  }
  if (!selnorm_is_descending_step_partition(
        REAL(z_lower), REAL(z_upper), B)) {
    Rf_error("Selection bounds do not define a descending step partition.");
  }

  validate_omega_matrix(REAL(omega), S, B);
  std::vector<int> seen(static_cast<std::size_t>(K), 0);
  int maximum_block_size = 0;
  for (int block_index = 0;
       block_index < Rf_length(dependency_blocks); ++block_index) {
    SEXP block = VECTOR_ELT(dependency_blocks, block_index);
    if (TYPEOF(block) != INTSXP || Rf_length(block) < 1) {
      Rf_error("'dependency_blocks' contains an invalid block.");
    }
    maximum_block_size = std::max(maximum_block_size, Rf_length(block));
    for (int j = 0; j < Rf_length(block); ++j) {
      const int row = INTEGER(block)[j] - 1;
      if (row < 0 || row >= K || seen[static_cast<std::size_t>(row)] != 0) {
        Rf_error("'dependency_blocks' must partition the observations.");
      }
      seen[static_cast<std::size_t>(row)] = 1;
    }
  }
  for (int k = 0; k < K; ++k) {
    if (seen[static_cast<std::size_t>(k)] != 1) {
      Rf_error("'dependency_blocks' must partition the observations.");
    }
    if (!std::isfinite(REAL(sei)[k]) || REAL(sei)[k] <= 0.0) {
      Rf_error("'sei' must contain finite positive values.");
    }
  }
  for (R_xlen_t i = 0; i < XLENGTH(means); ++i) {
    if (!std::isfinite(REAL(means)[i])) {
      Rf_error("'means' must contain finite values.");
    }
  }
  for (R_xlen_t i = 0; i < XLENGTH(covariance); ++i) {
    if (!std::isfinite(REAL(covariance)[i])) {
      Rf_error("'covariance' must contain finite values.");
    }
  }
  for (int s = 0; s < S; ++s) {
    const int mode = scalar_or_row_int(
      INTEGER(kernel_mode), mode_length, s, "kernel_mode"
    );
    if (mode != SELKERNEL_NORMAL && mode != SELKERNEL_STEP) {
      Rf_error("Exact selected response simulation requires a step selection kernel.");
    }
  }

  SEXP draws = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  std::fill(REAL(draws), REAL(draws) + static_cast<R_xlen_t>(S) * K,
            NA_REAL);
  std::vector<double> block_covariance(
    static_cast<std::size_t>(maximum_block_size * maximum_block_size)
  );
  std::vector<double> factor(
    static_cast<std::size_t>(maximum_block_size * maximum_block_size)
  );
  std::vector<double> eigenvalues(static_cast<std::size_t>(maximum_block_size));
  const int work_length = std::max(1, 3 * maximum_block_size - 1);
  std::vector<double> work(static_cast<std::size_t>(work_length));
  std::vector<double> rank_one(static_cast<std::size_t>(maximum_block_size));
  std::vector<double> cholesky(
    static_cast<std::size_t>(maximum_block_size * maximum_block_size)
  );
  std::vector<double> eigenvectors(
    static_cast<std::size_t>(maximum_block_size * maximum_block_size)
  );
  std::vector<double> normal(static_cast<std::size_t>(maximum_block_size));
  std::vector<double> candidate(static_cast<std::size_t>(maximum_block_size));

  int failure_code  = 0;
  int failure_size  = 0;

  GetRNGstate();
  for (int s = 0; s < S && failure_code == 0; ++s) {
    double maximum_weight = 0.0;
    for (int b = 0; b < B; ++b) {
      maximum_weight = std::max(maximum_weight, REAL(omega)[s + S * b]);
    }
    if (!(maximum_weight > 0.0)) {
      failure_code = 3;
      break;
    }
    const int mode = scalar_or_row_int(
      INTEGER(kernel_mode), mode_length, s, "kernel_mode"
    );

    for (int block_index = 0;
         block_index < Rf_length(dependency_blocks) && failure_code == 0;
         ++block_index) {
      SEXP block = VECTOR_ELT(dependency_blocks, block_index);
      const int n = Rf_length(block);
      block_covariance.resize(static_cast<std::size_t>(n * n));
      for (int column = 0; column < n; ++column) {
        const int source_column = INTEGER(block)[column] - 1;
        for (int row = 0; row < n; ++row) {
          const int source_row = INTEGER(block)[row] - 1;
          const double a = REAL(covariance)[
            s + S * (source_row + K * source_column)
          ];
          const double b = REAL(covariance)[
            s + S * (source_column + K * source_row)
          ];
          const double scale = std::max(std::abs(a), std::abs(b));
          const double operations =
            4.0 * static_cast<double>(std::max(1, n));
          const double scaled_epsilon = operations *
            std::numeric_limits<double>::epsilon();
          const double tolerance = scale == 0.0 ? 0.0 :
            scaled_epsilon / (1.0 - scaled_epsilon) * scale;
          if (std::abs(a - b) > tolerance) {
            failure_code  = 4;
            failure_size  = n;
            break;
          }
          block_covariance[static_cast<std::size_t>(row + n * column)] =
            b;
        }
        if (failure_code != 0) break;
      }
      if (failure_code != 0) break;
      if (!selnorm_mvn_sampling_factor(
            block_covariance, n, factor.data(), eigenvalues.data(),
            work.data(), work_length, rank_one.data(), cholesky.data(),
            eigenvectors.data())) {
        failure_code  = 1;
        failure_size  = n;
        break;
      }

      bool accepted = false;
      for (int attempt = 0;
           attempt < INTEGER(max_attempts)[0] && !accepted; ++attempt) {
        for (int row = 0; row < n; ++row) normal[row] = norm_rand();
        for (int column = 0; column < n; ++column) {
          double value = REAL(means)[
            s + S * (INTEGER(block)[column] - 1)
          ];
          for (int row = 0; row < n; ++row) {
            value += normal[row] * factor[row + n * column];
          }
          candidate[column] = value;
        }

        accepted = mode == SELKERNEL_NORMAL;
        if (!accepted) {
          double log_acceptance = 0.0;
          for (int row = 0; row < n; ++row) {
            const int observation = INTEGER(block)[row] - 1;
            const int bin = selnorm_mvn_step_bin(
              candidate[row], REAL(sei)[observation], INTEGER(sign)[0],
              REAL(z_lower), B
            );
            if (bin < 0) {
              log_acceptance = -std::numeric_limits<double>::infinity();
              break;
            }
            log_acceptance += std::log(
              REAL(omega)[s + S * bin] / maximum_weight
            );
          }
          accepted = std::log(unif_rand()) <= log_acceptance;
        }
        if (accepted) {
          for (int row = 0; row < n; ++row) {
            REAL(draws)[s + S * (INTEGER(block)[row] - 1)] = candidate[row];
          }
        }
      }
      if (!accepted) {
        failure_code  = 2;
        failure_size  = n;
      }
    }
  }
  PutRNGstate();

  SEXP output = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names  = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_VECTOR_ELT(output, 0, draws);
  SET_VECTOR_ELT(output, 1, Rf_ScalarInteger(failure_code));
  SET_VECTOR_ELT(output, 2, Rf_ScalarInteger(failure_size));
  SET_STRING_ELT(names, 0, Rf_mkChar("draws"));
  SET_STRING_ELT(names, 1, Rf_mkChar("failure_code"));
  SET_STRING_ELT(names, 2, Rf_mkChar("failure_size"));
  Rf_setAttrib(output, R_NamesSymbol, names);

  UNPROTECT(3);
  return output;
}
