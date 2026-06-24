context("VIF")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

test_that("VIF rejects data-only objects", {

  dat <- data.frame(
    yi  = c(0.10, 0.20, 0.30, 0.40),
    sei = c(0.10, 0.12, 0.11, 0.15),
    x   = c(0, 1, 0, 1)
  )

  object <- brma.norm(
    yi        = yi,
    sei       = sei,
    mods      = ~ x,
    data      = dat,
    only_data = TRUE
  )

  expect_error(
    vif(object),
    "data-only"
  )
})

.vif_full_gls_vcov_oracle <- function(X, covariance_samples) {

  K <- nrow(X)
  P <- ncol(X)

  if (is.matrix(covariance_samples)) {
    covariance_samples <- array(covariance_samples, dim = c(1L, K, K))
  }

  vcov_sum <- matrix(0, nrow = P, ncol = P)
  for (s in seq_len(dim(covariance_samples)[1L])) {
    M <- covariance_samples[s, , ]
    M <- (M + t(M)) / 2
    W <- solve(M)
    vcov_sum <- vcov_sum + solve(crossprod(X, W %*% X))
  }

  vcov <- vcov_sum / dim(covariance_samples)[1L]
  dimnames(vcov) <- list(colnames(X), colnames(X))

  return(vcov)
}

.vif_known_v_table_oracle <- function(object, covariance_samples) {

  X            <- .get_model_matrix(object)
  expected_cov <- .vif_full_gls_vcov_oracle(X, covariance_samples)

  .vif_table_from_vcov(
    vcov   = expected_cov,
    assign = attr(X, "assign"),
    object = object,
    X      = X
  )
}

test_that("full covariance VIF backend matches manual GLS oracle", {

  X <- cbind(
    intercept = 1,
    x         = c(-1, 0, 1, -1),
    z         = c(0, 1, 1, 2)
  )
  M1 <- matrix(c(
    0.08, 0.02, 0.00, 0.00,
    0.02, 0.09, 0.00, 0.00,
    0.00, 0.00, 0.07, 0.03,
    0.00, 0.00, 0.03, 0.10
  ), nrow = 4, byrow = TRUE)
  M2 <- M1 + diag(c(0.01, 0.02, 0.01, 0.03))

  covariance_samples <- array(NA_real_, dim = c(2, 4, 4))
  covariance_samples[1, , ] <- M1
  covariance_samples[2, , ] <- M2

  expected <- .vif_full_gls_vcov_oracle(X, covariance_samples)
  actual   <- .vif_vcov_from_covariance_samples(
    X                  = X,
    covariance_samples = covariance_samples
  )

  diagonal_covariance <- array(NA_real_, dim = c(2, 4, 4))
  diagonal_covariance[1, , ] <- diag(diag(M1))
  diagonal_covariance[2, , ] <- diag(diag(M2))
  diagonal_only <- .vif_vcov_from_covariance_samples(
    X                  = X,
    covariance_samples = diagonal_covariance
  )

  expect_equal(actual, expected, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(actual, diagonal_only, tolerance = 1e-8)))
})

test_that("known-V VIF uses sampling covariance alone without scalar tau", {

  object <- structure(
    list(
      summary = matrix(
        c(0.1, 0.2),
        nrow     = 2,
        dimnames = list(c("(mu) intercept", "(mu) x"), "Mean")
      ),
      data = structure(list(), scale = FALSE)
    ),
    class = c("brma.mv", "brma.norm", "brma")
  )
  V <- matrix(c(0.04, 0.01, 0.01, 0.05), nrow = 2)

  posterior_samples <- matrix(
    c(0.1, 0.2, 0.3, 0.4),
    nrow = 2,
    dimnames = list(NULL, c("mu", "x"))
  )
  covariance_samples <- .known_v_add_base_covariance(
    base_covariance    = V,
    covariance_samples = .known_v_diagonal_extra_covariance_samples(
      object            = object,
      posterior_samples = posterior_samples,
      K                 = 2
    )
  )

  expect_equal(.vif_scalar_tau_mean(object), 0)
  expect_equal(covariance_samples[1, , ], V)
})

test_that("fixed location coefficient extraction is strict only when requested", {

  samples <- matrix(
    seq_len(8),
    nrow     = 2,
    dimnames = list(
      NULL,
      c("(mu) intercept", "(mu) x", "(mu) var_frac(block)", "(tau) intercept")
    )
  )

  filtered <- .diagnostic_fixed_location_coefficient_samples(
    samples,
    require_mu_columns = TRUE
  )

  expect_equal(colnames(filtered), c("(mu) intercept", "(mu) x"))
  expect_equal(
    .diagnostic_fixed_location_coefficient_samples(samples[, 4, drop = FALSE]),
    samples[, 4, drop = FALSE]
  )
  expect_error(
    .diagnostic_fixed_location_coefficient_samples(
      samples[, 4, drop = FALSE],
      require_mu_columns = TRUE
    ),
    "Fixed location coefficient columns"
  )
  expect_error(
    .diagnostic_fixed_location_coefficient_samples(
      matrix(seq_len(4), nrow = 2),
      require_mu_columns = TRUE
    ),
    "Fixed location coefficient columns"
  )
})

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(vif_cases(), function(case) {
  test_that_case("VIF matches metafor", case, {
    expect_vif_matches_metafor(case)
  })
})

test_that("VIF returns optional posterior correlation diagnostics", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  fit_brma      <- fits[[name]]
  vif_with_post <- vif(fit_brma)
  post_cor      <- vif_with_post[["posterior_correlation"]]

  expect_s3_class(vif_with_post, "vif.brma")
  expect_vif_table(vif_with_post[["vif"]], 2)
  expect_type(post_cor, "double")
  expect_equal(nrow(post_cor), ncol(post_cor))
  expect_equal(unname(diag(post_cor)), rep(1, nrow(post_cor)), tolerance = 1e-12)
  expect_equal(post_cor, t(post_cor), tolerance = 1e-12)
  expect_true(all(is.finite(post_cor)))

  vif_without_post <- vif(fit_brma, posterior_correlation = FALSE)
  expect_null(vif_without_post[["posterior_correlation"]])
})

test_that("VIF rejects models without moderators", {

  model_names <- c(
    "bcg_meta-analysis",
    "konstantopoulos2011_3lvl",
    "nielweise2008_glmm",
    "dat.lehmann2018-PET",
    "dat.lehmann2018-3PSM",
    "brma.mv_block_mvn",
    "brma.mv_block_mvn_random_scale"
  )
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    expect_error(
      vif(fits[[name]]),
      "only meaningful for models with moderators",
      info = paste("vif rejects model without moderators for", name)
    )
  }
})

test_that("VIF supports brma.mv known-V GLS meta-regressions", {

  name <- "brma.mv_block_mvn_mods"
  skip_if_missing_fits(name)

  object  <- fits[[name]]
  covariance <- .known_v_marginal_covariance_samples(object)
  expected   <- .vif_known_v_table_oracle(
    object             = object,
    covariance_samples = covariance
  )
  actual <- vif(object, posterior_correlation = FALSE)[["vif"]]

  expect_vif_table(actual, 2, info = name)
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("VIF supports brma.mv random-formula marginal GLS covariance", {

  name <- "brma.mv_block_mvn_random_mods_scale"
  skip_if_missing_fits(name)

  object      <- fits[[name]]
  known_V     <- .data_known_v_data(object[["data"]])
  random_vcov <- BayesTools::random_effects_marginal_vcov(
    fit       = object[["fit"]],
    parameter = "mu"
  )
  covariance_samples <- .known_v_marginal_covariance_samples(object)

  expected <- .vif_known_v_table_oracle(
    object             = object,
    covariance_samples = covariance_samples
  )
  result <- vif(object)
  actual <- result[["vif"]]
  post_cor <- result[["posterior_correlation"]]

  expect_vif_table(actual, 2, info = name)
  expect_equal(actual, expected, tolerance = 1e-10)
  expect_equal(colnames(post_cor), c("(mu) intercept", "(mu) x", "(mu) z"))
  expect_false(any(grepl("tau", colnames(post_cor), fixed = TRUE)))
  expect_false(any(grepl("var_frac", colnames(post_cor), fixed = TRUE)))
  expect_equal(
    covariance_samples,
    .known_v_add_base_covariance(
      base_covariance    = known_V[["V"]],
      covariance_samples = random_vcov[["samples"]]
    )
  )
  expect_equal(random_vcov[["metadata"]][["n_rows"]], nobs(object))
  expect_true(length(random_vcov[["metadata"]][["included_blocks"]]) > 0)
  expect_equal(nrow(random_vcov[["metadata"]][["skipped_blocks"]]), 0)
  expect_true(any(abs(random_vcov[["samples"]]) > 0))
})

test_that("VIF supports brma.mv known-R marginal GLS covariance", {

  name <- "brma.mv_block_mvn_known_R"
  skip_if_missing_fits(name)

  object      <- fits[[name]]
  known_V     <- .data_known_v_data(object[["data"]])
  random_vcov <- BayesTools::random_effects_marginal_vcov(
    fit       = object[["fit"]],
    parameter = "mu"
  )
  covariance_samples <- .known_v_marginal_covariance_samples(object)

  expected <- .vif_known_v_table_oracle(
    object             = object,
    covariance_samples = covariance_samples
  )
  actual <- vif(object, posterior_correlation = FALSE)[["vif"]]
  group_covariance <- random_vcov[["metadata"]][["blocks"]][["study"]][["group_covariance"]]

  expect_vif_table(actual, 2, info = name)
  expect_equal(actual, expected, tolerance = 1e-10)
  expect_equal(
    covariance_samples,
    .known_v_add_base_covariance(
      base_covariance    = known_V[["V"]],
      covariance_samples = random_vcov[["samples"]]
    )
  )
  expect_equal(random_vcov[["metadata"]][["included_blocks"]], "study")
  expect_equal(group_covariance[["scale"]], "none")
  expect_equal(unname(diag(group_covariance[["kernel"]])), c(4, 9, 16))
  expect_true(any(abs(random_vcov[["samples"]]) > 0))
})

test_that("VIF chunks brma.mv known-V marginal covariance without changing results", {

  name <- "brma.mv_block_mvn_random_mods_scale"
  skip_if_missing_fits(name)

  object       <- fits[[name]]
  expected     <- vif(object, posterior_correlation = FALSE)[["vif"]]
  K            <- nobs(object)
  one_draw_mem <- .known_v_covariance_bytes(1L, K)
  old_options  <- options(
    RoBMA.known_v_covariance_max_bytes = 2 * one_draw_mem
  )
  on.exit(options(old_options), add = TRUE)

  result   <- vif(object, posterior_correlation = FALSE)
  actual   <- result[["vif"]]
  metadata <- attr(result, "known_v_diagnostic")
  attr(actual, "known_v_diagnostic") <- NULL

  expect_equal(actual, expected, tolerance = 1e-10)
  expect_true(metadata[["n_chunks"]] > 1L)
  expect_equal(metadata[["n_used_samples"]], metadata[["n_posterior_samples"]])
})

test_that("VIF exposes deterministic known-V draw thinning", {

  name <- "brma.mv_block_mvn_random_mods_scale"
  skip_if_missing_fits(name)

  object <- fits[[name]]

  result <- expect_warning(
    vif(object, posterior_correlation = FALSE, max_samples = 10),
    "deterministically thinned"
  )
  metadata <- attr(result, "known_v_diagnostic")

  expect_vif_table(result[["vif"]], 2, info = name)
  expect_equal(metadata[["n_used_samples"]], 10L)
  expect_true(metadata[["thinned"]])
})

test_that("Full known-V covariance helper respects allocation budget", {

  name <- "brma.mv_block_mvn_random_mods_scale"
  skip_if_missing_fits(name)

  object            <- fits[[name]]
  posterior_samples <- .get_posterior_samples(object[["fit"]])
  K                 <- nobs(object)
  one_draw_mem      <- .known_v_covariance_bytes(1L, K)

  expect_error(
    .known_v_marginal_covariance_samples(
      object            = object,
      posterior_samples = posterior_samples[seq_len(2L), , drop = FALSE],
      max_bytes         = one_draw_mem
    ),
    "exceeding the configured budget"
  )
})
