.covariance_underflow_fixture <- function() {

  smallest <- .Machine$double.xmin * .Machine$double.eps
  off_diagonal <- .8 * sqrt(smallest)
  matrix(c(1, off_diagonal, off_diagonal, smallest), 2L)
}

.covariance_cholesky_unavailable_message <- function(context) {

  paste0(context, " is positive definite, but a usable Cholesky factor is ",
         "unavailable at working precision. Rescale the outcome, covariance, ",
         "and associated priors consistently before fitting.")
}

test_that("positive-definite classification survives an unrepresentable Cholesky pivot", {

  covariance <- .covariance_underflow_fixture()
  # Positive diagonal and a correlation strictly inside (-1,1) establish SPD
  # without forming the determinant that is smaller than any positive double.
  expect_true(all(diag(covariance) > 0))
  expect_lt(abs(covariance[1L, 2L]) / sqrt(covariance[2L, 2L]), 1)
  factorization <- .covariance_factorization(covariance)
  expect_identical(factorization[["status"]], "positive_definite")
  expect_false(factorization[["singular"]])
  expect_equal(nrow(factorization[["sampling_factor"]]), 2L)
  expect_null(factorization[["cholesky"]])
  error <- tryCatch(.covariance_cholesky(factorization), error = identity)
  expect_identical(conditionMessage(error), .covariance_cholesky_unavailable_message("Covariance"))
  expect_null(conditionCall(error))

  # A subnormal variance with a representable pivot still works unchanged.
  diagonal <- diag(diag(covariance))
  expect_identical(.covariance_cholesky(.covariance_factorization(diagonal)), chol(diagonal))
})

test_that("sampling keeps the accepted full-rank factor without requiring Cholesky", {

  covariance <- .covariance_underflow_fixture()
  factorization <- .covariance_factorization(covariance)
  factor <- .covariance_sampling_factor(factorization)
  expect_identical(factor, factorization[["sampling_factor"]])
  expect_true(all(is.finite(factor)))
  withr::local_seed(27)
  expected <- matrix(stats::rnorm(6L), 3L, 2L, byrow = TRUE) %*% factor
  set.seed(27)
  actual <- .outcome_rng.norm_known_v_covariance(matrix(0, 3L, 2L), covariance)
  expect_identical(actual, expected)
  expect_true(all(is.finite(actual)))
  expect_true(all(actual[, 2L] != 0))
})

test_that("likelihood and projection boundaries identify unavailable Cholesky factors", {

  covariance <- .covariance_underflow_fixture()
  expect_error(.known_v_chol_covariance(covariance, "conditional"),
               .covariance_cholesky_unavailable_message("Known-V conditional covariance"), fixed = TRUE)
  expect_error(.known_v_is_numerically_positive_definite(covariance),
               .covariance_cholesky_unavailable_message("Known-V block-MVN covariance"), fixed = TRUE)
  expect_error(.vif_vcov_from_covariance_samples(matrix(1, 2L, 1L), array(covariance, c(1L, 2L, 2L))),
               .covariance_cholesky_unavailable_message("VIF marginal covariance"), fixed = TRUE)
  known_V <- .known_v_prepare(covariance, c(TRUE, TRUE), "block_mvn", warn_singular = FALSE)
  expect_error(.known_v_gls_projection_blocks(matrix(1, 2L, 1L), c(0, 0), known_V, c(0, 0)),
               .covariance_cholesky_unavailable_message("Known-V residual covariance"), fixed = TRUE)
  expect_error(.zplot_gaussian_conditional_variance(matrix(0, 2L, 2L), covariance),
               .covariance_cholesky_unavailable_message("Estimate-depth marginal covariance"), fixed = TRUE)
  expect_error(.evaluate.brma.known_v_blup.norm(matrix(0, 1L, 2L), matrix(0, 1L, 2L),
                                               c(0, 0), known_V),
               .covariance_cholesky_unavailable_message("Known-V BLUP covariance"), fixed = TRUE)
})

test_that("selection deletion does not silently invert a lower-rank approximation", {

  covariance <- diag(3L)
  covariance[1:2, 1:2] <- .covariance_underflow_fixture()
  expect_error(.selection_deleted_gaussian_block(covariance, c(0, 0, 0), 3L),
               .covariance_cholesky_unavailable_message("Retained selection-deletion covariance"), fixed = TRUE)

  # Genuine semidefinite support still uses the existing valid pseudoinverse.
  covariance <- diag(c(1, 0, 2))
  expect_null(.covariance_cholesky(.covariance_factorization(covariance)))
  result <- .selection_deleted_gaussian_block(covariance, c(.2, 0, 0), 3L)
  expect_identical(result[["mean"]], 0)
  expect_equal(result[["covariance"]], matrix(2, 1L), tolerance = 1e-15)
})

test_that("block-MVN fitting preflight checks the actual fixed integrated covariance", {

  covariance <- .covariance_underflow_fixture()
  fit <- function(sd = NULL) {
    if (is.null(sd)) {
      return(brma.mv(yi = c(0, 0), V = covariance,
        known_v_parameterization = "block_mvn", measure = "GEN",
        prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE))
    }
    brma.mv(yi = c(0, 0), V = covariance, random = ~ 1 | esid,
      data = data.frame(esid = 1:2), known_v_parameterization = "block_mvn",
      prior_heterogeneity = BayesTools::prior_random(sd = BayesTools::prior("point", list(sd))),
      measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  }
  message <- .covariance_cholesky_unavailable_message(
    "Known-V block-MVN covariance at retained rows {1, 2}")
  expect_error(fit(), message, fixed = TRUE)
  expect_error(fit(0), message, fixed = TRUE)
  rescued <- fit(.2)
  expect_s3_class(rescued, "brma.mv")
  expect_equal(.brma_mv_fixed_integrated_variance(rescued, 2L), rep(.2^2, 2L), tolerance = 0)
  expect_true(.known_v_is_numerically_positive_definite(covariance + diag(.2^2, 2L)))
})

test_that("underflowed directions cannot become whitening or latent structural zeros", {

  covariance <- .covariance_underflow_fixture()
  expect_error(.known_v_prepare(covariance, c(TRUE, TRUE), "whitened", warn_singular = FALSE),
    paste0("Known-V whitening covariance spectral variances are unavailable at working precision: ",
           "a retained positive direction has a zero or non-finite variance. ",
           "Rescale the outcome, covariance, and associated priors consistently before fitting."), fixed = TRUE)
  expect_error(.known_v_prepare(covariance, c(TRUE, TRUE), "latent", warn_singular = FALSE),
    paste0("Known-V latent decomposition is unavailable because a strictly ",
           "positive residual variance underflowed to zero. Rescale the outcome, ",
           "covariance, and associated priors consistently before fitting."), fixed = TRUE)

  extended <- matrix(0, 3L, 3L)
  extended[1:2, 1:2] <- covariance
  factorization <- .covariance_factorization(extended)
  expect_identical(factorization[["status"]], "positive_semidefinite")
  expect_error(.covariance_spectral_values(factorization),
               "a retained positive direction has a zero or non-finite variance", fixed = TRUE)
  expect_true(.known_v_nullspace_is_regularized(extended, c(FALSE, FALSE, TRUE)))
  expect_false(.known_v_nullspace_is_regularized(extended, c(FALSE, TRUE, FALSE)))
  ordinary <- .covariance_factorization(diag(c(1, 0, 2)))
  expect_identical(.covariance_spectral_values(ordinary), ordinary[["spectral_values"]])
})
