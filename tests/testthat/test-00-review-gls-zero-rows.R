.review_gls_replay <- function(V, y, x, sd = NULL, saturated = FALSE) {

  K <- length(y)
  dat <- data.frame(x = x, row = factor(seq_len(K)),
    study = factor(if (K == 3L) c("a", "a", "b") else seq_len(K)))
  args <- list(yi = y, V = V, mods = if (saturated) ~ row else ~ x,
    random = if (is.null(sd)) NULL else ~ 1 | study, data = dat,
    standardize_continuous_predictors = FALSE,
    prior_effect = prior("normal", list(0, 1)), prior_mods = prior("normal", list(0, 1)),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  if (!is.null(sd)) args$prior_heterogeneity <- prior_random(sd = prior("point", list(sd)))
  compiled <- review_compile_object(suppressWarnings(do.call(brma.mv, args)))
  priors <- c(compiled$compiled$mu$prior_list,
    .create_fit_priors(compiled$object$data, compiled$object$priors))
  names <- review_prior_columns(priors)
  samples <- matrix(0, 24L, length(names), dimnames = list(NULL, names))
  for (name in names(priors)) if (BayesTools::is.prior.point(priors[[name]])) {
    samples[, name] <- priors[[name]]$parameters$location
  }
  terms <- .formula_design_random_effects_by_mode(compiled$object$formula_design$mu, "sampled")
  for (term in terms) {
    names <- unlist(lapply(seq_len(term$n_columns), function(j) {
      paste0(term$parameter_stem, "_xRE_Zx[", seq_len(term$n_groups), ",", j, "]")
    }))
    samples <- cbind(samples, matrix(0, 24L, length(names), dimnames = list(NULL, names)))
  }
  review_attach_replay(compiled, samples)
}

test_that("known-V residual zeros use the covariance multiplier", {

  V <- matrix(c(1, .5, 0, .5, 1, 0, 0, 0, 1), 3L)
  object <- .review_gls_replay(V, c(0, 1, -1), c(1, 0, 0))
  for (depth in c("estimate", "marginal")) {
    actual <- rstandard(object, conditioning_depth = depth)
    expect_equal(as.numeric(actual[1, ]), c(.5, sqrt(.125), sqrt(2)), tolerance = 1e-12)
  }
  diagonal <- .review_gls_replay(diag(3), c(0, 1, -1), c(1, 0, 0))
  expect_identical(as.numeric(rstandard(diagonal)[1, ]), c(0, 0, 0))
  saturated <- .review_gls_replay(V, c(0, 1, -1), c(1, 0, 0), saturated = TRUE)
  expect_true(all(as.matrix(rstandard(saturated)) == 0))
})

test_that("marginal zero certification uses full random covariance beyond diagonal V", {

  object <- .review_gls_replay(diag(3), c(0, 1, -1), c(1, 0, 0), sd = sqrt(.5))
  actual <- rstandard(object, conditioning_depth = "marginal")
  expect_equal(as.numeric(actual[1, ]), c(1/3, sqrt(1/12), 2/sqrt(3)), tolerance = 1e-12)
  diagonal <- .review_gls_replay(diag(3), c(0, 1, -1), c(1, 0, 0), sd = 0)
  expect_identical(as.numeric(rstandard(diagonal)[1, ]), c(0, 0, 0))
  saturated <- .review_gls_replay(diag(3), c(0, 1, -1), c(1, 0, 0), sd = sqrt(.5), saturated = TRUE)
  expect_true(all(as.matrix(rstandard(saturated)) == 0))
})

test_that("zero certificates follow random covariance states and saturation", {

  plan <- list(sampling_covariance = diag(3),
    random_covariance_plans = list(list(type = "group", model_matrix = matrix(1, 3, 1),
      group_map = c(1L, 1L, 2L))),
    random_covariance_states = list(list(list(coefficient_factor = matrix(0))),
      list(list(coefficient_factor = matrix(.7)))))
  expect_identical(.known_v_isolated_projection_rows(plan, 1L, 1L), 1L)
  expect_identical(.known_v_isolated_projection_rows(plan, 2L, 1L), integer())
  expect_identical(.known_v_isolated_projection_rows(plan, 2L, 1:3), 1:3)
  expect_identical(.known_v_isolated_projection_rows(plan, 2L, 3L), 3L)
})
