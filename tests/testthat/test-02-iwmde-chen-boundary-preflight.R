testthat::test_that("finite declared omega conditioners retain endpoints and raw interior values", {
  context <- list(flat_prior_list = list(bias = BayesTools::prior_weightfunction(
    side = "one-sided", steps = .05, weights = BayesTools::wf_cumulative(c(1, 1)))),
    selection_spec = list(jags_omega = "omega"))
  values <- c(0, .2, .75, 1)
  transformed <- .iwmde_chen_transform_conditioning_column(
    context, values, rev(values), "omega[2]")
  testthat::expect_identical(transformed$fit, values)
  testthat::expect_identical(transformed$eval, rev(values))
  invalid <- .iwmde_chen_transform_conditioning_column(
    context, c(-.01, 1.01, NA, NaN, Inf), c(-Inf, NA), "omega[2]")
  testthat::expect_true(all(!is.finite(invalid$fit)))
  testthat::expect_true(all(!is.finite(invalid$eval)))
})

testthat::test_that("boundary conditioners use conditional normal weights without a forced fallback", {
  i <- seq_len(40L)
  omega <- rep(c(0, .2, .6, 1), 10L)
  samples <- cbind(mu = .7 * omega + sin(i) / 5, `omega[1]` = 1, `omega[2]` = omega)
  context <- list(posterior_samples = samples, indicator_names = character(),
    selection_spec = list(jags_omega = "omega"), priors = list(), formula_inputs = list(),
    flat_prior_list = list(mu = BayesTools::prior("normal", list(0, 1)),
      bias = BayesTools::prior_weightfunction("one-sided", .05, BayesTools::wf_cumulative(c(1, 1)))))
  context <- .iwmde_context_ensure_caches(context)
  spec <- list(type = "primitive", parameter = "mu")
  weight <- .iwmde_chen_log_weight(context, "mu", spec, i, samples[, "mu"], i,
    samples[, "mu"], c(-Inf, Inf))
  testthat::expect_identical(weight$method, "chen_conditional_normal")
  testthat::expect_identical(weight$fallback_count, 0L)
  conditioning <- tryCatch(.iwmde_chen_conditioning_matrix(context, "mu", spec,
    c(1L, 4L), 2:39), error = identity)
  testthat::expect_false(inherits(conditioning, "error"))
  if (!inherits(conditioning, "error")) {
    testthat::expect_identical(conditioning$columns, "omega[2]")
    testthat::expect_equal(as.numeric(conditioning$eval), omega[c(1L, 4L)])
  }
  constant <- context
  constant$posterior_samples[, "omega[2]"] <- 1
  constant <- .iwmde_context_ensure_caches(constant, reset = TRUE)
  control <- .iwmde_chen_conditioning_matrix(constant, "mu", spec, i, i)
  testthat::expect_length(control$columns, 0L)
})

testthat::test_that("omega change preserves lower-unbounded and nonomega transforms", {
  context <- list(flat_prior_list = list(
    bias = BayesTools::prior_weightfunction("one-sided", .05,
      BayesTools::wf_independent(BayesTools::prior("gamma", list(shape = 2, rate = 1)))),
    rho = BayesTools::prior("uniform", list(0, 1))), selection_spec = list(jags_omega = "omega"))
  a <- .iwmde_chen_transform_conditioning_column(context, c(0, .2, 1), 0, "omega[2]")
  testthat::expect_equal(a$fit, log1p(c(0, .2, 1)))
  testthat::expect_identical(a$eval, 0)
  a <- .iwmde_chen_transform_conditioning_column(context, c(.2, .8), .5, "rho")
  testthat::expect_equal(a$fit, qlogis(c(.2, .8)))
  testthat::expect_identical(a$eval, 0)
})

test_that("Chen boundary preflight retains constant-column semantics", {

  samples <- cbind(
    mu         = c(-.4, -.1, .2, .5),
    PET        = c(-1, -.2, .4, 1.1),
    rho        = rep(0, 4L),
    `omega[1]` = rep(1, 4L)
  )
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    selection_spec    = list(jags_omega = "omega"),
    priors            = list(),
    flat_prior_list   = list(mu = TRUE, PET = TRUE, rho = TRUE),
    formula_inputs    = list()
  )

  conditioning <- .iwmde_chen_conditioning_matrix(
    context        = context,
    parameter      = "mu",
    parameter_spec = list(type = "primitive", parameter = "mu"),
    active_rows    = 1:2,
    weight_rows    = 3:4
  )

  expect_identical(conditioning[["columns"]], "PET")
  expect_true(all(is.finite(conditioning[["fit"]])))
  expect_true(all(is.finite(conditioning[["eval"]])))
})


test_that("Chen boundary preflight does not reject evaluation-only boundaries", {

  samples <- cbind(
    mu         = c(-.4, -.1, .2, .5),
    PET        = c(-1, -.2, .4, 1.1),
    `omega[1]` = c(.2, .4, .6, 1)
  )
  context <- list(
    posterior_samples = samples,
    indicator_names   = character(),
    selection_spec    = list(jags_omega = "omega"),
    priors            = list(),
    flat_prior_list   = list(
      mu   = TRUE,
      PET  = TRUE,
      bias = BayesTools::prior_weightfunction(
        side    = "one-sided",
        steps   = .05,
        weights = BayesTools::wf_cumulative(c(1, 1))
      )
    ),
    formula_inputs    = list()
  )

  conditioning <- .iwmde_chen_conditioning_matrix(
    context        = context,
    parameter      = "mu",
    parameter_spec = list(type = "primitive", parameter = "mu"),
    active_rows    = 4L,
    weight_rows    = 1:3
  )

  expect_true(all(is.finite(conditioning[["fit"]])))
  expect_true(all(is.finite(conditioning[["eval"]])))
})



testthat::test_that("scalar and indexed marginalized local coordinates are excluded", {
  local <- c("gamma", "theta", "sampling_z", "pi", "phi", "gamma[1]", "theta[2]", "sampling_z[1]", "pi[3]", "phi[4]")
  testthat::expect_true(all(vapply(local, .iwmde_parameter_is_local_latent, logical(1))))
  other <- c("gamma_global", "theta2", "sampling_z_extra", "pilot", "phi_global", "mu")
  testthat::expect_false(any(vapply(other, .iwmde_parameter_is_local_latent, logical(1))))
  samples <- matrix(seq_len(48L), 3L, dimnames = list(NULL, c(local, other)))
  dropped <- .iwmde_drop_local_latent_sample_columns(samples)
  testthat::expect_identical(colnames(dropped), other)
})

test_that("declared formula coefficients remain global conditioning columns", {

  x <- seq(-1, 1, length.out = 40L)
  compiled <- BayesTools::JAGS_formula(~ theta, "mu", data.frame(theta = x),
    list(intercept = BayesTools::prior("point", list(0)),
      theta = BayesTools::prior("normal", list(0, 1))))
  expect_true("mu_theta" %in% names(compiled[["prior_list"]]))
  samples <- cbind(focal = sin(seq_along(x)), mu_theta = x,
    sampling_z = x + .1)
  context <- .iwmde_context_ensure_caches(list(posterior_samples = samples,
    flat_prior_list = c(compiled[["prior_list"]],
      list(focal = BayesTools::prior("normal", list(0, 1)),
        sampling_z = BayesTools::prior("normal", list(0, 1)))),
    indicator_names = character()))
  columns <- .iwmde_chen_conditioning_columns(context, "focal", list(type = "primitive"))
  expect_identical(columns, "mu_theta")
  expect_identical(colnames(.iwmde_drop_local_latent_sample_columns(samples)), c("focal", "mu_theta"))
})
