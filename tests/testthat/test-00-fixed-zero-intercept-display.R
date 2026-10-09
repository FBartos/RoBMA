context("Fixed-zero meta-regression intercept display")


.fixed_zero_intercept_test_data <- data.frame(
  yi    = c(-0.20, 0.05, 0.18, 0.31, 0.42),
  sei   = c(0.20, 0.18, 0.22, 0.19, 0.21),
  mod   = c(-1, -0.5, 0, 0.5, 1),
  study = factor(c(1, 1, 2, 2, 3))
)

.fixed_zero_intercept_point <- BayesTools::prior(
  "spike",
  parameters = list(0)
)

.fixed_zero_intercept_object <- function(mods, prior_effect = NULL) {

  args <- list(
    yi                        = .fixed_zero_intercept_test_data[["yi"]],
    sei                       = .fixed_zero_intercept_test_data[["sei"]],
    mods                      = mods,
    data                      = .fixed_zero_intercept_test_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  if (!is.null(prior_effect)) {
    args[["prior_effect"]] <- prior_effect
  }

  do.call(brma.norm, args)
}


test_that("only fixed-zero intercepts accompanying moderators are omitted", {

  no_intercept <- .fixed_zero_intercept_object(~ 0 + mod)
  fixed_zero   <- .fixed_zero_intercept_object(
    ~ mod,
    prior_effect = .fixed_zero_intercept_point
  )
  estimated    <- .fixed_zero_intercept_object(~ mod)
  intercept_only <- brma.norm(
    yi                        = yi,
    sei                       = sei,
    prior_effect              = .fixed_zero_intercept_point,
    data                      = .fixed_zero_intercept_test_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  random_model <- brma.mv(
    yi                        = yi,
    V                         = sei^2,
    mods                      = ~ 0 + mod,
    random                    = ~ (1 | study),
    data                      = .fixed_zero_intercept_test_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_true(.location_omit_fixed_zero_intercept(no_intercept))
  expect_true(.location_omit_fixed_zero_intercept(fixed_zero))
  expect_true(.location_omit_fixed_zero_intercept(random_model))
  expect_false(.location_omit_fixed_zero_intercept(estimated))
  expect_false(.location_omit_fixed_zero_intercept(intercept_only))
})


test_that("public parameter discovery follows fixed-zero intercept display", {

  fixed_zero <- .fixed_zero_intercept_object(
    ~ mod,
    prior_effect = .fixed_zero_intercept_point
  )
  estimated <- .fixed_zero_intercept_object(~ mod)

  fixed_catalog     <- .brma_parameter_catalog_unfitted(fixed_zero)
  estimated_catalog <- .brma_parameter_catalog_unfitted(estimated)

  expect_false("mu_intercept" %in% fixed_catalog[["parameter"]])
  expect_true("mu_mod" %in% fixed_catalog[["parameter"]])
  expect_true("mu_intercept" %in% estimated_catalog[["parameter"]])
})


test_that("marginal means drop the intercept from structural metadata", {

  fixed_zero <- .fixed_zero_intercept_object(
    ~ mod,
    prior_effect = .fixed_zero_intercept_point
  )
  estimated <- .fixed_zero_intercept_object(~ mod)
  formula   <- ~ mod
  terms     <- c("intercept", "mod")
  parameters <- c("mu_intercept", "mu_mod")

  fixed_setup <- .marginal_means_drop_fixed_zero_intercept(
    object     = fixed_zero,
    formula    = formula,
    terms      = terms,
    parameters = parameters
  )
  estimated_setup <- .marginal_means_drop_fixed_zero_intercept(
    object     = estimated,
    formula    = formula,
    terms      = terms,
    parameters = parameters
  )

  expect_identical(fixed_setup[["terms"]], "mod")
  expect_identical(fixed_setup[["parameters"]], "mu_mod")
  expect_identical(attr(stats::terms(fixed_setup[["formula"]]), "intercept"), 0L)
  expect_identical(estimated_setup[["terms"]], terms)
  expect_identical(estimated_setup[["parameters"]], parameters)
})


test_that("summary routes omit the fixed-zero location intercept", {

  object <- .fixed_zero_intercept_object(
    ~ mod,
    prior_effect = .fixed_zero_intercept_point
  )
  object[["fit"]] <- structure(
    list(),
    class      = "BayesTools_fit",
    prior_list = list()
  )
  observed <- NULL

  testthat::local_mocked_bindings(
    .location_omit_fixed_zero_intercept = function(object, ...) TRUE,
    .summary_estimates_pair = function(
        enabled, object, probs, include_mcmc_diagnostics, is_robma,
        conditional, main_args, conditional_args = main_args) {

      if (identical(main_args[["title"]], "Meta-Regression")) {
        observed <<- list(
          main        = main_args,
          conditional = conditional_args
        )
      }
      list(estimates = list(), conditional = list())
    },
    .package = "RoBMA"
  )

  summary.brma(object, include_mcmc_diagnostics = FALSE)

  expect_identical(observed[["main"]][["remove_parameters"]], "mu_intercept")
  expect_identical(
    observed[["conditional"]][["remove_parameters"]],
    "mu_intercept"
  )
})


test_that("stored summaries omit the fixed-zero location intercept", {

  observed <- NULL
  testthat::local_mocked_bindings(
    .location_omit_fixed_zero_intercept = function(object, ...) TRUE,
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    JAGS_estimates_table = function(...) {

      observed <<- list(...)
      data.frame(Mean = 0, row.names = "mu_mod")
    },
    .package = "BayesTools"
  )

  .object_summary(list(fit = structure(list(), class = "BayesTools_fit")))

  expect_true("mu_intercept" %in% observed[["remove_parameters"]])
  expect_true(all(c("theta", "gamma", "sampling_z", "pi", "phi") %in%
                    observed[["remove_parameters"]]))
})


test_that("intercept-only multivariate location coefficients use mu", {

  no_mods <- brma.mv(
    yi                        = yi,
    V                         = sei^2,
    random                    = ~ (1 | study),
    data                      = .fixed_zero_intercept_test_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  with_mods <- brma.mv(
    yi                        = yi,
    V                         = sei^2,
    mods                      = ~ mod,
    random                    = ~ (1 | study),
    data                      = .fixed_zero_intercept_test_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  no_mods_table <- data.frame(
    Mean      = 0,
    row.names = "intercept"
  )
  with_mods_table <- data.frame(
    Mean      = c(0, 1),
    row.names = c("intercept", "mod")
  )

  expect_identical(
    rownames(.summary_location_repair_row_labels(no_mods_table, no_mods)),
    "mu"
  )
  expect_identical(
    rownames(.summary_location_repair_row_labels(with_mods_table, with_mods)),
    c("intercept", "mod")
  )

  no_mods[["coefficients"]] <- c("(mu) intercept" = 0, "(mu) sd" = 1)
  expect_identical(names(coef(no_mods)), c("mu", "(mu) sd"))
})


# A03: the fixed zero is the intercept of the standardized coordinates. The
# original-scale intercept of a scaled predictor with center m and scale s is
# -b * m / s.
.fixed_zero_scaled_object <- function(x, standardize = TRUE, mods = ~ 0 + x) {

  brma(
    yi = yi, sei = sei, mods = mods,
    data = data.frame(yi = seq(-.1, .2, length.out = 6), sei = .1, x = x),
    measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("point", list(location = 0)),
    prior_mods = list(x = BayesTools::prior("normal", list(mean = 0, sd = 1))),
    prior_heterogeneity = BayesTools::prior("point", list(location = .1)),
    standardize_continuous_predictors = standardize,
    only_priors = TRUE, silent = TRUE
  )
}

test_that("the fixed-zero intercept gate follows the display scale", {

  off_center <- .fixed_zero_scaled_object(1:6)
  centered   <- .fixed_zero_scaled_object(c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5))
  unscaled   <- .fixed_zero_scaled_object(1:6, standardize = FALSE)

  # a scaled predictor with a nonzero center hides a nonzero original-scale
  # intercept; the standardized display has the fixed zero
  expect_false(.location_omit_fixed_zero_intercept(off_center))
  expect_false(.location_omit_fixed_zero_intercept(off_center, transform_scaled = TRUE))
  expect_true(.location_omit_fixed_zero_intercept(off_center, transform_scaled = FALSE))
  # centers exactly at zero, or no scaled predictor, keep the omission
  for (object in list(centered, unscaled)) {
    expect_true(.location_omit_fixed_zero_intercept(object))
    expect_true(.location_omit_fixed_zero_intercept(object, transform_scaled = FALSE))
  }
})

test_that("summaries show the original-scale intercept of a fixed-zero standardized model", {

  n <- 100L
  object <- .fixed_zero_scaled_object(1:6)
  setup  <- review_compile_object(object)
  slope  <- stats::qnorm((seq_len(n) - .5) / n) + .3
  samples <- cbind(slope, 0, .1)
  colnames(samples) <- c(
    review_coefficient_name(setup, "x"),
    review_coefficient_name(setup, "intercept"), "tau"
  )
  object <- review_attach_replay(setup, samples)
  object[["summary"]]      <- .object_summary(object)
  object[["coefficients"]] <- .object_coefficients(object)
  expected <- -mean(slope) * mean(1:6) / stats::sd(1:6)

  # stored summary and coefficients (original scale)
  expect_true("(mu) intercept" %in% rownames(object[["summary"]]))
  expect_equal(unname(object[["coefficients"]][["(mu) intercept"]]), expected,
               tolerance = 1e-10)
  expect_equal(unname(object[["coefficients"]][["(mu) x"]]),
               mean(slope) / stats::sd(1:6), tolerance = 1e-10)
  # public summary: original scale by default, the fixed zero when standardized
  original     <- summary(object, include_mcmc_diagnostics = FALSE)
  standardized <- summary(object, standardized_coefficients = TRUE,
                          include_mcmc_diagnostics = FALSE)
  expect_equal(original[["estimates_mods"]][["intercept", "Mean"]], expected,
               tolerance = 1e-10)
  expect_identical(rownames(standardized[["estimates_mods"]]), "x")
  expect_equal(standardized[["estimates_mods"]][["x", "Mean"]], mean(slope),
               tolerance = 1e-10)
})

test_that("marginal means of a fixed-zero standardized model stay on the fitted coordinates", {

  n <- 100L
  object <- .fixed_zero_scaled_object(1:6)
  setup  <- review_compile_object(object)
  slope  <- stats::qnorm((seq_len(n) - .5) / n) + .3
  samples <- cbind(slope, 0, .1)
  colnames(samples) <- c(
    review_coefficient_name(setup, "x"),
    review_coefficient_name(setup, "intercept"), "tau"
  )
  object <- review_attach_replay(setup, samples)

  # off-center predictor, nonzero original-scale intercept: the marginal means
  # are those of the structural-zero basis, without an intercept row
  means <- marginal_means(object)
  table <- `class<-`(summary(means), "data.frame")

  expect_identical(means[["parameters"]], "mu_x")
  expect_identical(attr(stats::terms(means[["formula"]]), "intercept"), 0L)
  expect_identical(rownames(table), c("x[-1SD]", "x[0SD]", "x[1SD]"))
  expect_equal(unname(table[, "Mean"]), c(-1, 0, 1) * mean(slope), tolerance = 1e-12)
  expect_equal(unname(table[, "SD"]), c(1, 0, 1) * stats::sd(slope), tolerance = 1e-12)
})

test_that("marginal means keep dropping the structural fixed-zero intercept", {

  off_center <- .fixed_zero_scaled_object(1:6, mods = ~ x)
  setup <- .marginal_means_drop_fixed_zero_intercept(
    object     = off_center,
    formula    = ~ x,
    terms      = c("intercept", "x"),
    parameters = c("mu_intercept", "mu_x")
  )

  expect_identical(setup[["terms"]], "x")
  expect_identical(setup[["parameters"]], "mu_x")
  expect_identical(attr(stats::terms(setup[["formula"]]), "intercept"), 0L)
})

test_that("parameter discovery lists the original-scale fixed-zero intercept", {

  off_center <- .fixed_zero_scaled_object(1:6, mods = ~ x)
  centered   <- .fixed_zero_scaled_object(c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5),
                                          mods = ~ x)

  expect_true("mu_intercept" %in% .brma_parameter_catalog_unfitted(off_center)[["parameter"]])
  expect_false("mu_intercept" %in% .brma_parameter_catalog_unfitted(centered)[["parameter"]])
})

test_that("catalogs and selectors offer the original-scale intercept of '~ 0 + x'", {

  centered <- .fixed_zero_scaled_object(c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5))
  unscaled <- .fixed_zero_scaled_object(1:6, standardize = FALSE)
  aliases  <- function(object) .brma_parameter_catalog(object)[["alias"]]

  # the formula has no intercept, but the standardized fit has the original-scale
  # intercept -b * 3.5 / sd(1:6) of its off-center predictor
  off_center <- .fixed_zero_scaled_object(1:6)
  expect_true("intercept" %in% aliases(off_center))
  expect_identical(
    .brma_parameter_select_entry(off_center, "intercept", component = "mods")[["parameter"]],
    "mu_intercept"
  )
  for (object in list(centered, unscaled)) {
    expect_false("intercept" %in% aliases(object))
    expect_error(.brma_parameter_select_entry(object, "intercept", component = "mods"))
  }

  n <- 100L
  setup <- review_compile_object(off_center)
  samples <- cbind(stats::qnorm((seq_len(n) - .5) / n) + .3, 0, .1)
  colnames(samples) <- c(
    review_coefficient_name(setup, "x"),
    review_coefficient_name(setup, "intercept"), "tau"
  )
  fitted <- review_attach_replay(setup, samples)
  expect_true("intercept" %in% aliases(fitted))
  expect_identical(
    .brma_parameter_select_entry(fitted, "intercept", component = "mods")[["parameter"]],
    "mu_intercept"
  )
  expect_true("mu_intercept" %in% hypothesis_quantities(fitted)[["parameter"]])
})
