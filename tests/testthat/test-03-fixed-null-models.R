context("Downstream methods for fully fixed models")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_no_fits()

fixed_model_names <- c(
  "fixed_null_brma",
  "fixed_nonzero_brma",
  "fixed_nonzero_brma_mv",
  "fixed_null_brma_glmm",
  "fixed_null_brma_glmm_pois",
  "fixed_null_bPET",
  "fixed_null_bPEESE",
  "fixed_null_bselmodel"
)

test_that("fixed mu and tau remain available to downstream methods", {

  skip_if_missing_fits(fixed_model_names)

  for (name in fixed_model_names) {
    fit      <- load_fit(name, validate = FALSE)
    info     <- load_info(name, validate = FALSE)
    draws    <- as.matrix(fit[["fit"]][["mcmc"]])
    log_lik  <- log_lik(fit)
    loo_fit  <- loo(fit)
    fit      <- suppressWarnings(add_waic(fit))
    waic_fit <- waic(fit)

    expect_true(all(c("mu", "tau") %in% colnames(draws)), info = name)
    expect_equal(unique(draws[, "mu"]), info[["mu"]], info = name)
    expect_equal(unique(draws[, "tau"]), info[["tau"]], info = name)
    expect_true(all(c("mu", "tau") %in% rownames(fit[["summary"]])), info = name)
    expect_equal(fit[["summary"]]["mu", "Mean"], info[["mu"]], info = name)
    expect_equal(fit[["summary"]]["tau", "Mean"], info[["tau"]], info = name)
    expect_s3_class(summary(fit), "summary.brma")
    expect_s3_class(predict(fit, type = "terms", quiet = TRUE), "brma_samples")
    expect_true(is.finite(logml(fit)), info = name)
    expect_true(all(apply(log_lik, 2L, stats::var) == 0), info = name)
    expect_equal(
      loo_fit[["estimates"]]["p_loo", "Estimate"],
      0,
      tolerance = 1e-12,
      info      = name
    )
    expect_equal(loo::pareto_k_values(loo_fit), rep(0, ncol(log_lik)), info = name)
    expect_equal(
      waic_fit[["estimates"]]["p_waic", "Estimate"],
      0,
      tolerance = 1e-12,
      info      = name
    )
  }
})

test_that("zero-dimensional marginal likelihoods are exact", {

  exact_names <- setdiff(
    fixed_model_names,
    c("fixed_null_brma_glmm", "fixed_null_brma_glmm_pois")
  )
  skip_if_missing_fits(exact_names)

  for (name in exact_names) {
    fit     <- load_fit(name, validate = FALSE)
    marglik <- fit[["marglik"]]
    error   <- tryCatch(bridge_sampler(fit), error = identity)

    expect_s3_class(error, "RoBMA_exact_marglik_no_bridge", info = name)
    expect_match(conditionMessage(error), "evaluated exactly", info = name)
    expect_identical(
      marglik[["aggregation"]][["rule"]],
      "exact_zero_dimensional",
      info = name
    )
    expect_identical(marglik[["aggregation"]][["n_repetitions"]], 0L, info = name)
    expect_identical(error[["reason"]], "exact_zero_dimensional", info = name)
    expect_equal(error[["logml"]], logml(fit), info = name)
  }

  for (name in c(
    "fixed_null_brma",
    "fixed_nonzero_brma",
    "fixed_null_bPET",
    "fixed_null_bPEESE",
    "fixed_null_bselmodel"
  )) {
    fit <- load_fit(name, validate = FALSE)
    expect_equal(logml(fit), sum(log_lik(fit)[1L, ]), tolerance = 1e-12, info = name)
  }

  fit_mv  <- load_fit("fixed_nonzero_brma_mv", validate = FALSE)
  info_mv <- load_info("fixed_nonzero_brma_mv", validate = FALSE)
  expected_mv <- mvtnorm::dmvnorm(
    info_mv[["yi"]],
    mean  = rep(info_mv[["mu"]], length(info_mv[["yi"]])),
    sigma = info_mv[["V"]] + diag(info_mv[["tau"]]^2, length(info_mv[["yi"]])),
    log   = TRUE
  )
  expect_equal(logml(fit_mv), expected_mv, tolerance = 1e-12)
})

test_that("fixed-model marginal likelihoods support hypothesis testing", {

  model_names <- c("fixed_null_brma", "fixed_nonzero_brma")
  skip_if_missing_fits(model_names)

  null        <- load_fit(model_names[[1L]], validate = FALSE)
  alternative <- load_fit(model_names[[2L]], validate = FALSE)
  log_bf      <- bf(null, alternative, log = TRUE)
  bf_result   <- bf(null, alternative)
  probability <- post_prob(
    null,
    alternative,
    model_names = c("null", "alternative")
  )

  expect_s3_class(log_bf, "bf_default")
  expect_s3_class(bf_result, "bf_default")
  expect_equal(log_bf[["bf"]], logml(null) - logml(alternative), tolerance = 1e-12)
  expect_equal(bf_result[["bf"]], exp(log_bf[["bf"]]), tolerance = 1e-12)
  expect_equal(sum(probability), 1, tolerance = sqrt(.Machine$double.eps))
  expect_equal(names(probability), c("null", "alternative"))
})
