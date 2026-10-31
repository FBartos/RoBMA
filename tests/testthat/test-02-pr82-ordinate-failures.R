source(testthat::test_path("common-functions.R"))

test_that("public point underflow retains log evidence and diagnostic rows", {

  name <- "brma.mv_block_mvn_fixed_random_null"
  skip_if_missing_fits(name)
  fit <- load_fit(name)
  input <- load_info(name)
  prior <- fit[["priors"]][["outcome"]][["mu"]]
  V <- input[["V"]]
  y <- input[["data"]][["yi"]]
  variance <- 1 / (1 / prior[["parameters"]][["sd"]]^2 + sum(solve(V, rep(1, length(y)))))
  mean <- variance * (prior[["parameters"]][["mean"]] / prior[["parameters"]][["sd"]]^2 + sum(solve(V, y)))
  exact_log_ordinate <- stats::dnorm(10, mean, sqrt(variance), log = TRUE)
  expect_true(is.finite(exact_log_ordinate))
  expect_identical(exp(exact_log_ordinate), 0)
  for (method in c("qCMDE", "IWMDE")) {
    failure <- tryCatch(hypothesis(fit, "mu = 10", density_method = method,
      density_control = list(samples = 40L, n_points = 20L)),
      RoBMA_density_ordinate_error = identity)
    expect_s3_class(failure, "RoBMA_density_ordinate_error")
    record <- density_diagnostics(failure)
    expect_equal(nrow(record), 1L)
    expect_identical(record[["numerical_status"]], "underflow")
    expect_identical(record[["ordinate"]], 0)
    expect_true(is.finite(record[["log_ordinate"]]))
    expect_identical(record[["requested_value"]], 10)
    expect_false(record[["bf_grade_met"]])
    expect_false(record[["precision_target_met"]])
    expect_true(is.na(record[["ess"]]))
    expect_match(conditionMessage(failure), "underflowed to zero", fixed = TRUE)
    expect_false(grepl("increas.*samples", conditionMessage(failure)))
    # IWMDE additionally estimates its weight function from the fixed sample,
    # so its log estimate has sampling error and is not compared here.
    if (method == "qCMDE") {
      # Error model. Every row's conditional law is this Normal posterior, so
      # the rows share one normalizer, contribute identical terms (no Monte
      # Carlo error), and carry the same log error at every value: the reported
      # truncation bound b = t / (1 - t) with t <= 1 - 0.999 (the default
      # 'normalization_prob') adds at most log(1 + b), and the trapezoid error
      # d is bounded by the nested-grid change c, log(1 + d) in
      # [log(1 - c), -log(1 - c)]. The underflowed record keeps b but not c,
      # so c comes from a finite ordinate of the same design.
      companion <- density_diagnostics(hypothesis(fit, "mu = 0",
        density_method = method, density_control = list(samples = 40L, n_points = 20L)))
      bound  <- record[["truncation_ordinate_bound"]]
      change <- companion[["ordinate_relative_change"]]
      expect_identical(record[["normalization_truncation_status"]], "exact")
      expect_lte(record[["normalization_truncation"]], (1 - .999) * (1 + 1e-9))
      expect_equal(bound, companion[["truncation_ordinate_bound"]])
      error <- record[["log_ordinate"]] - exact_log_ordinate
      expect_lte(error, log1p(bound) - log1p(-change))
      expect_gte(error, log1p(-change))
    }
  }
  mixed <- .iwmde_estimate(.iwmde_context(fit), "mu", "qCMDE",
    list(samples = 40L, n_points = 20L), outputs = "ordinate", values = c(0, 10))
  expect_false(mixed[["sampling_design"]][["bf_grade_met"]])
  expect_false(mixed[["sampling_design"]][["target_met"]])
  accepted <- .iwmde_posterior_ordinate_keep_values(mixed[["posterior_ordinate"]], 0)
  expect_true(accepted[["diagnostics"]][["bf_grade_met"]])
  expect_identical(mixed[["ordinate_failures"]][["requested_value"]], 10)
})
