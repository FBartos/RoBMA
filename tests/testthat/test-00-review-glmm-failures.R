test_that("AGHQ recovery preserves unrelated conditions even when their text matches", {

  unrelated <- structure(list(message = "Binomial AGHQ mode failed at sample 1, observation 1.",
    call = NULL, marker = "programming"), class = c("review_programming_error", "error", "condition"))
  result <- tryCatch(.glmm_aghq_dispatch(1L, 1L, "bin",
    evaluate_aghq = function(k) stop(unrelated),
    evaluate_grid = function(k) stop("fallback must not run")), error = identity)
  expect_identical(result, unrelated)
  typed <- structure(list(message = "wording does not select recovery", outcome_type = "bin"),
    class = c("RoBMA_glmm_aghq_nonconvergence", "RoBMA_glmm_numerical_error", "error", "condition"))
  expect_true(.glmm_aghq_numerical_failure(typed, "bin"))
  expect_false(.glmm_aghq_numerical_failure(typed, "pois"))
})


test_that("point AGHQ failures report their typed family and observed location", {

  condition <- tryCatch(.glmm_point_mode(function(x) rep(1, length(x)), 1L, "bin", 7L),
    error = identity)
  expect_s3_class(condition, "RoBMA_glmm_aghq_nonconvergence")
  expect_s3_class(condition, "RoBMA_glmm_numerical_error")
  expect_identical(condition[["outcome_type"]], "bin")
  expect_identical(condition[["sample"]], 1L)
  expect_identical(condition[["observation"]], 7L)
  expect_identical(conditionMessage(condition), paste0(
    "Binomial AGHQ was rejected by diagnostics: the mode could not be bracketed ",
    "within 20 bracket expansions at sample 1, observation 7. Inspect the data ",
    "and 'prior_baserate' before refitting the model."))
  expect_identical(condition[["reason"]], "mode_bracket")
  poisson <- tryCatch(.glmm_point_mode(function(x) rep(1, length(x)), 1L, "pois", 3L),
    error = identity)
  expect_s3_class(poisson, "RoBMA_glmm_aghq_nonconvergence")
  expect_match(conditionMessage(poisson), "Inspect the data and 'prior_lograte' before",
    fixed = TRUE)
})


test_that("prior-CDF failures retain their specific class and numerical family", {

  condition <- .glmm_grid_nonconvergence_condition("bin", 7L,
    c(2e-5, 3e-6), c(4e-7, 5e-8))
  expect_s3_class(condition, "RoBMA_glmm_grid_nonconvergence")
  expect_s3_class(condition, "RoBMA_glmm_numerical_error")
  expect_true(.glmm_grid_numerical_failure(condition, "bin"))
  expect_false(.glmm_grid_numerical_failure(condition, "pois"))
  expect_identical(conditionMessage(condition), paste0(
    "Binomial prior-CDF quadrature was rejected by diagnostics at observation 7: ",
    "last theta changes 2e-05, 3e-06; last nuisance changes 4e-07, 5e-08. ",
    "Inspect the data and 'prior_baserate' before refitting the model."))
})


test_that("native AGHQ failures transport diagnostics and explicit observation mapping", {

  control <- .glmm_aghq_control(orders = c(5L, 9L, 13L),
    tolerance = 1e-15, consecutive = 2L)
  cases <- list(
    bin = list(fun = .glmm_binom_aghq, args = list(ai = 7L, ci = 3L,
      n1i = 24L, n2i = 21L, mu_samples = matrix(.35), tau_within = matrix(.7),
      weights = NULL, prior_spec = c(alpha = 1.5, beta = 2.5), control = control)),
    pois = list(fun = .glmm_pois_aghq, args = list(x1i = 7L, x2i = 11L,
      t1i = 15, t2i = 18, mu_samples = matrix(.35), tau_within = matrix(.6),
      weights = NULL, prior_spec = c(mean = -1, sd = 1.3), control = control)))
  for (outcome in names(cases)) for (observation in list(NULL, 7L)) {
    args <- cases[[outcome]][["args"]]
    if (!is.null(observation)) args[["observation"]] <- observation
    condition <- tryCatch(do.call(cases[[outcome]][["fun"]], args), error = identity)
    expect_s3_class(condition, "RoBMA_glmm_aghq_nonconvergence")
    expect_s3_class(condition, "RoBMA_glmm_numerical_error")
    expect_identical(condition[["outcome_type"]], outcome)
    expect_identical(condition[["reason"]], "quadrature_nonconvergence")
    expect_identical(condition[["sample"]], 1L)
    expect_identical(condition[["observation"]], if (is.null(observation)) 1L else observation)
    expect_true(is.finite(condition[["change"]]))
    expect_identical(condition[["order"]], 13L)
  }
})
