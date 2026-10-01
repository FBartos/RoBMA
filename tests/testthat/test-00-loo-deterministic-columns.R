context("Exact LOO for deterministic columns")


test_that("deterministic log-likelihood columns are identified individually", {

  log_lik <- cbind(
    finite_constant   = rep(-1, 5),
    infinite_constant = rep(Inf, 5),
    varying           = seq_len(5)
  )

  expect_identical(
    unname(.deterministic_log_lik_columns(log_lik)),
    c(TRUE, FALSE, FALSE)
  )
})


test_that("user-supplied relative ESS is exact for deterministic columns", {

  deterministic <- c(TRUE, FALSE, TRUE)

  expect_equal(
    .loo_prepare_r_eff(0.25, 3, deterministic),
    c(1, 0.25, 1)
  )
  expect_equal(
    .loo_prepare_r_eff(c(0.25, 0.50, 0.75), 3, deterministic),
    c(1, 0.50, 1)
  )
  expect_equal(
    .loo_prepare_r_eff(rep(NA_real_, 3), 3, rep(TRUE, 3)),
    rep(1, 3)
  )
  expect_error(
    .loo_prepare_r_eff(rep(NA_real_, 3), 3, deterministic),
    "finite and positive",
    fixed = TRUE
  )
  expect_error(
    .loo_prepare_r_eff(c(0.5, 0.5), 3, deterministic),
    "one value per observation",
    fixed = TRUE
  )
})


test_that("finite varying LOO columns are invariant to extreme log offsets", {

  n_samples <- 1200L
  chain_id  <- rep(1:2, each = n_samples / 2L)
  ordinary  <- -1 + 0.01 * sin(seq_len(n_samples))
  baseline  <- 1e-6 * seq_len(n_samples)
  log_lik <- cbind(
    ordinary = ordinary,
    extreme  = baseline - 1000
  )
  ordinary_scale <- cbind(
    ordinary = ordinary,
    extreme  = baseline
  )
  centered <- .loo_center_finite_columns(log_lik)
  r_eff <- loo::relative_eff(
    exp(centered[["log_lik"]]),
    chain_id = chain_id
  )

  expect_true(all(is.finite(r_eff) & r_eff > 0))
  expect_lt(r_eff[[2L]], 0.01)
  expect_error(
    .loo_prepare_r_eff(c(r_eff[[1L]], NA_real_), 2L, c(FALSE, FALSE)),
    "finite and positive",
    fixed = TRUE
  )

  result <- .loo_with_deterministic_columns(
    log_lik       = log_lik,
    r_eff         = r_eff,
    deterministic = c(FALSE, FALSE),
    cores         = 1
  )
  reference <- .loo_with_deterministic_columns(
    log_lik       = ordinary_scale,
    r_eff         = r_eff,
    deterministic = c(FALSE, FALSE),
    cores         = 1
  )

  expect_equal(
    result[["pointwise"]]["extreme", "elpd_loo"] + 1000,
    reference[["pointwise"]]["extreme", "elpd_loo"]
  )
  expect_equal(
    result[["pointwise"]]["extreme", "looic"] - 2000,
    reference[["pointwise"]]["extreme", "looic"]
  )
  expect_equal(
    result[["pointwise"]][, c(
      "mcse_elpd_loo", "p_loo", "influence_pareto_k"
    )],
    reference[["pointwise"]][, c(
      "mcse_elpd_loo", "p_loo", "influence_pareto_k"
    )]
  )
  expect_true(is.finite(
    result[["pointwise"]]["extreme", "mcse_elpd_loo"]
  ))
  expect_equal(
    exp(stats::weights(result[["psis_object"]])),
    exp(stats::weights(reference[["psis_object"]]))
  )
})


test_that("precomputed centering gives the unchanged exact LOO result", {

  set.seed(321)
  log_lik <- cbind(
    exact     = rep(-2, 1000L),
    variable1 = stats::rnorm(1000L, -1, 0.2),
    variable2 = stats::rnorm(1000L, -3, 0.5)
  )
  deterministic <- .deterministic_log_lik_columns(log_lik)
  variable      <- !deterministic
  centered <- .loo_center_finite_columns(
    log_lik[, variable, drop = FALSE]
  )
  r_eff <- c(1, 0.8, 0.9)

  result <- .loo_with_deterministic_columns(
    log_lik           = log_lik,
    r_eff             = r_eff,
    deterministic     = deterministic,
    cores             = 1,
    centered_variable = centered
  )
  reference <- .loo_with_deterministic_columns(
    log_lik       = log_lik,
    r_eff         = r_eff,
    deterministic = deterministic,
    cores         = 1
  )

  expect_identical(result, reference)
})


test_that("mixed deterministic LOO columns have exact saved diagnostics", {

  set.seed(123)
  n_samples <- 1000L
  log_lik <- cbind(
    exact    = rep(-1, n_samples),
    variable = stats::rnorm(n_samples, mean = -1, sd = 0.5)
  )
  deterministic <- .deterministic_log_lik_columns(log_lik)

  result <- .loo_with_deterministic_columns(
    log_lik       = log_lik,
    r_eff         = c(1, 0.8),
    deterministic = deterministic,
    cores         = 1
  )
  reference <- loo::loo(
    log_lik[, "variable", drop = FALSE],
    r_eff     = 0.8,
    save_psis = TRUE,
    cores     = 1
  )

  expect_s3_class(result, "psis_loo")
  expect_equal(
    result[["pointwise"]]["exact", ],
    c(
      elpd_loo           = -1,
      mcse_elpd_loo      = 0,
      p_loo              = 0,
      looic              = 2,
      influence_pareto_k = 0
    )
  )
  expect_equal(
    result[["pointwise"]]["variable", ],
    reference[["pointwise"]]["variable", ]
  )
  expect_equal(result[["diagnostics"]][["pareto_k"]][1], 0)
  expect_equal(result[["diagnostics"]][["n_eff"]][1], n_samples)
  expect_equal(result[["diagnostics"]][["r_eff"]][1], 1)

  psis <- result[["psis_object"]]
  expect_equal(psis[["diagnostics"]][["pareto_k"]][1], 0)
  expect_equal(psis[["diagnostics"]][["n_eff"]][1], n_samples)
  expect_equal(attr(psis, "r_eff", exact = TRUE)[1], 1)
  expect_equal(unname(attr(psis, "tail_len", exact = TRUE)[1]), 0)
  expect_equal(psis[["log_weights"]][, "exact"], rep(0, n_samples))
  expect_equal(
    exp(stats::weights(psis)[, "exact"]),
    rep(1 / n_samples, n_samples)
  )
  expect_equal(
    psis[["log_weights"]][, "variable"],
    reference[["psis_object"]][["log_weights"]][, 1]
  )
})


test_that("extreme constant LOO columns remain exact and finite", {

  n_samples <- 500L
  log_lik <- matrix(
    -1000,
    nrow     = n_samples,
    ncol     = 1,
    dimnames = list(NULL, "extreme")
  )

  expect_warning(
    result <- .loo_with_deterministic_columns(
      log_lik       = log_lik,
      r_eff         = 1,
      deterministic = TRUE,
      cores         = 1
    ),
    NA
  )

  expect_equal(
    result[["pointwise"]]["extreme", ],
    c(
      elpd_loo           = -1000,
      mcse_elpd_loo      = 0,
      p_loo              = 0,
      looic              = 2000,
      influence_pareto_k = 0
    )
  )
  expect_equal(loo::mcse_loo(result), 0)
  expect_true(all(is.finite(result[["estimates"]][, "Estimate"])))
  expect_equal(
    exp(stats::weights(result[["psis_object"]])[, "extreme"]),
    rep(1 / n_samples, n_samples)
  )
})


test_that("mixed LOO preserves warnings from varying columns", {

  set.seed(123)
  n_samples <- 4000L
  log_lik <- cbind(
    exact    = rep(-1, n_samples),
    variable = -stats::rexp(n_samples, rate = 1)
  )
  warnings <- character()

  result <- withCallingHandlers(
    .loo_with_deterministic_columns(
      log_lik       = log_lik,
      r_eff         = c(1, 1),
      deterministic = c(TRUE, FALSE),
      cores         = 1
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("Pareto k", warnings, fixed = TRUE)))
  expect_false(any(grepl("tail values are the same", warnings, fixed = TRUE)))
  expect_equal(result[["diagnostics"]][["pareto_k"]][1], 0)
  expect_gt(result[["diagnostics"]][["pareto_k"]][2], 0.7)
})
