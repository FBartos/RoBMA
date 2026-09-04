test_that("selection reweighting metrics retain their exact definitions", {

  constant <- .selection_approximation_weight_metrics(rep(-1000, 32L))
  expect_equal(constant[["ess_fraction"]], 1)
  expect_equal(constant[["ess_fraction_mcse"]], 0)
  expect_equal(constant[["total_variation"]], 0)
  expect_equal(constant[["total_variation_mcse"]], 0)
  expect_equal(constant[["log_weight_iqr"]], 0)

  log_weights <- log(c(0.2, 0.5, 1, 0.7, 0.1, 0.9))
  shifted     <- .selection_approximation_weight_metrics(log_weights - 800)
  direct      <- exp(log_weights)
  expected_ess <- mean(direct)^2 / mean(direct^2)
  expected_tv  <- 0.5 * mean(abs(direct / mean(direct) - 1))

  expect_equal(shifted[["ess_fraction"]], expected_ess)
  expect_equal(shifted[["total_variation"]], expected_tv)
  expect_equal(
    shifted[["log_weight_iqr"]],
    unname(stats::IQR(log_weights, type = 8))
  )
  expect_true(shifted[["ess_fraction_mcse"]] > 0)
  expect_true(shifted[["total_variation_mcse"]] > 0)
})


test_that("selection latent draws preserve every covariance factor type", {

  latent_samples <- 100000L
  K              <- 4L
  X <- matrix(
    c(1, 0.2, 1, -0.4, 0, 0.6, 0.5, 0),
    nrow = K,
    ncol = 2L,
    byrow = TRUE
  )
  coefficient_factor <- matrix(c(1, 0.2, 0, 0.7), 2L, 2L)
  group_map           <- c(1L, 1L, 2L, 2L)
  basis               <- X %*% coefficient_factor

  plans <- list(
    group = list(
      type                  = "group",
      model_matrix          = X,
      group_map             = group_map,
      coefficient_structure = "dense"
    ),
    row_group = list(
      type                  = "row_group",
      model_matrix          = X,
      group_map             = group_map,
      coefficient_structure = "dense"
    ),
    known_group = list(
      type                  = "known_group",
      model_matrix          = X,
      group_map             = group_map,
      coefficient_structure = "dense",
      group_covariance      = matrix(c(1, 0.3, 0.3, 1.4), 2L, 2L)
    )
  )
  states <- list(
    group = list(coefficient_factor = coefficient_factor),
    row_group = list(
      coefficient_factor = coefficient_factor,
      row_scale           = c(0.5, 1, 1.5, 2)
    ),
    known_group = list(coefficient_factor = coefficient_factor)
  )
  expected <- list(
    group = tcrossprod(basis) * outer(group_map, group_map, "=="),
    row_group = {
      scaled_basis <- basis * states[["row_group"]][["row_scale"]]
      tcrossprod(scaled_basis) * outer(group_map, group_map, "==")
    },
    known_group = tcrossprod(basis) *
      plans[["known_group"]][["group_covariance"]][
        group_map,
        group_map,
        drop = FALSE
      ]
  )

  for (type in names(plans)) {
    set.seed(440 + match(type, names(plans)))
    draws <- .selection_approximation_factor_draws(
      plan           = plans[[type]],
      state          = states[[type]],
      rows           = seq_len(K),
      latent_samples = latent_samples,
      K              = K
    )
    expect_lt(max(abs(stats::cov(draws) - expected[[type]])), 0.025)
  }

  dense_covariance <- tcrossprod(matrix(
    c(1, 0.2, 0, 0.5, 1, 0.4, 0, 0.4),
    nrow = K,
    ncol = 2L
  ))
  set.seed(445)
  dense_draws <- .selection_approximation_factor_draws(
    plan           = list(type = "dense"),
    state          = list(covariance = dense_covariance),
    rows           = seq_len(K),
    latent_samples = latent_samples,
    K              = K
  )
  expect_lt(max(abs(stats::cov(dense_draws) - dense_covariance)), 0.025)
})


test_that("selection approximation summaries expose posterior and MC uncertainty", {

  metrics <- list(
    ess_fraction = matrix(c(0.2, 0.4, 0.6, 0.8), 2L, 2L),
    ess_fraction_mcse = matrix(c(0.01, 0.02, 0.03, 0.04), 2L, 2L),
    total_variation = matrix(c(0.1, 0.2, 0.3, 0.4), 2L, 2L),
    total_variation_mcse = matrix(c(0.02, 0.03, 0.04, 0.05), 2L, 2L),
    log_weight_iqr = matrix(c(1, 2, 3, 4), 2L, 2L)
  )
  out <- .selection_approximation_summarize(
    metrics,
    list(1:2, 3:5)
  )

  expect_s3_class(out, "data.frame")
  expect_identical(out[["rows"]], c("1,2", "3,4,5"))
  expect_identical(out[["n_estimates"]], c(2L, 3L))
  expect_equal(out[["ess_fraction_median"]], c(0.3, 0.7))
  expect_equal(out[["ess_fraction_mcse_max"]], c(0.02, 0.04))
  expect_equal(out[["total_variation_mcse_max"]], c(0.03, 0.05))
})


test_that("selection approximation notifications use posterior-median TV", {

  expect_silent(.selection_approximation_notify(c(0.01, 0.05)))
  expect_message(
    .selection_approximation_notify(c(0.01, 0.064)),
    paste0(
      "Approximate selection-likelihood diagnostic: the largest block-level ",
      "posterior-median total variation distance was 6.4%. Inspect ",
      "'selection_approximation_diagnostics()' for block-level results."
    ),
    fixed = TRUE
  )
  expect_warning(
    .selection_approximation_notify(c(0.01, 0.10)),
    paste0(
      "The approximate selection likelihood showed substantial ",
      "latent-distribution reweighting: the largest block-level ",
      "posterior-median total variation distance was 10.0%. The approximate ",
      "and exact likelihoods may yield different inference. Consider ",
      "refitting with 'selection_likelihood' set to 'exact'."
    ),
    fixed = TRUE
  )
})
