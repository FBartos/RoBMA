test_that("QQ limits contain one-sided residual ranges and reject missing residuals", {

  z <- c(0.5, 1, 2)
  testthat::local_mocked_bindings(
    residuals.brma = function(...) z,
    .package = "RoBMA"
  )
  get_data <- function() .qqnorm_data(
    x = list(), type = "LOO-PIT", unit = "estimate", conditioning_depth = "marginal",
    envelope = FALSE, conf_level = 95, bonferroni = FALSE, reps = 10,
    smooth = FALSE, max_samples = 10, xlim = NULL, ylim = NULL,
    xlab = NULL, ylab = NULL, dots = list()
  )
  out <- get_data()
  expect_lt(out$ylim[1L], min(z))
  expect_gt(out$ylim[2L], max(z))
  z <- -z
  out <- get_data()
  expect_lt(out$ylim[1L], min(z))
  expect_gt(out$ylim[2L], max(z))
  z[1L] <- NA_real_
  expect_error(get_data(), "finite standardized residuals", fixed = TRUE)
})

test_that("radial confidence limits survive scientific-notation probabilities", {

  values <- matrix(seq(-1, 1, length.out = 100), ncol = 1L,
                   dimnames = list(NULL, "mu"))
  testthat::local_mocked_bindings(
    pooled_effect = function(x, probs) .new_effect_brma_samples(
      samples = values, n_chains = 1L, n_iter = 100L, title = "test", probs = probs
    ),
    .outcome_data_yi = function(...) c(-0.2, 0.1, 0.4),
    .outcome_data_vi = function(...) c(0.04, 0.09, 0.16),
    .get_radial_tau_rows = function(...) rep(0.1, 3),
    .package = "RoBMA"
  )
  out <- .radial_data(
    x = list(), center = FALSE, xlim = NULL, zlim = NULL, xlab = NULL,
    zlab = NULL, atz = NULL, aty = NULL, steps = 5, level = 99.99,
    digits = 2, transf = NULL, targs = NULL, dots = .set_dots_radial(arc.res = 2)
  )
  expect_equal(out$ci_values[c(1L, 3L)],
               unname(stats::quantile(values, c(0.00005, 0.99995))))
  expect_gte(nrow(out$ci_arc), 2L)
})
