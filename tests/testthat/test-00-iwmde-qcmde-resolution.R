.resolution_seed_fit <- function(divisor, prior_effect, offset_fraction = .25, tau_broad = 1) {
  control <- .density_control_normalize("qCMDE")
  control <- .hypothesis_plan_density_control(control)
  N <- control[["normalization_points"]]
  S <- control[["samples"]]
  K <- 2L
  q <- stats::qnorm((1 - control[["normalization_prob"]]) / 2, lower.tail = FALSE)
  h <- 2 * q / (N - 1L) / divisor
  sei <- sqrt(K * (h^2 / (1 - h^2)))
  precision <- function(tau) K / (sei^2 + tau^2)
  A <- precision(c(0, tau_broad))
  conditional_sds <- 1 / sqrt(1 + A)
  mean_factors    <- A / (1 + A)
  step   <- 2 * q * conditional_sds[[2L]] / (N - 1L)
  offset <- offset_fraction * step
  ybar   <- offset / diff(rev(mean_factors))
  conditional_means <- mean_factors * ybar
  object <- BMA.norm(
    yi = rep(ybar, K), sei = rep(sei, K), measure = "GEN",
    prior_effect = prior_effect,
    prior_effect_null = NULL,
    prior_heterogeneity = BayesTools::prior("cauchy",
      list(location = 0, scale = 1), truncation = list(lower = 0, upper = Inf)),
    prior_heterogeneity_null = BayesTools::prior("point", list(location = 0)),
    only_priors = TRUE, silent = TRUE
  )
  samples <- cbind(
    mu  = c(conditional_means[[1L]] + conditional_sds[[1L]] *
              stats::qnorm((seq_len(S - 1L) - .5) / (S - 1L)),
            conditional_means[[2L]]),
    tau = c(rep(0, S - 1L), tau_broad)
  )
  priors     <- .create_fit_priors(object[["data"]], object[["priors"]])
  components <- attr(priors[["tau"]], "components")
  samples <- cbind(samples,
    tau_indicator = c(rep(which(components == "null"), S - 1L),
                      which(components == "alternative")),
    mu_indicator  = rep(1, S))
  samples <- samples[append(seq_len(S - 1L), S, after = 249L), , drop = FALSE]
  fit <- coda::mcmc.list(coda::mcmc(samples))
  attr(fit, "formula_design") <- object[["formula_design"]]
  attr(fit, "prior_list")     <- priors
  object[["fit"]] <- as_bayestools_fit(fit)
  list(object = object, control = control, S = S, K = K, ybar = ybar, sei = sei, tau_broad = tau_broad,
       value = conditional_means[[1L]], step = step, sds = conditional_sds)
}

# Reference: the mixture of each row's conditional of mu given tau, normalized
# by adaptive quadrature around the kernel (the Cauchy prior is evaluated
# exactly), never by the code under test.
.resolution_row_density <- function(seed, tau, x, log_prior) {
  v <- seed$sei^2 + tau^2
  prec <- seed$K / v
  center <- seed$ybar
  s <- 1 / sqrt(prec)
  lq <- function(m) -.5 * prec * (m - center)^2 + log_prior(m)
  ref <- lq(center)
  Z <- stats::integrate(function(m) exp(lq(m) - ref), center - 60 * s, center + 60 * s,
                        rel.tol = 1e-13, subdivisions = 5000L)$value
  exp(lq(x) - ref) / Z
}




testthat::test_that("quarter-step lattice under-resolution reaches the existing public gate", {
  prior <- BayesTools::prior("cauchy", list(location = 0, scale = 1))
  seed <- .resolution_seed_fit(8, prior, -1.85, 1)
  reference <- ((seed$S - 1) * .resolution_row_density(seed, 0, seed$value, function(x) dcauchy(x, 0, 1, log = TRUE)) +
    .resolution_row_density(seed, 1, seed$value, function(x) dcauchy(x, 0, 1, log = TRUE))) / seed$S
  testthat::expect_identical(seed$S, 500L)
  testthat::expect_true(is.finite(reference) && reference > 0)
  result <- tryCatch(RoBMA::hypothesis(seed$object, paste0("mu = ", sprintf("%.17g", seed$value)),
    density_method = "qCMDE", columns = "all"), error = identity)
  if (!inherits(result, "error")) {
    raw <- unname(attr(result, "raw_BF"))
    estimate <- dcauchy(seed$value) / raw
    cat("D5_FROZEN_FALSE_ACCEPT", "reference", reference, "estimate", estimate,
      "relative_error", estimate / reference - 1, "BF", raw, "\n")
    print(as.data.frame(RoBMA::density_diagnostics(result)))
  }
  testthat::expect_true(inherits(result, "RoBMA_density_ordinate_error"))
  if (inherits(result, "RoBMA_density_ordinate_error")) {
    diagnostic <- as.data.frame(result[["density_diagnostics"]])
    context <- .iwmde_context(seed$object)
    estimated <- .iwmde_estimate(context, "mu", "qCMDE", seed$control,
      outputs = "ordinate", values = seed$value)
    density <- estimated[["diagnostics"]][["ordinate"]][["iwmde"]]
    plan <- estimated[["plan"]]
    execution <- .iwmde_plan_row_execution(context, plan)
    tau <- context[["posterior_samples"]][execution[["active_rows"]], "tau"]
    # Independent physical likelihood width for the two observations. This
    # fixture has an identity chart, so its actual nested step is physical.
    width <- sqrt((seed$sei^2 + tau^2) / seed$K)
    step <- diff(density[["normalization_refined_range"]]) /
      (density[["normalization_refined_points"]] - 1L)
    e <- 2 * exp(-2 * pi^2 * (width / step)^2)
    row_error <- e / (1 - e)
    log_q <- .iwmde_log_q_grid(context, "mu", seed$value, execution[["row_states"]], plan[["replacement"]])
    log_contribution <- as.numeric(log_q) - density[["log_normalizer"]]
    shares <- exp(log_contribution - max(log_contribution))
    expected_resolution <- sum(shares * row_error) / sum(shares)
    expect_lte(abs(diagnostic[["lattice_resolution_relative_error"]] - expected_resolution), 1e-12)
    expect_lte(abs(diagnostic[["stability_relative_error"]] - diagnostic[["ordinate_relative_change"]] -
      diagnostic[["truncation_ordinate_bound"]] - diagnostic[["lattice_resolution_relative_error"]]), 1e-12)
    expect_identical(diagnostic[["stability_metric"]], "ordinate_error_estimate")
    expect_match(conditionMessage(result), "error estimate", fixed = TRUE)
  }
})

testthat::test_that("the resolved divisor-four window remains publicly available", {
  prior <- BayesTools::prior("cauchy", list(location = 0, scale = 1))
  seed <- .resolution_seed_fit(4, prior, -1.85, 1)
  result <- tryCatch(RoBMA::hypothesis(seed$object, paste0("mu = ", sprintf("%.17g", seed$value)),
    density_method = "qCMDE", columns = "all"), error = identity)
  testthat::expect_false(inherits(result, "error"))
  if (!inherits(result, "error")) testthat::expect_true(is.finite(unname(attr(result, "raw_BF"))))
})

testthat::test_that("the sharp Cauchy core retains its existing diagnostic refusal", {
  # Exact original report08/p06 control: the likelihood-kernel width can miss
  # prior curvature, so the additional resolution estimate is not a theorem.
  seed <- .resolution_seed_fit(1, BayesTools::prior("cauchy", list(location = 0, scale = .01)), .25, 1)
  testthat::expect_error(RoBMA::hypothesis(seed$object,
    paste0("mu = ", sprintf("%.17g", seed$value)), density_method = "qCMDE"),
    class = "RoBMA_density_ordinate_error")
})

testthat::test_that("Normal-prior closed forms retain the exact conditional mixture", {
  seed <- .resolution_seed_fit(8, BayesTools::prior("normal", list(mean = 0, sd = 1)), -1.85, 1)
  result <- RoBMA::hypothesis(seed$object, paste0("mu = ", sprintf("%.17g", seed$value)),
    density_method = "qCMDE", columns = "all")
  precision <- seed$K / (seed$sei^2 + c(0, 1)^2)
  mean <- precision / (1 + precision) * seed$ybar
  sd <- sqrt(1 / (1 + precision))
  reference <- ((seed$S - 1) * dnorm(seed$value, mean[1L], sd[1L]) +
    dnorm(seed$value, mean[2L], sd[2L])) / seed$S
  estimate <- dnorm(seed$value) / unname(attr(result, "raw_BF"))
  testthat::expect_lte(abs(estimate / reference - 1), 1e-11)
  diagnostic <- as.data.frame(RoBMA::density_diagnostics(result))
  testthat::expect_identical(diagnostic$normalization_truncation_status, "exact")
  testthat::expect_identical(diagnostic$truncation_ordinate_bound, 0)
  expect_identical(diagnostic[["lattice_resolution_relative_error"]], 0)
})

test_that("lattice widths use existing nested maxima and log-chart Jacobians", {

  z <- seq(-.1, .1, by = .001)
  widths <- c(.0004, .01)
  grid <- list(z = z, log_jacobian = rep(0, length(z)),
    log_q = vapply(widths, function(sd) -.5 * (z / sd)^2, numeric(length(z))))
  curvature <- .iwmde_qcmde_row_resolution(grid, NULL, .001, 1:2)
  expect_lte(max(abs(curvature[["width"]] / widths - 1)), 1e-12)
  e <- 2 * exp(-2 * pi^2 * (widths / .001)^2)
  expect_lte(max(abs(curvature[["relative_error"]] - e / (1 - e))), 1e-12)
  z <- seq(-5, 1, by = .1)
  grid <- list(z = z, log_jacobian = z, log_q = matrix(-.5 * 1e4 * (exp(z) + .01)^2, ncol = 1L))
  gaussian <- .iwmde_qcmde_row_resolution(grid, list(quadratic = 1e4), .1, 1L)
  maximum <- which.max(as.numeric(grid[["log_q"]]) + z)
  expect_equal(gaussian[["width"]], .01 / exp(z[maximum]), tolerance = 1e-14)
  expect_true(is.finite(gaussian[["relative_error"]]))
  flat <- list(log_q = matrix(0, 5L, 1L), log_jacobian = rep(0, 5L))
  expect_true(is.na(.iwmde_qcmde_row_resolution(flat, NULL, 1, 1L)[["relative_error"]]))
  expect_identical(.iwmde_qcmde_row_resolution(flat, list(quadratic = 1e12), 1, 1L)[["relative_error"]], Inf)
})

test_that("resolution errors use actual ordinate shares and preserve unavailable widths", {

  contributions <- rbind(c(1, 9), c(3, 1))
  normalizers <- log(c(2, 3))
  log_q <- sweep(log(contributions), 2L, normalizers, "+")
  resolution <- list(relative_error = c(.1, .4))
  error <- .iwmde_qcmde_resolution_ordinate_error(log_q, normalizers, resolution)
  expect_equal(error, c(.37, .175), tolerance = 1e-15)
  expect_true(all(is.na(.iwmde_qcmde_resolution_ordinate_error(log_q, normalizers,
    list(relative_error = c(NA_real_, .4))))))
  expect_identical(.iwmde_qcmde_resolution_ordinate_error(matrix(c(-Inf, 0), nrow = 1), c(0, 0),
    list(relative_error = c(NA_real_, .4))), .4)
})
