# Independent declared-model and rounded-input references are documented beside
# the exact binary fixture. These helpers check model input identities and
# the actual replay's physical mean/variance perturbation; they do not evaluate an
# alternative likelihood or change the original output-scale allowance.
.ordered_gaussian_reference <- function(case) {

  readRDS(testthat::test_path("fixtures", "ordered-gaussian-reference.rds"))[[case]]
}


.ordered_gaussian_reference_inputs <- function(reference, block) {

  input <- reference[["inputs"]]
  expect_identical(block[["setup"]][["yi"]], input[["yi"]])
  for (field in c("mu0", "mu_a", "mu_b", "mu_c", "variance")) {
    expect_identical(as.numeric(block[[field]]), as.numeric(input[[field]]))
  }
  expect_identical(as.numeric(block[["weights"]]), as.numeric(input[["powers"]]))
}


.ordered_gaussian_replay_matches <- function(reference, mean, sd) {

  identical(unname(mean), reference[["replay_mean"]]) && identical(unname(sd), reference[["replay_sd"]])
}


.ordered_gaussian_rounding_bound <- function(reference, weights, variance, mean, sd, selection = NULL) {

  eps <- .Machine$double.eps
  # Four eps cover subtraction, absolute value, radius addition and the outward
  # multiplication. The exact rational means are independently enclosed in the
  # fixture; neither their replay rounding nor the BLAS sum is prescribed.
  upper <- function(value) value * (1 + 4 * eps)
  midpoint <- reference[["declared_mean"]]
  radius <- reference[["declared_mean_radius"]]
  delta <- upper(upper(abs(mean - midpoint)) + radius)
  yi <- matrix(reference[["inputs"]][["yi"]], nrow(mean), ncol(mean), byrow = TRUE)
  residual <- upper(upper(abs(yi - midpoint)) + radius)
  replay_residual <- upper(abs(yi - mean))
  v <- matrix(variance, nrow(delta), ncol(delta), byrow = TRUE)
  squared_sd <- sd^2
  # A full eps bounds rounding of the exact square of the actual input SD.
  square_error <- upper(eps * squared_sd)
  dv <- upper(upper(abs(squared_sd - v)) + square_error)
  lower_v <- pmin(v, (squared_sd - square_error) * (1 - 4 * eps))
  stopifnot(all(is.finite(delta)), all(is.finite(dv)), all(lower_v > 0))
  # |L(m+delta,v)-L(m,v)| plus the variance change at the actual replay mean.
  gaussian <- residual * delta / v + delta^2 / (2 * v) +
    .5 * (log1p(dv / lower_v) + replay_residual^2 * dv / (v * lower_v))
  selected <- 0
  if (!is.null(selection)) {
    omega <- as.numeric(selection[["omega"]])
    variation <- sum(abs(diff(omega)))
    minimum <- min(omega)
    # Z >= min(omega); telescoping step-CDF sensitivities bound both changes.
    selected <- variation * delta / (sqrt(2 * pi * lower_v) * minimum) +
      variation * dv / (2 * lower_v * sqrt(2 * pi * exp(1)) * minimum)
  }
  # The finite operation margin covers the per-cell arithmetic and weighted
  # row reduction. It propagates input rounding only, with no output-scale
  # evaluation floor added to the original likelihood allowance.
  rowSums(sweep(gaussian + selected, 2L, weights, "*")) * (1 + (64 + 2 * ncol(delta)) * eps)
}


# Finite-range arithmetic accounting for the literal Normal source tests.
# The named eight-eps relative log-CDF allowance is an ordinary libm allowance, not a
# rigorous enclosure of arbitrary pnorm implementations or extreme inputs.
.ordered_normal_interval_arithmetic_allowance <- function(lower, upper, mean, sd) {

  eps <- .Machine$double.eps
  gamma2 <- 2 * eps / (1 - 2 * eps)
  stopifnot(is.finite(mean), is.finite(sd), sd > 0, lower < upper)
  zero <- c(standardization = 0, log_cdf_libm = 0, difference_rounding = 0, exponential_subtraction = 0)
  if (identical(c(lower, upper), c(-Inf, Inf))) return(list(log_error = 0, candidates = list(left = zero, right = zero)))
  endpoints <- c(lower, upper)
  z <- (endpoints - mean) / sd
  perturbation <- rep(0, 2L)
  finite <- is.finite(endpoints)
  perturbation[finite] <- gamma2 * ((abs(endpoints[finite]) + abs(mean)) / sd + abs(z[finite]))
  maximum_pdf <- vapply(seq_along(z), function(i) {
    if (!finite[[i]]) return(0)
    closest <- if (z[[i]] - perturbation[[i]] <= 0 && z[[i]] + perturbation[[i]] >= 0) 0 else
      min(abs(z[[i]] - perturbation[[i]]), abs(z[[i]] + perturbation[[i]]))
    stats::dnorm(closest)
  }, numeric(1))
  cdf_perturbation <- maximum_pdf * perturbation
  candidates <- list()
  for (tail in c("left", "right")) {
    logs <- stats::pnorm(endpoints, mean, sd, lower.tail = tail == "left", log.p = TRUE)
    order <- if (tail == "left") c(2L, 1L) else c(1L, 2L)
    log_x <- logs[[order[[1L]]]]
    log_y <- logs[[order[[2L]]]]
    libm <- rep(0, 2L)
    libm[finite] <- 8 * eps * abs(logs[finite])
    delta <- log_y - log_x
    gap <- -expm1(delta)
    ratio <- exp(delta)
    log_mass <- log_x + log(gap)
    stopifnot(is.finite(log_mass), is.finite(gap), gap > 0)
    terms <- c(standardization = sum(cdf_perturbation) / exp(log_mass),
      log_cdf_libm = (libm[[order[[1L]]]] + ratio * libm[[order[[2L]]]]) / gap,
      difference_rounding = if (is.finite(log_y)) eps * (abs(log_x) + abs(log_y)) * ratio / gap else 0,
      exponential_subtraction = eps * ratio / gap)
    stopifnot(all(is.finite(terms)), all(terms >= 0), sum(terms) < 1)
    candidates[[tail]] <- c(terms, log_error = -log1p(-sum(terms)))
  }
  list(log_error = max(vapply(candidates, `[[`, numeric(1), "log_error")),
    candidates = candidates, standardized = z, endpoint_perturbation = perturbation)
}
