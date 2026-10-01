.selection_corrected_sampler_regression <- function(
    sampler_stats = function() RoBMA::selection_sampler_info(), draws = 5000L) {

  stopifnot(is.function(sampler_stats), length(draws) == 1L,
    is.finite(draws), draws >= 2000L, draws == as.integer(draws))
  stopifnot(requireNamespace("RoBMA", quietly = TRUE),
    requireNamespace("rjags", quietly = TRUE),
    requireNamespace("coda", quietly = TRUE),
    requireNamespace("BayesTools", quietly = TRUE))
  option_names <- c("selection.sampler", "selection.coarse_grid",
                    "selection.cache_max_bytes", "selection.coarse_max_rules")
  old_options <- RoBMA::RoBMA.options()
  stopifnot(all(option_names %in% names(old_options)))
  on.exit(do.call(RoBMA::RoBMA.options, old_options[option_names]), add = TRUE)
  grid <- c(mean = .35, variance = .10, log_weight = .25)
  configure <- function(value, capacity) {

    RoBMA::RoBMA.options(selection.sampler = "coarse_corrected",
      selection.coarse_grid = value, selection.cache_max_bytes = capacity,
      selection.coarse_max_rules = 3L)
    invisible(NULL)
  }

  rho <- .4
  omega <- .2
  covariance <- matrix(rho, 3L, 3L)
  diag(covariance) <- 1
  rules <- RoBMA:::.selection_joint_cluster_quadrature_rules(
    c(7L, 15L, 31L, 63L, 127L, 255L, 511L, 1023L))
  design <- BayesTools::selection_qmc_design(6L, 512L, 8L, 9049L)
  data <- list(y = rep(0, 3L), covariance = covariance[lower.tri(covariance, diag = TRUE)],
    sei = rep(1, 3L), omega = c(1, omega), obs_bin = rep(1, 3L),
    z_lower = c(0, -1e300), z_upper = c(1e300, 0),
    qmc = as.numeric(design), nqmc = length(design),
    nodes = rules$nodes, log_weights = rules$log_weights, orders = rules$orders,
    nq = length(rules$nodes), nr = length(rules$orders))
  syntax <- "model {
    mu ~ dt(0, 1, 4)
    for (j in 1:3) { means[j] <- mu }
    y[1:3] ~ dselnorm_mnorm_step(means[1:3], covariance[1:6], sei[1:3],
      omega[1:2], z_lower[1:2], z_upper[1:2], obs_bin[1:3],
      1, 1, 1, qmc[1:nqmc], 512, 8, 0.005, 0,
      nodes[1:nq], log_weights[1:nq], orders[1:nr])
  }"
  create <- function(seed) {

    connection <- textConnection(syntax)
    on.exit(close(connection), add = TRUE)
    fit <- rjags::jags.model(connection, data = data,
      inits = list(mu = -.25, .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed),
      n.chains = 1L, n.adapt = 0L, quiet = TRUE)
    assigned <- rjags::list.samplers(fit)
    stopifnot(identical(names(assigned), "RoBMA::CoarseCorrectedSlice"),
      length(assigned[[1L]]) == 1L,
      as.character(assigned[[1L]]) %in% c("mu", "means[1]"))
    # JAGS may name the common stochastic node through its means[1] alias.
    attr(fit, "sampled_node") <- as.character(assigned[[1L]])
    fit
  }
  adapt <- function(fit) {

    for (batch in seq_len(5L)) {
      ready <- rjags::adapt(fit, 100L, end.adaptation = FALSE)
    }
    stopifnot(isTRUE(ready))
    stopifnot(isTRUE(rjags::adapt(fit, 0L, end.adaptation = TRUE)))
    invisible(NULL)
  }
  sample <- function(fit, n) {

    values <- as.matrix(rjags::coda.samples(fit, c("mu", "means"), n.iter = n,
      thin = 1L, progress.bar = "none"))
    stopifnot(all(vapply(1:3, function(j) identical(as.numeric(values[, "mu"]),
      as.numeric(values[, paste0("means[", j, "]")]), num.eq = FALSE), logical(1L))))
    as.numeric(values[, "mu"])
  }

  # Two models own the same original policy/RNG state. A later option change
  # changes future factory construction and shared cache capacity only.
  configure(grid, 4 * 1024^2)
  first <- create(18431L)
  second <- create(18431L)
  adapt(first)
  adapt(second)
  first_draws <- sample(first, 100L)
  configure(c(mean = .8, variance = .4, log_weight = .7), 0)
  second_draws <- sample(second, 100L)
  stopifnot(identical(first_draws, second_draws, num.eq = FALSE))
  rm(first, second)

  # Independent Gaussian-factor integral. It uses base-R QUADPACK, not the
  # native normalizer, its covariance envelope, or its GH rule implementation.
  inner_error <- 0
  acceptance <- function(mu) {

    integral <- stats::integrate(function(z) {
      stats::dnorm(z) * (omega + (1 - omega) * stats::pnorm(
        (mu + sqrt(rho) * z) / sqrt(1 - rho)))^3
    }, -Inf, Inf, rel.tol = 1e-11, abs.tol = 1e-12, subdivisions = 200L)
    stopifnot(identical(integral$message, "OK"), integral$value > 0)
    inner_error <<- max(inner_error, integral$abs.error)
    integral$value
  }
  half_sum <- (1 + omega) / 2
  zero_reference <- half_sum^3 + half_sum * (1 - omega)^2 * 3 * asin(rho) / (2 * pi)
  stopifnot(abs(acceptance(0) - zero_reference) <= 5e-11)
  information <- 3 / (1 + 2 * rho)
  moment <- function(observable, subdivisions = 200L) {

    stats::integrate(function(mu) vapply(mu, function(value) {
      core <- exp(stats::dt(value, df = 4, log = TRUE) - information * value^2 / 2)
      if (core == 0) return(0)
      observable(value) * core / acceptance(value)
    }, numeric(1L)), -Inf, Inf, rel.tol = 1e-10, abs.tol = 1e-12,
      subdivisions = subdivisions)
  }
  integrals <- lapply(0:2, function(power) moment(function(value) value^power))
  stopifnot(all(vapply(integrals, function(x) identical(x$message, "OK"), logical(1L))))
  values <- vapply(integrals, `[[`, numeric(1L), "value")
  errors <- vapply(integrals, `[[`, numeric(1L), "abs.error")
  reference_mean <- values[2L] / values[1L]
  reference_second <- values[3L] / values[1L]
  reference_variance <- reference_second - reference_mean^2
  # A missing/sign-reversed correction can leave a periodic sawtooth in log
  # density. Its first sine harmonic can be visible when low moments look good.
  # The actual smooth posterior's expectation is integrated, never set to zero.
  mean_grid_width <- grid[["mean"]] * max(data$sei)
  oscillatory <- function(value) sin(2 * pi * value / mean_grid_width)
  oscillatory_integral <- moment(oscillatory, subdivisions = 500L)
  stopifnot(identical(oscillatory_integral$message, "OK"))
  reference_oscillatory <- oscillatory_integral$value / values[1L]
  # QUADPACK error estimates remain estimates. They are reported separately
  # from MCSE and are far below the intended stochastic comparison scale.
  stopifnot(values[1L] > errors[1L], reference_second > 0, inner_error < omega^3)
  relative_inner_error <- inner_error / omega^3
  moment_error <- (errors[2:3] + abs(c(reference_mean, reference_second)) * errors[1L]) /
    (values[1L] - errors[1L])
  moment_error <- moment_error + 2 * relative_inner_error / (1 - relative_inner_error) *
    c(sqrt(reference_second), reference_second)
  oscillatory_error <- (oscillatory_integral$abs.error +
    abs(reference_oscillatory) * errors[1L]) / (values[1L] - errors[1L]) +
    2 * relative_inner_error / (1 - relative_inner_error)
  reference_error <- c(moment_error[1L], moment_error[2L] +
    2 * abs(reference_mean) * moment_error[1L] + moment_error[1L]^2,
    oscillatory_error)
  stopifnot(reference_variance > 0, all(is.finite(reference_error)))

  configure(grid, 4 * 1024^2)
  fit <- create(74891L)
  adapt(fit)
  for (batch in seq_len(5L)) stats::update(fit, n.iter = 100L, progress.bar = "none")
  counters <- function() {

    values <- sampler_stats()
    required <- c("proposed", "accepted", "correction_rejected")
    stopifnot(all(required %in% names(values)))
    if (is.data.frame(values)) stopifnot(nrow(values) == 1L)
    stats::setNames(as.numeric(unlist(values[required], use.names = FALSE)), required)
  }
  before <- counters()
  values_drawn <- numeric(draws)
  for (start in seq.int(1L, draws, by = 100L)) {
    count <- min(100L, draws - start + 1L)
    values_drawn[seq.int(start, length.out = count)] <- sample(fit, count)
  }
  after <- counters()
  transitions <- after[names(before)] - before
  stopifnot(transitions[["proposed"]] == draws,
    transitions[["accepted"]] + transitions[["correction_rejected"]] == draws,
    transitions[["correction_rejected"]] > 0)
  oscillatory_draws <- oscillatory(values_drawn)
  estimate <- c(mean(values_drawn), stats::var(values_drawn), mean(oscillatory_draws))
  observables <- cbind(mean = values_drawn,
    variance = (values_drawn - mean(values_drawn))^2, sin_grid = oscillatory_draws)
  ess <- as.numeric(coda::effectiveSize(coda::mcmc(observables)))
  mcse <- sqrt(apply(observables, 2L, stats::var) / ess)
  reference <- c(reference_mean, reference_variance, reference_oscillatory)
  comparison <- data.frame(quantity = c("mean", "variance", "sin_grid"), estimate = estimate,
    reference = reference, difference = estimate - reference,
    MCSE = mcse, ESS = ess, reference_error_estimate = reference_error,
    standardized_error = (estimate - reference) / mcse)
  stopifnot(all(is.finite(mcse)), all(mcse > 0), all(ess >= 200),
    all(abs(estimate - reference) <= 6 * mcse + 2 * reference_error))
  invisible(list(comparison = comparison, transitions = transitions,
    immutable_grid_and_cache_capacity_bitwise = TRUE,
    independent_inner_error_estimate = inner_error,
    reference_normalizer = values[1L], grid = grid, draws = draws,
    oscillatory_reference = list(mean_grid_width = mean_grid_width,
      numerator = oscillatory_integral$value,
      numerator_absolute_error_estimate = oscillatory_integral$abs.error,
      normalized_expectation = reference_oscillatory,
      normalized_error_estimate = oscillatory_error)))
}


test_that("corrected selection sampling preserves priors, state and target", {

  skip_on_cran()
  result <- .selection_corrected_sampler_regression()
  expect_true(result$immutable_grid_and_cache_capacity_bitwise)
  expect_gt(result$transitions[["correction_rejected"]], 0)
  expect_true(all(abs(result$comparison$standardized_error) < 6))
})
