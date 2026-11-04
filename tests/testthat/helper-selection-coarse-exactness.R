# Exactness of the corrected sampler on cluster and factor blocks. The model is
# a location mu ~ t(4) with y = 0 observed on the rows of one selection block,
# whose rows share the location mu, so the exact posterior of mu is
#   p(mu | y) ~ t4(mu) exp(-I mu^2 / 2) / A(mu),   I = 1' S^-1 1,
# with S = diag(sd^2) + L L' the Gaussian covariance of the block and
#   A(mu) = E_z prod_i P_i(mu + L_i z)
# its selection normalizer (the rows are observed in the bin of weight one, so
# the observed weights are constants). The references integrate A with QUADPACK
# (nested over the latent axes) and the surrogate-only target with the
# normalizer of the Gauss-Hermite rule the sampler's surrogate uses, both in
# base R.

# crude schedules whose rule 3 has three nodes (per axis)
.coarse_exactness_setup <- function(kind) {

  omega <- c(1, .05)
  if (kind == "cluster") {
    # equicorrelated, rho = .8
    loading     <- matrix(sqrt(.8), 4L, 1L)
    residual_sd <- rep(sqrt(.2), 4L)
    orders      <- c(1L, 2L, 3L, 7L, 15L, 31L, 63L, 127L, 255L, 511L, 1023L)
  } else {
    # two rows load on the first axis only, two on both: non-parallel rows
    loading     <- cbind(c(.7, .7, .6, .6), c(0, 0, .5, .5))
    residual_sd <- c(.5, .5, .45, .45)
    orders      <- c(1L, 2L, 3L, 5L, 7L, 9L, 11L, 15L, 21L, 31L, 43L, 63L, 95L)
  }
  covariance <- diag(residual_sd^2) + loading %*% t(loading)
  list(kind = kind, omega = omega, loading = loading, residual_sd = residual_sd,
       rows = nrow(loading), rank = ncol(loading), orders = orders,
       rules = .coarse_rules(orders),
       information = sum(solve(covariance, rep(1, nrow(loading)))))
}


# Row probabilities at the latent points `latent` (points x rank) of location mu,
# multiplied over the rows (identical rows once, to their multiplicity).
.coarse_exactness_product <- function(setup, mu, latent) {

  keys    <- apply(cbind(setup$loading, setup$residual_sd), 1L, paste, collapse = " ")
  unique_rows <- which(!duplicated(keys))
  product <- rep(1, nrow(latent))
  for (i in unique_rows) {
    count <- sum(keys == keys[i])
    probability <- .coarse_reference_probability(
      mu + as.numeric(latent %*% setup$loading[i, ]), setup$residual_sd[i], 1,
      setup$omega, c(0, -1e300), c(1e300, 0), 1)
    product <- product * probability^count
  }
  product
}


# log A_j(mu): the tensor Gauss-Hermite rule of order `order` per axis.
.coarse_exactness_rule_log_normalizer <- function(setup, mu, order) {

  rule   <- .coarse_reference_rule(order)
  grid   <- as.matrix(expand.grid(rep(list(seq_len(order)), setup$rank)))
  weight <- Reduce(`*`, lapply(seq_len(setup$rank), function(axis) rule$weights[grid[, axis]]))
  latent <- matrix(rule$nodes[grid], nrow(grid), setup$rank)
  log(sum(weight * .coarse_exactness_product(setup, mu, latent)))
}


# log A(mu) by QUADPACK, nested over the axes; `tracker` collects the largest
# relative error estimate of the integrals.
.coarse_exactness_quadpack_log_normalizer <- function(setup, mu, tracker) {

  integral <- function(f) {
    # Beyond 10 standard deviations the latent normal has mass below 1e-22.
    result <- stats::integrate(f, -10, 10, rel.tol = 1e-7, abs.tol = 1e-12,
                               subdivisions = 200L)
    stopifnot(identical(result$message, "OK"), result$value > 0)
    tracker$relative_error <- max(tracker$relative_error, result$abs.error / result$value)
    result$value
  }
  if (setup$rank == 1L) {
    value <- integral(function(z) {
      stats::dnorm(z) * .coarse_exactness_product(setup, mu, matrix(z, ncol = 1L))
    })
  } else {
    value <- integral(function(z1) {
      vapply(z1, function(first) {
        stats::dnorm(first) * integral(function(z2) {
          stats::dnorm(z2) * .coarse_exactness_product(setup, mu, cbind(first, z2))
        })
      }, numeric(1L))
    })
  }
  log(value)
}


# Mean and variance of the posterior of mu with the given log normalizer (a function
# of mu), by QUADPACK over mu, with the error allowances of .coarse_moments_1d.
.coarse_exactness_moments <- function(setup, log_normalizer, tracker) {

  cache <- new.env()
  normalizer <- function(mu) {
    key <- sprintf("%.17g", mu)
    if (is.null(cache[[key]])) cache[[key]] <- log_normalizer(mu)
    cache[[key]]
  }
  log_density <- function(mu) {
    stats::dt(mu, df = 4, log = TRUE) - setup$information * mu^2 / 2 - normalizer(mu)
  }
  .coarse_moments_1d(log_density, -14, 8, tracker, rel.tol = 1e-6, abs.tol = 1e-9)
}


# The rjags model of one selection block with the pieces a test varies: the
# `definitions` (priors and deterministic nodes the arguments refer to), the
# syntax of the residual SD and loading arguments, the constants of the kernel
# mode and vector rule, and the data. The corrected sampler is switched on while
# the model compiles and the options are restored. Returns the compiled model and
# the samplers it was assigned.
.coarse_hand_model <- function(kind, rows, rank, rules, data, definitions,
                               residual_sd = paste0("residual_sd[1:", rows, "]"),
                               loading = if (rank == 1L) paste0("loading[1:", rows, "]") else
                                 paste0("loading[1:", rows, ", 1:", rank, "]"),
                               mode = 1, vector_rule = 0, seed = 74891L, inits = list()) {

  old <- RoBMA::RoBMA.options()[c("selection.sampler", "selection.coarse_max_rules")]
  on.exit(do.call(RoBMA::RoBMA.options, old), add = TRUE)
  RoBMA::RoBMA.options(selection.sampler = "coarse_corrected", selection.coarse_max_rules = 3L)
  design <- BayesTools::selection_qmc_design(2L * rank, 512L, 8L, 9049L)
  base <- list(y = rep(0, rows), sei = rep(1, rows), omega = c(1, .05),
    obs_bin = rep(1, rows), z_lower = c(0, -1e300), z_upper = c(1e300, 0),
    qmc = as.numeric(design), nqmc = length(design), nodes = rules$nodes,
    log_weights = rules$log_weights, orders = rules$orders, nq = length(rules$nodes),
    nr = length(rules$orders))
  data <- utils::modifyList(base, data)
  distribution <- if (kind == "cluster") "dselnorm_cluster_step" else "dselnorm_factor_step"
  syntax <- paste0("model {\n", definitions, "\n  y[1:", rows, "] ~ ", distribution,
    "(means[1:", rows, "], ", residual_sd, ", ", loading, ", sei[1:", rows, "], omega[1:2], ",
    "z_lower[1:2], z_upper[1:2], obs_bin[1:", rows, "], 1, 1, ", mode, ", nodes[1:nq], ",
    "log_weights[1:nq], orders[1:nr], qmc[1:nqmc], 64, 512, 8, 0.005, ", vector_rule, ")\n}")
  connection <- textConnection(syntax)
  on.exit(close(connection), add = TRUE)
  model <- rjags::jags.model(connection, data = data,
    inits = c(inits, list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed)),
    n.chains = 1L, n.adapt = 0L, quiet = TRUE)
  list(model = model, samplers = names(rjags::list.samplers(model)))
}


# Adapts, burns in and draws `draws` values of `monitor`; returns them (a matrix for
# several monitored nodes) with the sampler counters of the draws.
.coarse_hand_chain <- function(hand, monitor, draws) {

  model <- hand$model
  for (batch in seq_len(5L)) ready <- rjags::adapt(model, 100L, end.adaptation = FALSE)
  stopifnot(isTRUE(ready), isTRUE(rjags::adapt(model, 0L, end.adaptation = TRUE)))
  for (batch in seq_len(5L)) stats::update(model, n.iter = 100L, progress.bar = "none")
  RoBMA::selection_sampler_info(clear = TRUE)
  values <- as.matrix(rjags::coda.samples(model, monitor, n.iter = draws,
    thin = 1L, progress.bar = "none"))
  list(draws = if (length(monitor) == 1L) as.numeric(values) else values,
       samplers = hand$samplers, counters = RoBMA::selection_sampler_info())
}


# A short corrected chain of mu on the rjags model of the block, with the
# schedule of the setup. Returns the draws and the sampler counters of them.
.coarse_exactness_chain <- function(setup, draws, seed = 74891L) {

  rows <- setup$rows
  rank <- setup$rank
  data <- list(residual_sd = setup$residual_sd, loading = setup$loading, omega = setup$omega)
  if (rank == 1L) data$loading <- as.numeric(setup$loading)
  hand <- .coarse_hand_model(setup$kind, rows, rank, setup$rules, data,
    definitions = paste0("  mu ~ dt(0, 1, 4)\n  for (j in 1:", rows, ") { means[j] <- mu }"),
    seed = seed, inits = list(mu = -.25))
  .coarse_hand_chain(hand, "mu", draws)
}


# Mean and variance of a density on a finite interval given by its log (a function
# of the parameter, up to a constant), by QUADPACK, with the error allowances of
# .coarse_moment_summary (the relative error of the normalizer integrals is read
# from `tracker$relative_error` after the integration).
.coarse_moments_1d <- function(log_density, lower, upper, tracker,
                               rel.tol = 1e-8, abs.tol = 1e-12) {

  raw <- .coarse_raw_moments(log_density, lower, upper, log_density((lower + upper) / 2),
                             rel.tol = rel.tol, abs.tol = abs.tol)
  .coarse_moment_summary(raw$values, raw$errors, tracker$relative_error)
}


# Gaussian log density of x for the covariance diag(sd^2) + L L' around mean.
.coarse_gaussian_log_density <- function(x, mean, sd, loading) {

  loading    <- as.matrix(loading)
  covariance <- diag(sd^2, length(x)) + loading %*% t(loading)
  factor     <- chol(covariance)
  centered   <- backsolve(factor, x - mean, transpose = TRUE)
  -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(factor))) + sum(centered^2))
}


# For each state mu: the log normalizer the exact kernel used (numerator minus
# exact density through the probe), its distance from the rule the surrogate
# integrates over (rule `max_rules`) and the first rule of the schedule whose
# independent Gauss-Hermite normalizer it equals (NA when none does: a
# randomized or averaged result).
.coarse_exactness_rungs <- function(setup, states, max_rules = 3L, agreement = 1e-9) {

  rows <- setup$rows
  do.call(rbind, lapply(states, function(mu) {
    block <- .coarse_block(setup$kind, mean = rep(mu, rows), residual_sd = setup$residual_sd,
      loading = setup$loading, sei = rep(1, rows), omega = setup$omega,
      obs_bin = rep(1L, rows), x = rep(0, rows), rules = setup$rules)
    probe  <- .coarse_probe(block, max_rules)
    exact_log_a <- probe[["numerator"]] - probe[["exact"]]
    rung_log_a  <- vapply(setup$orders, function(order) {
      .coarse_exactness_rule_log_normalizer(setup, mu, order)
    }, numeric(1L))
    matched <- which(abs(rung_log_a - exact_log_a) < agreement)
    data.frame(mu = mu, exact_log_a = exact_log_a, surrogate_log_a = probe[["log_normalizer"]],
      surrogate_vs_rule = abs(probe[["log_normalizer"]] - rung_log_a[max_rules]),
      gap = abs(exact_log_a - rung_log_a[max_rules]),
      accepted_rule = if (length(matched)) matched[1L] else NA_integer_)
  }))
}


# Raw moments (integrals of exp(log_density - reference) times x^0, x^1, x^2 over
# [lower, upper]) by QUADPACK, with the absolute error estimates of the three.
.coarse_raw_moments <- function(log_density, lower, upper, reference,
                                rel.tol = 1e-8, abs.tol = 1e-12) {

  density <- function(x) exp(vapply(x, log_density, numeric(1L)) - reference)
  integrals <- lapply(0:2, function(power) {
    result <- stats::integrate(function(x) density(x) * x^power, lower, upper,
                               rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = 200L)
    stopifnot(identical(result$message, "OK"))
    result
  })
  list(values = vapply(integrals, `[[`, numeric(1L), "value"),
       errors = vapply(integrals, `[[`, numeric(1L), "abs.error"))
}


# log A of the one-sided best-weight rule of a rank-one block with two bins cut at
# zero (weights 1 and omega_2) for rows of location mu: the weight of a vector is
# that of its largest z, so A = 1 - (1 - omega_2) P(every row at or below zero),
# and that probability is a one-dimensional integral over the shared latent axis.
.coarse_exactness_best_log_normalizer <- function(setup, mu, tracker) {

  result <- stats::integrate(function(z) {
    probability <- rep(1, length(z))
    for (i in seq_len(setup$rows)) {
      probability <- probability * stats::pnorm(-(mu + setup$loading[i, 1L] * z) / setup$residual_sd[i])
    }
    stats::dnorm(z) * probability
  }, -10, 10, rel.tol = 1e-9, abs.tol = 1e-14, subdivisions = 200L)
  stopifnot(identical(result$message, "OK"))
  # the error of A is the error of the probability, scaled, relative to A itself
  normalizer <- 1 - (1 - setup$omega[2L]) * result$value
  tracker$relative_error <- max(tracker$relative_error,
                                (1 - setup$omega[2L]) * result$abs.error / normalizer)
  log(normalizer)
}


# Mean, variance and their error allowances from the raw moments (values and
# absolute error estimates of the integrals of 1, x and x^2, possibly summed over
# components) and the largest relative error of the normalizers inside them.
.coarse_moment_summary <- function(values, errors, relative_error) {

  mean   <- values[2L] / values[1L]
  second <- values[3L] / values[1L]
  moment_error <- (errors[2:3] + abs(c(mean, second)) * errors[1L]) / (values[1L] - errors[1L]) +
    2 * relative_error / (1 - relative_error) * c(sqrt(second), second)
  list(mean = mean, variance = second - mean^2, mean_error = moment_error[1L],
       variance_error = moment_error[2L] + 2 * abs(mean) * moment_error[1L] + moment_error[1L]^2)
}
