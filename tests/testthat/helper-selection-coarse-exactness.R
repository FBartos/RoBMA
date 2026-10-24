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
  weight <- apply(matrix(rule$weights[grid], nrow(grid), setup$rank), 1L, prod)
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


# Mean and variance of the posterior with the given log normalizer (a function of
# mu), by QUADPACK over mu, and their error estimates in the manner of
# test-00-selection-corrected-sampler.R (QUADPACK absolute errors of the moment
# integrals, plus the largest relative error of the normalizer carried through).
.coarse_exactness_moments <- function(setup, log_normalizer, tracker) {

  cache <- new.env()
  normalizer <- function(mu) {
    vapply(mu, function(value) {
      key <- sprintf("%.17g", value)
      if (is.null(cache[[key]])) { cache[[key]] <- log_normalizer(value); tracker$evaluations <- (if (is.null(tracker$evaluations)) 0 else tracker$evaluations) + 1 }
      cache[[key]]
    }, numeric(1L))
  }
  centre  <- -1
  density <- function(mu) {
    exp(stats::dt(mu, df = 4, log = TRUE) - setup$information * mu^2 / 2 -
          normalizer(mu) - (stats::dt(centre, df = 4, log = TRUE) -
          setup$information * centre^2 / 2 - normalizer(centre)))
  }
  moment <- function(power) {
    result <- stats::integrate(function(mu) density(mu) * mu^power, -14, 8,
                               rel.tol = 1e-6, abs.tol = 1e-9, subdivisions = 200L)
    stopifnot(identical(result$message, "OK"))
    result
  }
  integrals <- lapply(0:2, moment)
  values <- vapply(integrals, `[[`, numeric(1L), "value")
  errors <- vapply(integrals, `[[`, numeric(1L), "abs.error")
  mean     <- values[2L] / values[1L]
  variance <- values[3L] / values[1L] - mean^2
  inner    <- tracker$relative_error
  moment_error <- (errors[2:3] + abs(c(mean, values[3L] / values[1L])) * errors[1L]) /
    (values[1L] - errors[1L]) + 2 * inner / (1 - inner) *
    c(sqrt(values[3L] / values[1L]), values[3L] / values[1L])
  list(mean = mean, variance = variance, mean_error = moment_error[1L],
       variance_error = moment_error[2L] + 2 * abs(mean) * moment_error[1L] +
         moment_error[1L]^2)
}


# A short corrected chain of mu on the rjags model of the block, with the
# schedule of the setup. Returns the draws and the sampler counters of them.
.coarse_exactness_chain <- function(setup, draws, seed = 74891L) {

  old <- RoBMA::RoBMA.options()[c("selection.sampler", "selection.coarse_max_rules")]
  on.exit(do.call(RoBMA::RoBMA.options, old), add = TRUE)
  RoBMA::RoBMA.options(selection.sampler = "coarse_corrected", selection.coarse_max_rules = 3L)
  rows   <- setup$rows
  rank   <- setup$rank
  design <- BayesTools::selection_qmc_design(2L * rank, 512L, 8L, 9049L)
  data <- list(y = rep(0, rows), residual_sd = setup$residual_sd, loading = setup$loading,
    sei = rep(1, rows), omega = setup$omega, obs_bin = rep(1, rows),
    z_lower = c(0, -1e300), z_upper = c(1e300, 0), qmc = as.numeric(design),
    nqmc = length(design), nodes = setup$rules$nodes, log_weights = setup$rules$log_weights,
    orders = setup$rules$orders, nq = length(setup$rules$nodes), nr = length(setup$rules$orders))
  if (rank == 1L) data$loading <- as.numeric(setup$loading)
  distribution <- if (rank == 1L) "dselnorm_cluster_step" else "dselnorm_factor_step"
  loading_term <- if (rank == 1L) paste0("loading[1:", rows, "]") else
    paste0("loading[1:", rows, ", 1:", rank, "]")
  syntax <- paste0("model {
    mu ~ dt(0, 1, 4)
    for (j in 1:", rows, ") { means[j] <- mu }
    y[1:", rows, "] ~ ", distribution, "(means[1:", rows, "], residual_sd[1:", rows, "], ",
    loading_term, ", sei[1:", rows, "], omega[1:2], z_lower[1:2], z_upper[1:2],
      obs_bin[1:", rows, "], 1, 1, 1, nodes[1:nq], log_weights[1:nq], orders[1:nr],
      qmc[1:nqmc], 64, 512, 8, 0.005, 0)
  }")
  connection <- textConnection(syntax)
  on.exit(close(connection), add = TRUE)
  model <- rjags::jags.model(connection, data = data,
    inits = list(mu = -.25, .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed),
    n.chains = 1L, n.adapt = 0L, quiet = TRUE)
  assigned <- names(rjags::list.samplers(model))
  for (batch in seq_len(5L)) ready <- rjags::adapt(model, 100L, end.adaptation = FALSE)
  stopifnot(isTRUE(ready), isTRUE(rjags::adapt(model, 0L, end.adaptation = TRUE)))
  for (batch in seq_len(5L)) stats::update(model, n.iter = 100L, progress.bar = "none")
  RoBMA::selection_sampler_info(clear = TRUE)
  values <- as.numeric(as.matrix(rjags::coda.samples(model, "mu", n.iter = draws,
    thin = 1L, progress.bar = "none")))
  list(draws = values, samplers = assigned, counters = RoBMA::selection_sampler_info())
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
