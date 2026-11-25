# Fixtures and independent references for the corrected block sampler
# (RoBMA::CoarseCorrectedBlock, src/samplers/CoarseCorrectedBlock.cc): blocks of
# scalar parameters declared by BayesTools' dbt_block markers, with a selection
# child whose surrogate is crude, in hand-written rjags models.
#
# The models share one structure. Members a and b enter the selection children
# only through their sum mu = a + b (a shared deterministic predictor), so the
# surrogate's error is a function of a + b (not separable in the members):
#   a ~ N(0, sa^2), w ~ N(a, 1) observed (an ordinary child),
#   b ~ N(0, s^2) with s ~ U(slo, shi) a member (a member that is a child of a
#     member) or b ~ N(0, sb^2) with sb fixed,
#   y_c ~ selection block with mean mu on every row, y_c = 0 observed (c = 1..K).
# Integrating a out, the joint posterior of (mu, s) is
#   p(mu, s | data) ~ U(s) N(mu; ma, va + s^2) exp(-K I mu^2 / 2) / A(mu)^K,
# with va = 1 / (1 / sa^2 + 1), ma = va w, I = 1' S^-1 1 of one block and A its
# selection normalizer, and a | mu, s is normal with variance va s^2 / (va + s^2)
# and mean (ma s^2 + mu va) / (va + s^2). The references integrate over mu (and s)
# with QUADPACK, A by QUADPACK (helper-selection-coarse-exactness.R).


# The native test hook of the block sampler: the fixed-order mutant (frozen into
# block samplers built afterwards) and the evaluation trace; returns the
# settings, the exhausted slice intervals and the trace records since the last call.
.block_test_hook <- function(fixed_order = NULL, trace = NULL) {

  .Call("RoBMA_corrected_block_test", fixed_order, trace, PACKAGE = "RoBMA")
}


# The trace records as a list of evaluations (kind 0 surrogate / 1 full, chain,
# member values, log density without the Jacobian).
.block_trace_records <- function(records) {

  out <- list()
  i <- 1L
  while (i <= length(records)) {
    d <- records[i + 2L]
    out[[length(out) + 1L]] <- list(kind = records[i], chain = records[i + 1L],
                                    x = records[i + 2L + seq_len(d)], value = records[i + 3L + d])
    i <- i + 4L + d
  }
  out
}


# One cluster selection child of `rows` rows on `means` (syntax), with the data
# names of .block_hand_model().
.block_cluster_child <- function(target, means, rows, vector_rule = 0, mode = 1) {

  paste0("  ", target, " ~ dselnorm_cluster_step(", means, ", residual_sd[1:", rows, "], loading[1:", rows,
         "], sei[1:", rows, "], omega[1:2], z_lower[1:2], z_upper[1:2], obs_bin[1:", rows, "], 1, 1, ", mode,
         ", nodes[1:nq], log_weights[1:nq], orders[1:nr], qmc[1:nqmc], 64, 512, 8, 0.005, ", vector_rule, ")")
}


# The data of cluster children with the crude schedule of
# .coarse_exactness_setup("cluster") (rule 3 has three nodes).
.block_cluster_data <- function(setup = .coarse_exactness_setup("cluster")) {

  design <- BayesTools::selection_qmc_design(2L, 512L, 8L, 9049L)
  rows   <- setup$rows
  list(residual_sd = setup$residual_sd, loading = as.numeric(setup$loading), sei = rep(1, rows),
       omega = setup$omega, obs_bin = rep(1, rows), z_lower = c(0, -1e300), z_upper = c(1e300, 0),
       qmc = as.numeric(design), nqmc = length(design), nodes = setup$rules$nodes,
       log_weights = setup$rules$log_weights, orders = setup$rules$orders,
       nq = length(setup$rules$nodes), nr = length(setup$rules$orders))
}


# Markers declaring `members` as block `id` (syntax and observed zeros).
.block_markers <- function(members, id = 1L) {

  names <- paste0("marker_", gsub("[^A-Za-z0-9]", "_", members))
  list(syntax = paste0("  ", names, " ~ dbt_block(", members, ", ", id, ")", collapse = "\n"),
       data = stats::setNames(as.list(rep(0, length(members))), names))
}


# Compiles a hand-written model with the corrected sampler options and the
# block sampler's test flags, restoring both; the per-fit counters and the test
# hook start empty. Returns the model and the sampler of each node.
.block_hand_model <- function(syntax, data, inits = list(), min_acceptance = 0, sampler = "slice.coarse",
                              seed = 74891L, chains = 1L, fixed_order = FALSE, grid = NULL) {

  BayesTools::BayesTools_load_JAGS_module(quiet = TRUE)
  old <- RoBMA::RoBMA.options()[c("selection.sampler", "selection.coarse_max_rules",
                                  "selection.coarse_min_acceptance", "selection.coarse_grid")]
  on.exit(do.call(RoBMA::RoBMA.options, old), add = TRUE)
  RoBMA::RoBMA.options(selection.sampler = sampler, selection.coarse_max_rules = 3L,
                       selection.coarse_min_acceptance = min_acceptance)
  if (!is.null(grid)) RoBMA::RoBMA.options(selection.coarse_grid = grid)
  RoBMA:::.selection_sampler_report_read()
  .block_test_hook(fixed_order = fixed_order, trace = FALSE)
  on.exit(.block_test_hook(fixed_order = FALSE, trace = FALSE), add = TRUE)
  connection <- textConnection(syntax)
  on.exit(close(connection), add = TRUE)
  inits <- lapply(seq_len(chains), function(chain) {
    c(inits, list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed + chain - 1L))
  })
  model <- rjags::jags.model(connection, data = data, inits = inits, n.chains = chains,
                             n.adapt = 0L, quiet = TRUE)
  samplers <- rjags::list.samplers(model)
  assigned <- unlist(lapply(seq_along(samplers), function(i) {
    stats::setNames(rep(names(samplers)[i], length(samplers[[i]])), samplers[[i]])
  }))
  list(model = model, samplers = assigned)
}


# Adapts for `adapt` iterations, ends the adaptation (each block decides its
# kernel there), burns in and draws `draws` values of `monitor`; returns the draws
# (a matrix, or a list of matrices for several chains), the process counters of
# the draws, the per-fit report of burn-in and draws, and the test hook's
# exhausted intervals and trace records of the draws.
.block_chain <- function(hand, monitor, draws, adapt = 0L, burnin = 500L) {

  model <- hand$model
  if (adapt > 0L) rjags::adapt(model, adapt, end.adaptation = FALSE)
  rjags::adapt(model, 0L, end.adaptation = TRUE)
  if (burnin > 0L) stats::update(model, n.iter = burnin, progress.bar = "none")
  RoBMA::selection_sampler_info(clear = TRUE)
  .block_test_hook()
  values <- rjags::coda.samples(model, monitor, n.iter = draws, thin = 1L, progress.bar = "none")
  hook <- .block_test_hook()
  list(draws = if (length(values) == 1L) as.matrix(values[[1L]]) else lapply(values, as.matrix),
       samplers = hand$samplers, counters = RoBMA::selection_sampler_info(),
       report = RoBMA:::.selection_sampler_report_read(), exhausted = hook$exhausted, trace = hook$records)
}


# Runs a compiled model: `adapt` adaptive updates, the end of adaptation when
# `end_adaptation` (otherwise the whole run adapts: rjags' update() alone does
# not end it) and `iterations` further updates; returns the error message, ""
# without one.
.block_run_message <- function(model, iterations, end_adaptation, adapt = 20L) {

  tryCatch({
    rjags::adapt(model, adapt, end.adaptation = FALSE)
    if (end_adaptation) rjags::adapt(model, 0L, end.adaptation = TRUE)
    stats::update(model, n.iter = iterations, progress.bar = "none")
    ""
  }, error = function(error) conditionMessage(error))
}


# ---------------------------------------------------------------------------- #
# The power fixture: members a, b, s (bounded, the parent of b) and r ~ Beta(rb,
# rb) (two-sided, mass at both bounds; independent of the rest), K crude cluster
# children on a + b.
# ---------------------------------------------------------------------------- #

.block_power_setup <- function(children = 3L) {

  cluster <- .coarse_exactness_setup("cluster")
  list(cluster = cluster, children = children, rows = cluster$rows, sa = 2, w = .5,
       slo = .5, shi = 2.5, rb = .3)
}

.block_power_model <- function(setup, members = c("a", "b", "s", "r"),
                               inits = list(a = .2, b = -.3, s = 1.1, r = .4), ...) {

  rows    <- setup$rows
  markers <- .block_markers(members)
  syntax  <- paste0("model {\n",
    "  a ~ dnorm(0, 1 / (sa * sa))\n  w ~ dnorm(a, 1)\n",
    "  s ~ dunif(slo, shi)\n  b ~ dnorm(0, 1 / (s * s))\n  r ~ dbeta(rb, rb)\n",
    "  for (j in 1:", rows, ") { means[j] <- a + b }\n",
    "  for (c in 1:K) {\n", .block_cluster_child(paste0("y[c, 1:", rows, "]"), paste0("means[1:", rows, "]"), rows),
    "\n  }\n", markers$syntax, "\n}")
  data <- c(.block_cluster_data(setup$cluster), markers$data,
            list(y = matrix(0, setup$children, rows), K = setup$children, sa = setup$sa, w = setup$w,
                 slo = setup$slo, shi = setup$shi, rb = setup$rb))
  .block_hand_model(syntax, data, inits = inits, ...)
}



# Moments of a, b, s and r under the power fixture's posterior with the log
# normalizer `log_a(mu)` of one child: E and Var of a, b, s, r and E[a b], with
# error allowances from the QUADPACK error estimates (the outer integrals over
# mu, the inner ones over s, and the largest relative error of A in `tracker`).
.block_power_reference <- function(setup, log_a, tracker) {

  K  <- setup$children
  I  <- K * setup$cluster$information
  va <- 1 / (1 / setup$sa^2 + 1)
  ma <- va * setup$w
  log_a_cache <- new.env()
  log_a_at <- function(mu) {
    key <- sprintf("%.17g", mu)
    if (is.null(log_a_cache[[key]])) log_a_cache[[key]] <- log_a(mu)
    log_a_cache[[key]]
  }
  inner_cache <- new.env()
  # integrals over s of U(s) N(mu; ma, va + s^2) times 1, E[a | mu, s],
  # E[a^2 | mu, s], s, s^2
  inner <- function(mu) {
    key <- sprintf("%.17g", mu)
    if (!is.null(inner_cache[[key]])) return(inner_cache[[key]])
    terms <- function(s, k) {
      total <- va + s^2
      density <- stats::dnorm(mu, ma, sqrt(total))
      mean_a  <- (ma * s^2 + mu * va) / total
      switch(k, density, density * mean_a, density * (va * s^2 / total + mean_a^2),
             density * s, density * s^2)
    }
    value <- vapply(1:5, function(k) {
      result <- stats::integrate(terms, setup$slo, setup$shi, k = k, rel.tol = 1e-11, abs.tol = 0,
                                 subdivisions = 200L)
      stopifnot(identical(result$message, "OK"))
      tracker$relative_error <- max(tracker$relative_error, result$abs.error / abs(result$value))
      result$value
    }, numeric(1L))
    inner_cache[[key]] <- value
    value
  }
  reference <- -I * 0^2 / 2 - K * log_a_at(0)
  combinations <- list(
    one = function(v, mu) v[1L], a = function(v, mu) v[2L], a2 = function(v, mu) v[3L],
    b = function(v, mu) mu * v[1L] - v[2L], b2 = function(v, mu) mu^2 * v[1L] - 2 * mu * v[2L] + v[3L],
    ab = function(v, mu) mu * v[2L] - v[3L], s = function(v, mu) v[4L], s2 = function(v, mu) v[5L])
  outer <- lapply(combinations, function(combination) {
    integrand <- function(mu) vapply(mu, function(m) {
      exp(-I * m^2 / 2 - K * log_a_at(m) - reference) * combination(inner(m), m)
    }, numeric(1L))
    result <- stats::integrate(integrand, -10, 8, rel.tol = 1e-9, abs.tol = 0, subdivisions = 400L)
    stopifnot(identical(result$message, "OK"))
    result
  })
  value <- vapply(outer, `[[`, numeric(1L), "value")
  error <- vapply(outer, `[[`, numeric(1L), "abs.error")
  moment <- value[-1L] / value[["one"]]
  # quadrature error of the ratio plus twice the relative error of the integrands
  # (A to the power K and the inner integrals) on the scale of the moment
  relative <- (K + 1) * tracker$relative_error
  scale <- c(a = sqrt(moment[["a2"]]), a2 = moment[["a2"]], b = sqrt(moment[["b2"]]), b2 = moment[["b2"]],
             ab = sqrt(moment[["a2"]] * moment[["b2"]]), s = setup$shi, s2 = setup$shi^2)
  moment_error <- (error[-1L] + abs(moment) * error[["one"]]) / value[["one"]] + 2 * relative * scale[names(moment)]
  variance <- function(first, second) {
    c(moment[[second]] - moment[[first]]^2,
      moment_error[[second]] + 2 * abs(moment[[first]]) * moment_error[[first]] + moment_error[[first]]^2)
  }
  rb <- setup$rb
  estimate <- c(a = moment[["a"]], b = moment[["b"]], s = moment[["s"]], r = .5,
                var_a = variance("a", "a2")[1L], var_b = variance("b", "b2")[1L],
                var_s = variance("s", "s2")[1L], var_r = 1 / (4 * (2 * rb + 1)), ab = moment[["ab"]])
  allowance <- c(a = moment_error[["a"]], b = moment_error[["b"]], s = moment_error[["s"]], r = 0,
                 var_a = variance("a", "a2")[2L], var_b = variance("b", "b2")[2L],
                 var_s = variance("s", "s2")[2L], var_r = 0, ab = moment_error[["ab"]])
  list(estimate = estimate, allowance = allowance,
       correlation_ab = (moment[["ab"]] - moment[["a"]] * moment[["b"]]) /
         sqrt(estimate[["var_a"]] * estimate[["var_b"]]))
}


# The observables of the power fixture from draws of a, b, s, r: their means
# estimate the reference's means, variances (centred at the draws' means) and E[a b].
.block_power_observables <- function(draws) {

  centred <- function(x) (x - mean(x))^2
  cbind(a = draws[, "a"], b = draws[, "b"], s = draws[, "s"], r = draws[, "r"],
        var_a = centred(draws[, "a"]), var_b = centred(draws[, "b"]), var_s = centred(draws[, "s"]),
        var_r = centred(draws[, "r"]), ab = draws[, "a"] * draws[, "b"])
}


# Estimates, Monte Carlo standard errors (from the effective sample size of each
# observable) and the standardized distance from a reference.
.block_compare <- function(observables, reference) {

  ess      <- as.numeric(coda::effectiveSize(coda::mcmc(observables)))
  mcse     <- sqrt(apply(observables, 2L, stats::var) / ess)
  estimate <- colMeans(observables)
  # the variance observables are centred at the draws' mean: their mean times
  # n / (n - 1) is the sample variance
  is_var <- startsWith(colnames(observables), "var_")
  estimate[is_var] <- estimate[is_var] * nrow(observables) / (nrow(observables) - 1)
  data.frame(observable = colnames(observables), estimate = unname(estimate),
             reference = unname(reference$estimate[colnames(observables)]),
             allowance = unname(reference$allowance[colnames(observables)]),
             mcse = mcse, ess = ess,
             z = unname((estimate - reference$estimate[colnames(observables)]) / mcse),
             stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------------- #
# Blocks of a and b with a fixed SD of b (the dense fixture and the stationary
# start): p(mu) ~ N(mu; ma, va + sb^2) exp(-K I mu^2 / 2) / A(mu)^K.
# ---------------------------------------------------------------------------- #

# Moments of a and b (E, Var, E[a b]) given the log of the unnormalized density
# of mu, by QUADPACK over mu.
.block_ab_reference <- function(log_mu_density, va, ma, sb, tracker, lower = -10, upper = 8) {

  v  <- va * sb^2 / (va + sb^2)
  ca <- function(mu) (ma * sb^2 + mu * va) / (va + sb^2)
  reference <- log_mu_density(0)
  functions <- list(one = function(mu) 1, a = ca, a2 = function(mu) v + ca(mu)^2,
                    b = function(mu) mu - ca(mu), b2 = function(mu) mu^2 - 2 * mu * ca(mu) + v + ca(mu)^2,
                    ab = function(mu) mu * ca(mu) - v - ca(mu)^2)
  cache <- new.env()
  density <- function(mu) vapply(mu, function(m) {
    key <- sprintf("%.17g", m)
    if (is.null(cache[[key]])) cache[[key]] <- exp(log_mu_density(m) - reference)
    cache[[key]]
  }, numeric(1L))
  outer <- lapply(functions, function(f) {
    result <- stats::integrate(function(mu) density(mu) * f(mu), lower, upper, rel.tol = 1e-10,
                               abs.tol = 0, subdivisions = 400L)
    stopifnot(identical(result$message, "OK"))
    result
  })
  value <- vapply(outer, `[[`, numeric(1L), "value")
  error <- vapply(outer, `[[`, numeric(1L), "abs.error")
  moment <- value[-1L] / value[["one"]]
  scale <- c(a = sqrt(moment[["a2"]]), a2 = moment[["a2"]], b = sqrt(moment[["b2"]]), b2 = moment[["b2"]],
             ab = sqrt(moment[["a2"]] * moment[["b2"]]))
  moment_error <- (error[-1L] + abs(moment) * error[["one"]]) / value[["one"]] +
    2 * tracker$relative_error * scale[names(moment)]
  variance <- function(first, second) {
    c(moment[[second]] - moment[[first]]^2,
      moment_error[[second]] + 2 * abs(moment[[first]]) * moment_error[[first]] + moment_error[[first]]^2)
  }
  list(estimate = c(a = moment[["a"]], b = moment[["b"]], var_a = variance("a", "a2")[1L],
                    var_b = variance("b", "b2")[1L], ab = moment[["ab"]]),
       allowance = c(a = moment_error[["a"]], b = moment_error[["b"]], var_a = variance("a", "a2")[2L],
                     var_b = variance("b", "b2")[2L], ab = moment_error[["ab"]]))
}

.block_ab_observables <- function(draws) {

  centred <- function(x) (x - mean(x))^2
  cbind(a = draws[, "a"], b = draws[, "b"], var_a = centred(draws[, "a"]), var_b = centred(draws[, "b"]),
        ab = draws[, "a"] * draws[, "b"])
}


# ---------------------------------------------------------------------------- #
# The stationary start: n independent blocks (a_i, b_i) with b_i ~ N(0, sb^2),
# one crude cluster child each, started at exact posterior draws.
# ---------------------------------------------------------------------------- #

.block_stationary_setup <- function(sb = 1.2) {

  cluster <- .coarse_exactness_setup("cluster")
  va <- 1 / (1 / 2^2 + 1)
  list(cluster = cluster, rows = cluster$rows, sa = 2, w = .5, sb = sb, va = va, ma = va * .5)
}

# log of the unnormalized density of mu = a + b with the log normalizer log_a
.block_stationary_log_mu <- function(setup, log_a) {

  function(mu) {
    stats::dnorm(mu, setup$ma, sqrt(setup$va + setup$sb^2), log = TRUE) -
      setup$cluster$information * mu^2 / 2 - log_a(mu)
  }
}

.block_stationary_reference <- function(setup, tracker) {

  .block_ab_reference(.block_stationary_log_mu(setup, function(mu) {
    .coarse_exactness_quadpack_log_normalizer(setup$cluster, mu, tracker)
  }), setup$va, setup$ma, setup$sb, tracker)
}

# Exact draws of (a, b) by rejection. The Gaussian part of the density of mu is
# N(mu; m, V); the normalizer A(mu) increases with mu (every row's probability
# of selection increases with its mean), so on a cell [t_k, t_k+1) of a grid
# 1 / A(mu) <= 1 / A(t_k), and below the grid 1 / A >= omega_low^rows bounds it.
# A cell is drawn with probability proportional to its Gaussian mass times its
# bound, mu from the truncated normal in it, and accepted with probability
# (1 / A(mu)) / bound: the accepted mu have the exact density. a | mu is normal.
.block_stationary_draws <- function(setup, n) {

  cluster  <- setup$cluster
  tracker  <- new.env()
  tracker$relative_error <- 0
  log_a    <- function(mu) .coarse_exactness_quadpack_log_normalizer(cluster, mu, tracker)
  total    <- setup$va + setup$sb^2
  precision <- 1 / total + cluster$information
  m        <- setup$ma / total / precision
  sd       <- sqrt(1 / precision)
  edges    <- m + sd * seq(-8, 8, length.out = 161L)
  bound    <- c(-cluster$rows * log(min(cluster$omega)), -vapply(edges, log_a, numeric(1L)))
  cuts     <- stats::pnorm(c(-Inf, edges, Inf), m, sd)
  log_mass <- log(diff(cuts))
  weight   <- exp(log_mass + bound - max(log_mass + bound))
  mu       <- numeric()
  proposals <- 0L
  while (length(mu) < n) {
    cell <- sample.int(length(weight), 1L, prob = weight)
    candidate <- stats::qnorm(stats::runif(1L, cuts[cell], cuts[cell + 1L]), m, sd)
    proposals <- proposals + 1L
    if (log(stats::runif(1L)) < -log_a(candidate) - bound[cell]) mu <- c(mu, candidate)
  }
  stopifnot(tracker$relative_error < 1e-6)
  v <- setup$va * setup$sb^2 / total
  a <- stats::rnorm(n, (setup$ma * setup$sb^2 + mu * setup$va) / total, sqrt(v))
  list(a = a, b = mu - a, acceptance = n / proposals)
}

# n replicate blocks started at `start` (a list with a and b).
.block_stationary_model <- function(setup, start, ...) {

  rows <- setup$rows
  n    <- length(start$a)
  syntax <- paste0("model {\n  for (i in 1:n) {\n",
    "    a[i] ~ dnorm(0, 1 / (sa * sa))\n    w[i] ~ dnorm(a[i], 1)\n    b[i] ~ dnorm(0, 1 / (sb * sb))\n",
    "    for (j in 1:", rows, ") { means[i, j] <- a[i] + b[i] }\n",
    "  ", .block_cluster_child(paste0("y[i, 1:", rows, "]"), paste0("means[i, 1:", rows, "]"), rows), "\n",
    "    marker_a[i] ~ dbt_block(a[i], i)\n    marker_b[i] ~ dbt_block(b[i], i)\n  }\n}")
  data <- c(.block_cluster_data(setup$cluster),
            list(n = n, y = matrix(0, n, rows), w = rep(setup$w, n), sa = setup$sa, sb = setup$sb,
                 marker_a = rep(0, n), marker_b = rep(0, n)))
  .block_hand_model(syntax, data, inits = list(a = start$a, b = start$b), ...)
}


# ---------------------------------------------------------------------------- #
# The dense fixture: the block (a, b) with b ~ N(0, sb^2) on one dense
# multivariate child (dselnorm_mnorm_step, .coarse_dense_setup()) with the crude
# anchor grid of that setup.
# ---------------------------------------------------------------------------- #

.block_dense_setup <- function(sb = 1.2) {

  dense <- .coarse_dense_setup()
  va <- 1 / (1 / 2^2 + 1)
  list(dense = dense, rows = dense$rows, sa = 2, w = .5, sb = sb, va = va, ma = va * .5)
}

.block_dense_model <- function(setup, inits = list(a = .2, b = -.3), ...) {

  dense <- setup$dense
  rows  <- dense$rows
  covariance <- matrix(dense$rho, rows, rows)
  diag(covariance) <- 1
  design  <- BayesTools::selection_qmc_design(6L, 512L, 8L, 9049L)
  markers <- .block_markers(c("a", "b"))
  syntax <- paste0("model {\n",
    "  a ~ dnorm(0, 1 / (sa * sa))\n  w ~ dnorm(a, 1)\n  b ~ dnorm(0, 1 / (sb * sb))\n",
    "  for (j in 1:", rows, ") { means[j] <- a + b }\n",
    "  y[1:", rows, "] ~ dselnorm_mnorm_step(means[1:", rows, "], covariance[1:", rows * (rows + 1L) / 2L,
    "], sei[1:", rows, "], omega[1:2], z_lower[1:2], z_upper[1:2], obs_bin[1:", rows, "], 1, 1, 1, ",
    "qmc[1:nqmc], 512, 8, 0.005, 0, nodes[1:nq], log_weights[1:nq], orders[1:nr])\n", markers$syntax, "\n}")
  data <- c(list(y = rep(0, rows), covariance = covariance[lower.tri(covariance, diag = TRUE)],
                 sei = rep(1, rows), omega = c(1, dense$omega), obs_bin = rep(1, rows),
                 z_lower = c(0, -1e300), z_upper = c(1e300, 0), qmc = as.numeric(design), nqmc = length(design),
                 nodes = dense$rules$nodes, log_weights = dense$rules$log_weights, orders = dense$rules$orders,
                 nq = length(dense$rules$nodes), nr = length(dense$rules$orders),
                 sa = setup$sa, w = setup$w, sb = setup$sb), markers$data)
  .block_hand_model(syntax, data, inits = inits, grid = dense$grid, ...)
}

.block_dense_reference <- function(setup, tracker) {

  .block_ab_reference(function(mu) {
    stats::dnorm(mu, setup$ma, sqrt(setup$va + setup$sb^2), log = TRUE) - setup$dense$information * mu^2 / 2 -
      .coarse_dense_log_normalizer(setup$dense, mu, setup$dense$omega, tracker)
  }, setup$va, setup$ma, setup$sb, tracker)
}


# The surrogate evaluations of a traced run, keyed by the member values (hex).
.block_surrogate_trace <- function(records) {

  records <- Filter(function(record) record$kind == 0, .block_trace_records(records))
  data.frame(key = vapply(records, function(record) paste(sprintf("%a", record$x), collapse = " "), ""),
             value = vapply(records, `[[`, numeric(1L), "value"), stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------------- #
# A block started where the surrogate is poor: the power fixture's (a, b) block
# (s and r keep their stock samplers) started at a = 2, b = -6, where the crude
# surrogate's error varies steeply (a block that stays at the start while
# adapting accepts about 6% of its proposals; one that moves, about half of its
# 500 adaptive proposals), adapted for 500 updates, then 10,000 draws. Returns
# the adaptive acceptance, the widths the block learned and the widths an exact
# sampler learns at stationarity (2.5 sqrt of the eigenvalues of the posterior
# covariance of (a, b), from the reference moments), the share of the last 200
# adaptive updates after which the block's state differs from its entry state
# (the entry states from the evaluation trace: a full evaluation followed by a
# surrogate evaluation at the same point opens a corrected update), and the
# comparison of the draws with the reference.
# ---------------------------------------------------------------------------- #

.block_poor_start <- function(min_acceptance, draws = 10000L) {

  setup   <- .block_power_setup()
  tracker <- new.env()
  tracker$relative_error <- 0
  reference <- .block_power_reference(setup, function(mu) {
    .coarse_exactness_quadpack_log_normalizer(setup$cluster, mu, tracker)
  }, tracker)
  e <- reference$estimate
  covariance <- matrix(c(e[["var_a"]], rep(e[["ab"]] - e[["a"]] * e[["b"]], 2L), e[["var_b"]]), 2L)
  hand  <- .block_power_model(setup, members = c("a", "b"), inits = list(a = 2, b = -6, s = 1.1, r = .4),
                              min_acceptance = min_acceptance)
  model <- hand$model
  RoBMA::selection_sampler_info(clear = TRUE)
  .block_test_hook(trace = TRUE)
  rjags::adapt(model, 500L, end.adaptation = FALSE)
  records <- .block_trace_records(.block_test_hook(trace = FALSE)$records)
  adaptive <- RoBMA::selection_sampler_info()
  # the entry state of every adaptive update and the state after the last
  kind  <- vapply(records, `[[`, numeric(1), "kind")
  x     <- do.call(rbind, lapply(records, `[[`, "x"))
  n     <- length(records)
  entry <- which(kind[-n] == 1 & kind[-1L] == 0 & x[-n, 1L] == x[-1L, 1L] & x[-n, 2L] == x[-1L, 2L])
  state <- model$state(internal = FALSE)[[1L]]
  after <- rbind(x[entry[-1L], , drop = FALSE], c(state[["a"]], state[["b"]]))
  late  <- utils::tail(seq_along(entry), 200L)
  moved <- mean(rowSums(after[late, , drop = FALSE] != x[entry[late], , drop = FALSE]) > 0)
  rjags::adapt(model, 0L, end.adaptation = TRUE)
  # the widths of the block built last (the hook lists the live blocks in
  # construction order: chain, d, widths)
  widths <- .block_test_hook()$widths
  values <- rjags::coda.samples(model, c("a", "b", "s", "r"), n.iter = draws, progress.bar = "none")
  list(samplers = hand$samplers, adaptive_rate = adaptive[["accepted"]] / adaptive[["proposed"]],
       adaptive_updates = length(entry), moved_late = moved, widths = sort(utils::tail(widths, 2L)),
       reference_widths = sort(2.5 * sqrt(eigen(covariance, symmetric = TRUE)$values)),
       report = RoBMA:::.selection_sampler_report_read(),
       comparison = .block_compare(.block_power_observables(as.matrix(values[[1L]])), reference))
}
