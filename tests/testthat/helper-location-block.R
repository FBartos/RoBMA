# Fixtures and independent references for the location block sampler
# (RoBMA::LocationBlock, src/samplers/LocationBlock.cc) in hand-written rjags
# models that declare the sampler's terminal tags themselves.


# The native test hook: the mutant and the probe (frozen into samplers built
# afterwards); returns the previous settings, the last factory decision and the
# probe records since the last call.
.location_hook <- function(mutant = NULL, probe = NULL) {

  .Call("RoBMA_location_block_test", if (is.null(mutant)) NULL else as.integer(mutant), probe,
        PACKAGE = "RoBMA")
}


# The native switches (LocationBlockReport.h, frozen into samplers built
# afterwards): the refresh variant (0 default, 1 R_j, 2 one elliptical slice
# step E, 3 two E steps, 4 E and R_j), the Newton fits' rule (0 the coarsest
# rung with its safeguard, 1 always the coarsest rung, 2 the sampler's rule)
# and the timing of the update's parts; negative values keep a setting.
# Returns the previous settings.
.location_config <- function(config = c(refresh = -1L, fit_rule = -1L, timing = -1L)) {

  .Call("RoBMA_location_block_test_config", as.integer(config), PACKAGE = "RoBMA")
}


# The probe records as a data frame (LocationBlockReport.h).
.location_probe_records <- function(records) {

  if (length(records) == 0L) {
    return(data.frame(type = integer(), selection = logical(), kind = numeric(), stage1 = numeric(), full = numeric(),
                      selection_part = numeric(), log_odds = numeric()))
  }
  m <- matrix(records, ncol = 6L, byrow = TRUE)
  data.frame(type = c("flip_on", "flip_off", "refresh", "study", "swap", "joint", "ellipse",
                      "family", "slope", "label_swap", "label", "restoration", "label_record")[m[, 1L] + 1L],
             selection = m[, 2L] == 1, kind = m[, 2L], stage1 = m[, 3L], full = m[, 4L], selection_part = m[, 5L],
             log_odds = m[, 6L])
}


# Batch-means comparison of observables (columns) with reference values: the
# estimate, its MCSE (50 batches), the z-score and the rule fixed beforehand
# (every |z| < 4 and at most one above 3).
.location_compare <- function(observables, reference, batches = 50L) {

  n <- nrow(observables)
  size <- floor(n / batches)
  out <- lapply(colnames(observables), function(name) {
    x <- observables[seq_len(size * batches), name]
    means <- colMeans(matrix(x, size, batches))
    estimate <- mean(x)
    mcse <- stats::sd(means) / sqrt(batches)
    data.frame(observable = name, estimate = estimate, reference = unname(reference[[name]]), mcse = mcse,
               z = (estimate - reference[[name]]) / mcse)
  })
  out <- do.call(rbind, out)
  attr(out, "passes") <- all(abs(out$z) < 4) && sum(abs(out$z) > 3) <= 1L
  out
}


# ---------------------------------------------------------------------------- #
# The Gaussian fixture: an ungated intercept, three gated dummies that are
# collinear with it (between-study, one study each, the fourth study the
# reference), standardized study effects with a fixed scale, known sampling
# variances. Every child is a plain normal: the sampler's Gaussian label.
# ---------------------------------------------------------------------------- #

.location_gaussian_setup <- function() {

  studies <- rep(1:4, each = 3L)
  x <- sapply(1:3, function(j) as.numeric(studies == j))
  list(
    y        = c(0.52, 0.31, 0.44, 0.71, 0.88, 0.63, 0.05, 0.22, 0.17, 0.35, 0.12, 0.29),
    variance = rep(c(0.04, 0.06, 0.05), 4L),
    study    = studies,
    x        = x,
    scale    = 0.3,
    prior_intercept = 1,  # precisions
    prior_slab      = 4,
    weights  = c(1, 1)    # gate prior weights (spike, slab)
  )
}


.location_gaussian_model <- function(setup) {

  n <- length(setup$y)
  S <- max(setup$study)
  syntax <- c(
    "model {",
    "  b0 ~ dnorm(0, prior_intercept)",
    "  for (j in 1:3) {",
    "    g[j] ~ dcat(gate_weights)",
    "    slab[j] ~ dnorm(0, prior_slab)",
    "    b[j] <- 0 * (g[j] == 1) + slab[j] * (g[j] == 2)",
    "  }",
    "  for (s in 1:S) {",
    "    Z[s] ~ dnorm(0, 1)",
    "    u[s] <- scale * Z[s]",
    "  }",
    "  for (i in 1:n) {",
    "    mu[i] <- b0 + b[1] * x[i, 1] + b[2] * x[i, 2] + b[3] * x[i, 3] + u[study[i]]",
    "    y[i] ~ dnorm(mu[i], 1 / variance[i])",
    "  }",
    "  tag_group_1 <- robma_location_group(b0, b0, 0, b0, 1)",
    paste0("  tag_group_", 2:4, " <- robma_location_group(g[", 1:3, "], slab[", 1:3, "], 2, b[", 1:3, "], ", 2:4, ")"),
    paste0("  tag_latent_", 1:S, " <- robma_location_latent(Z[", 1:S, "], ", 1:S, ")"),
    paste0("  tag_setup <- robma_location_setup(b0, 4, ", S, ")"),
    "}")
  data <- list(y = setup$y, variance = setup$variance, study = setup$study, x = setup$x, n = n, S = S,
               scale = setup$scale, prior_intercept = setup$prior_intercept, prior_slab = setup$prior_slab,
               gate_weights = setup$weights)
  list(syntax = paste(syntax, collapse = "\n"), data = data)
}


# Closed-form posterior: every gate pattern's marginal likelihood
# N(y; 0, Sigma + scale^2 H H' + X_a V_a X_a') and the posterior means of the
# physical coefficients given the pattern.
.location_gaussian_reference <- function(setup) {

  n <- length(setup$y)
  H <- outer(setup$study, seq_len(max(setup$study)), `==`) * 1
  S <- diag(setup$variance) + setup$scale^2 * H %*% t(H)
  patterns <- as.matrix(expand.grid(g1 = 0:1, g2 = 0:1, g3 = 0:1))
  log_prior <- log(setup$weights / sum(setup$weights))
  rows <- lapply(seq_len(nrow(patterns)), function(k) {
    on <- which(patterns[k, ] == 1)
    X <- cbind(1, setup$x[, on, drop = FALSE])
    V <- diag(c(1 / setup$prior_intercept, rep(1 / setup$prior_slab, length(on))), nrow = 1L + length(on))
    C <- S + X %*% V %*% t(X)
    R <- chol(C)
    log_ml <- -sum(log(diag(R))) - 0.5 * n * log(2 * pi) - 0.5 * sum(backsolve(R, setup$y, transpose = TRUE)^2)
    mean <- V %*% t(X) %*% solve(C, setup$y)
    b <- numeric(4)
    b[c(1, 1 + on)] <- mean
    c(log_weight = log_ml + sum(log_prior[patterns[k, ] + 1]), b = b)
  })
  rows <- do.call(rbind, rows)
  w <- exp(rows[, "log_weight"] - max(rows[, "log_weight"]))
  w <- w / sum(w)
  c(g1 = sum(w * patterns[, 1]), g2 = sum(w * patterns[, 2]), g3 = sum(w * patterns[, 3]),
    b0 = sum(w * rows[, "b1"]), b1 = sum(w * rows[, "b2"]), b2 = sum(w * rows[, "b3"]), b3 = sum(w * rows[, "b4"]))
}


# Compiles a hand-written location model (one chain, the given mutant frozen
# into its sampler), runs adapt + burnin + draws and returns the draws, the
# sampler of each node, the factory decision and the per-fit counters.
.location_hand_chain <- function(model, monitor, draws, burnin = 500L, adapt = 100L, mutant = 0L,
                                 probe = FALSE, seed = 4127L, sampler = "slice.coarse", min_acceptance = 0,
                                 config = NULL) {

  BayesTools::BayesTools_load_JAGS_module(quiet = TRUE)
  old <- RoBMA::RoBMA.options()[c("selection.sampler", "selection.coarse_max_rules",
                                  "selection.coarse_min_acceptance")]
  on.exit(do.call(RoBMA::RoBMA.options, old), add = TRUE)
  RoBMA::RoBMA.options(selection.sampler = sampler, selection.coarse_max_rules = 3L,
                       selection.coarse_min_acceptance = min_acceptance)
  .location_hook(mutant = mutant, probe = probe)
  on.exit(.location_hook(mutant = 0L, probe = FALSE), add = TRUE)
  if (!is.null(config)) {
    previous <- .location_config(config)
    on.exit(.location_config(previous), add = TRUE)
  }
  RoBMA:::.location_block_report_read()
  connection <- textConnection(model$syntax)
  on.exit(close(connection), add = TRUE)
  jags <- rjags::jags.model(connection, data = model$data, n.chains = 1L, n.adapt = 0L, quiet = TRUE,
                            inits = list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed))
  hook <- .location_hook()
  samplers <- rjags::list.samplers(jags)
  assigned <- unlist(lapply(seq_along(samplers), function(i) {
    stats::setNames(rep(names(samplers)[i], length(samplers[[i]])), samplers[[i]])
  }))
  if (adapt > 0L) rjags::adapt(jags, adapt, end.adaptation = FALSE)
  rjags::adapt(jags, 0L, end.adaptation = TRUE)
  if (burnin > 0L) stats::update(jags, n.iter = burnin, progress.bar = "none")
  RoBMA:::.location_block_report_read()
  values <- rjags::coda.samples(jags, monitor, n.iter = draws, progress.bar = "none")
  list(draws = as.matrix(values[[1L]]), samplers = assigned, decision = hook$decision,
       counters = RoBMA:::.location_block_report_read(), probe = .location_probe_records(.location_hook()$records))
}


# Compiles a hand-written location model (one chain, no adaptation run) with
# the mutant and switches frozen into its sampler and the selection options of
# .location_hand_chain.
.location_compile <- function(model, inits = list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 4127L),
                              mutant = 0L, config = NULL, min_acceptance = 0) {

  BayesTools::BayesTools_load_JAGS_module(quiet = TRUE)
  old <- RoBMA::RoBMA.options()[c("selection.sampler", "selection.coarse_max_rules",
                                  "selection.coarse_min_acceptance")]
  on.exit(do.call(RoBMA::RoBMA.options, old), add = TRUE)
  RoBMA::RoBMA.options(selection.sampler = "slice.coarse", selection.coarse_max_rules = 3L,
                       selection.coarse_min_acceptance = min_acceptance)
  .location_hook(mutant = mutant)
  on.exit(.location_hook(mutant = 0L), add = TRUE)
  if (!is.null(config)) {
    previous <- .location_config(config)
    on.exit(.location_config(previous), add = TRUE)
  }
  connection <- textConnection(model$syntax)
  on.exit(close(connection), add = TRUE)
  rjags::jags.model(connection, data = model$data, n.chains = 1L, n.adapt = 0L, quiet = TRUE, inits = inits)
}


.location_gaussian_observables <- function(draws) {

  cbind(g1 = draws[, "g[1]"] == 2, g2 = draws[, "g[2]"] == 2, g3 = draws[, "g[3]"] == 2,
        b0 = draws[, "b0"], b1 = draws[, "b[1]"], b2 = draws[, "b[2]"], b3 = draws[, "b[3]"]) * 1
}


# Short single-chain fits: muffle only the R-hat warning of one chain.
.location_one_chain <- function(expr) {

  withCallingHandlers(expr, warning = function(w) {
    if (grepl("Only one chain was run", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
}


# Data of the public location block fits: 8 studies of 3 estimates with
# independent sampling errors (several likelihood blocks share each study
# effect), a three-level factor, and a numeric moderator with overwhelming
# evidence; no heterogeneity.
.location_test_data <- function() {

  data.frame(
    study = rep(1:8, each = 3L),
    group = factor(rep(c("a", "b", "c"), 8L)),
    x     = rep(c(-1, 0, 1), 8L) + rep(seq(-0.2, 0.2, length.out = 8L), each = 3L),
    yi    = 0.15 + 0.8 * (rep(c(-1, 0, 1), 8L) + rep(seq(-0.2, 0.2, length.out = 8L), each = 3L)) +
      rep(c(0.03, -0.02, 0.05, -0.04, 0.01, 0.02, -0.03, 0.04), 3L),
    vi    = rep(c(0.002, 0.003, 0.0025), 8L)
  )
}


# ---------------------------------------------------------------------------- #
# The selection fixture: three clusters of 4, 3 and 3 estimates with
# strongly correlated sampling errors (rank-one cluster selection blocks,
# a crude coarse schedule whose rule 3 has three nodes), a gated intercept with
# a nonzero slab mean, a gated between-cluster moderator, a gated moderator
# collinear with the intercept, standardized cluster effects with a fixed
# scale and a fixed one-sided step selection (weight omega below zero).
# ---------------------------------------------------------------------------- #

.location_selection_setup <- function(omega = .05, rho = .85) {

  se <- c(.20, .22, .25, .21, .18, .24, .19, .23, .20, .22)
  list(
    y           = c(.25, .32, .18, .40, .45, .30, .52, .12, .28, .22),
    se          = se,
    residual_sd = se * sqrt(1 - rho),
    loading     = se * sqrt(rho),
    cluster     = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L),
    x1          = c(0, 0, 0, 0, 1, 1, 1, .5, .5, .5),
    x2          = c(1, 1, 1, 1, 1, 1, 1, 0, 1, 1),
    scale       = .15,
    prior_mean  = c(.2, 0, 0),
    prior_precision = c(4, 4, 4),
    weights     = c(1, 1),
    omega       = c(1, omega),
    rules       = .coarse_rules(c(1L, 2L, 3L, 7L, 15L, 31L, 63L, 127L, 255L, 511L, 1023L))
  )
}


.location_selection_model <- function(setup) {

  clusters <- split(seq_along(setup$y), setup$cluster)
  child <- function(s) {
    n <- length(clusters[[s]])
    paste0("  y", s, "[1:", n, "] ~ dselnorm_cluster_step(m", s, "[1:", n, "], rsd", s, "[1:", n, "], ld", s,
           "[1:", n, "], se", s, "[1:", n, "], omega[1:2], z_lower[1:2], z_upper[1:2], bin", s, "[1:", n,
           "], 1, 1, 1, nodes[1:nq], log_weights[1:nq], orders[1:nr], qmc[1:nqmc], 64, 512, 8, 0.005, 0)")
  }
  means <- unlist(lapply(seq_along(clusters), function(s) {
    rows <- clusters[[s]]
    paste0("  m", s, "[", seq_along(rows), "] <- b[1] + b[2] * x1[", rows, "] + b[3] * x2[", rows,
           "] + scale * Z[", s, "]")
  }))
  syntax <- c(
    "model {",
    "  for (j in 1:3) {",
    "    g[j] ~ dcat(gate_weights)",
    "    slab[j] ~ dnorm(prior_mean[j], prior_precision[j])",
    "    b[j] <- 0 * (g[j] == 1) + slab[j] * (g[j] == 2)",
    "  }",
    "  for (s in 1:3) { Z[s] ~ dnorm(0, 1) }",
    means,
    vapply(seq_along(clusters), child, character(1)),
    paste0("  tag_group_", 1:3, " <- robma_location_group(g[", 1:3, "], slab[", 1:3, "], 2, b[", 1:3, "], ", 1:3, ")"),
    paste0("  tag_latent_", 1:3, " <- robma_location_latent(Z[", 1:3, "], ", 1:3, ")"),
    "  tag_setup <- robma_location_setup(g[1], 3, 3)",
    "}")
  design <- BayesTools::selection_qmc_design(2L, 512L, 8L, 9049L)
  data <- list(x1 = setup$x1, x2 = setup$x2, scale = setup$scale, prior_mean = setup$prior_mean,
               prior_precision = setup$prior_precision, gate_weights = setup$weights, omega = setup$omega,
               z_lower = c(0, -1e300), z_upper = c(1e300, 0), qmc = as.numeric(design), nqmc = length(design),
               nodes = setup$rules$nodes, log_weights = setup$rules$log_weights, orders = setup$rules$orders,
               nq = length(setup$rules$nodes), nr = length(setup$rules$orders))
  for (s in seq_along(clusters)) {
    rows <- clusters[[s]]
    data[[paste0("y", s)]]   <- setup$y[rows]
    data[[paste0("rsd", s)]] <- setup$residual_sd[rows]
    data[[paste0("ld", s)]]  <- setup$loading[rows]
    data[[paste0("se", s)]]  <- setup$se[rows]
    data[[paste0("bin", s)]] <- ifelse(setup$y[rows] > 0, 1L, 2L)
  }
  list(syntax = paste(syntax, collapse = "\n"), data = data, qmc = design)
}


# The exact log likelihood of cluster s at the rows of `means` (one row per
# point), by the package's batch cluster kernel with the model's controls.
.location_selection_cluster_loglik <- function(setup, model, s, means) {

  rows <- which(setup$cluster == s)
  n <- nrow(means)
  out <- .Call("RoBMA_selnorm_cluster_step_loglik_batch", setup$y[rows], means,
               matrix(rep(setup$residual_sd[rows], each = n), n), matrix(rep(setup$loading[rows], each = n), n),
               setup$se[rows], matrix(rep(setup$omega, each = n), n), c(0, -Inf), c(Inf, 0),
               ifelse(setup$y[rows] > 0, 1L, 2L), 1L, TRUE, 1L, setup$rules$nodes, setup$rules$log_weights,
               as.double(setup$rules$orders), as.double(model$qmc), 64L, 512L, 8L, .005, FALSE, 0L,
               PACKAGE = "RoBMA")
  out$log_density
}


# The reference: per gate pattern, the active coefficients integrated on a
# Gauss-Hermite product grid adapted to the mode and curvature of their
# marginal density, the cluster effects by Gauss-Hermite rules against their
# N(0, 1) priors; the gate probabilities, the first two moments of the
# physical coefficients, and the mean of the first cluster's effect Z[1] and
# of its product with the intercept (a joint latent and coefficient moment).
.location_selection_reference <- function(setup, model, coefficient_nodes = 16L, latent_nodes = 32L) {

  latent <- .coarse_reference_rule(latent_nodes)
  X <- cbind(1, setup$x1, setup$x2)
  # log of the likelihood integrated over the cluster effects, at coefficient
  # points (rows of B, three columns), with the conditional mean of Z[1]
  log_marginal <- function(B) {
    total <- numeric(nrow(B))
    for (s in 1:3) {
      rows  <- which(setup$cluster == s)
      base  <- B %*% t(X[rows, , drop = FALSE])
      means <- base[rep(seq_len(nrow(B)), each = latent_nodes), , drop = FALSE] +
        setup$scale * latent$nodes[rep(seq_len(latent_nodes), nrow(B))]
      ll  <- matrix(.location_selection_cluster_loglik(setup, model, s, means), nrow = latent_nodes)
      top <- apply(ll, 2L, max)
      weights <- latent$weights * exp(sweep(ll, 2L, top))
      total <- total + top + log(colSums(weights))
      if (s == 1L) z1 <- colSums(latent$nodes * weights) / colSums(weights)
    }
    attr(total, "z1") <- z1
    total
  }
  patterns <- as.matrix(expand.grid(g1 = 0:1, g2 = 0:1, g3 = 0:1))
  log_gate <- log(setup$weights / sum(setup$weights))
  results <- lapply(seq_len(nrow(patterns)), function(k) {
    on <- which(patterns[k, ] == 1)
    log_prior_gate <- sum(log_gate[patterns[k, ] + 1])
    if (length(on) == 0L) {
      none <- log_marginal(matrix(0, 1, 3))
      return(list(log_mass = log_prior_gate + as.numeric(none), m1 = rep(0, 3), m2 = rep(0, 3),
                  z1 = attr(none, "z1"), b1z1 = 0))
    }
    full <- function(beta) {
      beta <- matrix(beta, ncol = length(on))
      B <- matrix(0, nrow(beta), 3)
      B[, on] <- beta
      prior <- rep(0, nrow(beta))
      for (i in seq_along(on)) {
        prior <- prior + stats::dnorm(beta[, i], setup$prior_mean[on[i]], 1 / sqrt(setup$prior_precision[on[i]]),
                                      log = TRUE)
      }
      value <- log_marginal(B)
      out <- as.numeric(value) + prior
      attr(out, "z1") <- attr(value, "z1")
      out
    }
    fit <- stats::optim(setup$prior_mean[on], function(beta) -full(beta), method = "BFGS", hessian = TRUE)
    L <- t(chol(solve(fit$hessian)))
    rule <- .coarse_reference_rule(coefficient_nodes)
    grid <- as.matrix(expand.grid(rep(list(seq_len(coefficient_nodes)), length(on))))
    x <- matrix(rule$nodes[grid], ncol = length(on))
    w <- apply(matrix(rule$weights[grid], ncol = length(on)), 1L, prod)
    beta <- sweep(x %*% t(L), 2L, fit$par, `+`)
    # the N(0, I) weight of the rule is divided out
    value <- full(beta)
    log_f <- as.numeric(value) + 0.5 * rowSums(x^2)
    top <- max(log_f)
    weight <- w * exp(log_f - top)
    mass <- sum(weight)
    m1 <- m2 <- rep(0, 3)
    m1[on] <- colSums(weight * beta) / mass
    m2[on] <- colSums(weight * beta^2) / mass
    intercept <- if (1L %in% on) beta[, match(1L, on)] else 0
    list(log_mass = log_prior_gate + top + log(mass) + sum(log(diag(L))) + 0.5 * length(on) * log(2 * pi),
         m1 = m1, m2 = m2, z1 = sum(weight * attr(value, "z1")) / mass,
         b1z1 = sum(weight * intercept * attr(value, "z1")) / mass)
  })
  log_mass <- vapply(results, `[[`, numeric(1), "log_mass")
  p  <- exp(log_mass - max(log_mass))
  p  <- p / sum(p)
  m1 <- colSums(p * t(vapply(results, `[[`, numeric(3), "m1")))
  m2 <- colSums(p * t(vapply(results, `[[`, numeric(3), "m2")))
  c(g1 = sum(p * patterns[, 1]), g2 = sum(p * patterns[, 2]), g3 = sum(p * patterns[, 3]),
    b1 = m1[[1]], b2 = m1[[2]], b3 = m1[[3]], b1_2 = m2[[1]], b2_2 = m2[[2]], b3_2 = m2[[3]],
    z1 = sum(p * vapply(results, `[[`, numeric(1), "z1")), b1z1 = sum(p * vapply(results, `[[`, numeric(1), "b1z1")))
}


# The observables of the selection fixture; with the cluster effects
# monitored also Z[1] and its product with the intercept.
.location_selection_observables <- function(draws) {

  b <- draws[, c("b[1]", "b[2]", "b[3]")]
  out <- cbind(g1 = draws[, "g[1]"] == 2, g2 = draws[, "g[2]"] == 2, g3 = draws[, "g[3]"] == 2,
               b1 = b[, 1], b2 = b[, 2], b3 = b[, 3], b1_2 = b[, 1]^2, b2_2 = b[, 2]^2, b3_2 = b[, 3]^2) * 1
  if ("Z[1]" %in% colnames(draws)) out <- cbind(out, z1 = draws[, "Z[1]"], b1z1 = b[, 1] * draws[, "Z[1]"])
  out
}


# ---------------------------------------------------------------------------- #
# The step-switch fixture: 14 independent dselnorm_step_switch rows (the
# children of RoBMA() meta-regressions), a one-sided step selection, row
# observation weights, one bias-exempt row (bin 0), three gates (one collinear
# with the intercept); `direction` -1 mirrors the effects, the slab means and
# the selection (effect sign -1), so that the posterior of the physical
# coefficients is mirrored. Its reference writes the row density in R.
# ---------------------------------------------------------------------------- #

.location_step_setup <- function(direction = 1) {

  se <- c(.20, .22, .25, .21, .18, .24, .19, .23, .20, .22, .26, .17, .21, .24)
  y  <- c(.45, .52, .10, .61, .70, .48, .05, .44, .58, .02, .66, .49, .30, .75)
  bin <- ifelse(y / se > 1.96, 1L, 2L)
  bin[13] <- 0L
  list(direction = direction, y = direction * y, se = se, sd = sqrt(se^2 + .1^2),
       x1 = c(-1, -.6, -.2, .2, .6, 1, -1, -.5, 0, .5, 1, -.8, .3, .9),
       x2 = c(1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1),
       power = c(1, 1, .5, 1, 1.5, 1, 1, 1, 2, 1, 1, .7, 1, 1),
       omega = c(1, .3), z_lower = c(1.96, -1e300), z_upper = c(1e300, 1.96), bin = bin,
       prior_mean = direction * c(.2, 0, 0))
}


.location_step_model <- function(setup) {

  syntax <- c(
    "model {",
    "  for (j in 1:3) {",
    "    g[j] ~ dcat(gate_weights)",
    "    slab[j] ~ dnorm(prior_mean[j], 4)",
    "    b[j] <- 0 * (g[j] == 1) + slab[j] * (g[j] == 2)",
    "  }",
    "  for (i in 1:n) {",
    "    mu[i] <- b[1] + b[2] * x1[i] + b[3] * x2[i]",
    "    y[i] ~ dselnorm_step_switch(mu[i], sd[i], se[i], power[i], omega[1:2], z_lower[1:2], z_upper[1:2], bin[i], direction, 1, 1)",
    "  }",
    paste0("  tag_group_", 1:3, " <- robma_location_group(g[", 1:3, "], slab[", 1:3, "], 2, b[", 1:3, "], ", 1:3, ")"),
    "  tag_setup <- robma_location_setup(g[1], 3, 0)",
    "}")
  data <- setup[c("y", "sd", "se", "power", "omega", "z_lower", "z_upper", "bin", "x1", "x2", "prior_mean",
                  "direction")]
  data$n <- length(setup$y)
  data$gate_weights <- c(1, 1)
  list(syntax = paste(syntax, collapse = "\n"), data = data)
}


# The row log density written independently: power * [log N(y; mu, sd) +
# log omega_bin - log A(mu)], A(mu) = sum_b omega_b P(lower_b < s Y / se <
# upper_b), Y ~ N(mu, sd^2), s the effect direction; an exempt row power *
# log N. The reference as .location_selection_reference without latents.
.location_step_reference <- function(setup, nodes = 30L) {

  s  <- setup$direction
  lo <- ifelse(setup$z_lower <= -1e299, -Inf, setup$z_lower)
  up <- ifelse(setup$z_upper >= 1e299, Inf, setup$z_upper)
  loglik <- function(B) {
    mu <- B %*% t(cbind(1, setup$x1, setup$x2))
    total <- numeric(nrow(B))
    for (i in seq_along(setup$y)) {
      m <- mu[, i]
      v <- stats::dnorm(setup$y[i], m, setup$sd[i], log = TRUE)
      if (setup$bin[i] > 0) {
        A <- 0
        for (b in seq_along(setup$omega)) {
          # P(lo < s Y / se < up) for the sign s
          bounds <- sort(s * c(lo[b], up[b]) * setup$se[i])
          A <- A + setup$omega[b] * (stats::pnorm(bounds[2], m, setup$sd[i]) - stats::pnorm(bounds[1], m, setup$sd[i]))
        }
        v <- v + log(setup$omega[setup$bin[i]]) - log(A)
      }
      total <- total + setup$power[i] * v
    }
    total
  }
  patterns <- as.matrix(expand.grid(g1 = 0:1, g2 = 0:1, g3 = 0:1))
  results <- lapply(seq_len(nrow(patterns)), function(k) {
    on <- which(patterns[k, ] == 1)
    lp <- 3 * log(0.5)
    if (length(on) == 0L) return(list(log_mass = lp + loglik(matrix(0, 1, 3)), m1 = rep(0, 3), m2 = rep(0, 3)))
    full <- function(beta) {
      beta <- matrix(beta, ncol = length(on))
      B <- matrix(0, nrow(beta), 3)
      B[, on] <- beta
      prior <- rep(0, nrow(beta))
      for (i in seq_along(on)) prior <- prior + stats::dnorm(beta[, i], setup$prior_mean[on[i]], 0.5, log = TRUE)
      loglik(B) + prior
    }
    fit <- stats::optim(setup$prior_mean[on], function(beta) -full(beta), method = "BFGS", hessian = TRUE)
    L <- t(chol(solve(fit$hessian)))
    rule <- .coarse_reference_rule(nodes)
    grid <- as.matrix(expand.grid(rep(list(seq_len(nodes)), length(on))))
    x <- matrix(rule$nodes[grid], ncol = length(on))
    w <- apply(matrix(rule$weights[grid], ncol = length(on)), 1L, prod)
    beta <- sweep(x %*% t(L), 2L, fit$par, `+`)
    log_f <- full(beta) + 0.5 * rowSums(x^2)
    top <- max(log_f)
    weight <- w * exp(log_f - top)
    mass <- sum(weight)
    m1 <- m2 <- rep(0, 3)
    m1[on] <- colSums(weight * beta) / mass
    m2[on] <- colSums(weight * beta^2) / mass
    list(log_mass = lp + top + log(mass) + sum(log(diag(L))) + 0.5 * length(on) * log(2 * pi), m1 = m1, m2 = m2)
  })
  log_mass <- vapply(results, `[[`, numeric(1), "log_mass")
  p <- exp(log_mass - max(log_mass))
  p <- p / sum(p)
  M1 <- colSums(p * t(vapply(results, `[[`, numeric(3), "m1")))
  M2 <- colSums(p * t(vapply(results, `[[`, numeric(3), "m2")))
  c(g1 = sum(p * patterns[, 1]), g2 = sum(p * patterns[, 2]), g3 = sum(p * patterns[, 3]),
    b1 = M1[[1]], b2 = M1[[2]], b3 = M1[[3]], b1_2 = M2[[1]], b2_2 = M2[[2]], b3_2 = M2[[3]])
}


# ---------------------------------------------------------------------------- #
# The two-coordinate slab fixture: the selection fixture's clusters with a
# gated intercept and one gated group of two coordinates (dmnorm slab with a
# correlated prior precision) for the two moderators.
# ---------------------------------------------------------------------------- #

.location_factor_model <- function(setup, slab_mean = c(.1, -.1), slab_precision = matrix(c(4, 1.5, 1.5, 4), 2L)) {

  clusters <- split(seq_along(setup$y), setup$cluster)
  child <- function(s) {
    n <- length(clusters[[s]])
    paste0("  y", s, "[1:", n, "] ~ dselnorm_cluster_step(m", s, "[1:", n, "], rsd", s, "[1:", n, "], ld", s,
           "[1:", n, "], se", s, "[1:", n, "], omega[1:2], z_lower[1:2], z_upper[1:2], bin", s, "[1:", n,
           "], 1, 1, 1, nodes[1:nq], log_weights[1:nq], orders[1:nr], qmc[1:nqmc], 64, 512, 8, 0.005, 0)")
  }
  means <- unlist(lapply(seq_along(clusters), function(s) {
    rows <- clusters[[s]]
    paste0("  m", s, "[", seq_along(rows), "] <- b1 + b2[1] * x1[", rows, "] + b2[2] * x2[", rows,
           "] + scale * Z[", s, "]")
  }))
  syntax <- c(
    "model {",
    "  g1 ~ dcat(gate_weights)",
    "  slab1 ~ dnorm(0.2, 4)",
    "  b1 <- 0 * (g1 == 1) + slab1 * (g1 == 2)",
    "  g2 ~ dcat(gate_weights)",
    "  slab2[1:2] ~ dmnorm(slab_mean[1:2], slab_precision[1:2, 1:2])",
    "  b2[1:2] <- zero2[1:2] * (g2 == 1) + slab2[1:2] * (g2 == 2)",
    "  for (s in 1:3) { Z[s] ~ dnorm(0, 1) }",
    means,
    vapply(seq_along(clusters), child, character(1)),
    "  tag_group_1 <- robma_location_group(g1, slab1, 2, b1, 1)",
    "  tag_group_2 <- robma_location_group(g2, slab2[1:2], 2, b2[1:2], 2)",
    paste0("  tag_latent_", 1:3, " <- robma_location_latent(Z[", 1:3, "], ", 1:3, ")"),
    "  tag_setup <- robma_location_setup(g1, 2, 3)",
    "}")
  design <- BayesTools::selection_qmc_design(2L, 512L, 8L, 9049L)
  data <- list(x1 = setup$x1, x2 = setup$x2, scale = setup$scale, slab_mean = slab_mean,
               slab_precision = slab_precision, zero2 = c(0, 0), gate_weights = setup$weights, omega = setup$omega,
               z_lower = c(0, -1e300), z_upper = c(1e300, 0), qmc = as.numeric(design), nqmc = length(design),
               nodes = setup$rules$nodes, log_weights = setup$rules$log_weights, orders = setup$rules$orders,
               nq = length(setup$rules$nodes), nr = length(setup$rules$orders))
  for (s in seq_along(clusters)) {
    rows <- clusters[[s]]
    data[[paste0("y", s)]]   <- setup$y[rows]
    data[[paste0("rsd", s)]] <- setup$residual_sd[rows]
    data[[paste0("ld", s)]]  <- setup$loading[rows]
    data[[paste0("se", s)]]  <- setup$se[rows]
    data[[paste0("bin", s)]] <- ifelse(setup$y[rows] > 0, 1L, 2L)
  }
  list(syntax = paste(syntax, collapse = "\n"), data = data, qmc = design, slab_mean = slab_mean,
       slab_precision = slab_precision)
}


# Per gate pattern the active physical coefficients (b1; b2[1:2]) on a
# Laplace-adapted Gauss-Hermite grid, the cluster effects by Gauss-Hermite
# against N(0, 1), the package's exact cluster likelihood.
.location_factor_reference <- function(setup, model, coefficient_nodes = 16L, latent_nodes = 32L) {

  latent <- .coarse_reference_rule(latent_nodes)
  X <- cbind(1, setup$x1, setup$x2)
  m2 <- model$slab_mean
  P2 <- model$slab_precision
  log_marginal <- function(B) {
    total <- numeric(nrow(B))
    for (s in 1:3) {
      rows <- which(setup$cluster == s)
      base <- B %*% t(X[rows, , drop = FALSE])
      mm <- base[rep(seq_len(nrow(B)), each = latent_nodes), , drop = FALSE] +
        setup$scale * latent$nodes[rep(seq_len(latent_nodes), nrow(B))]
      ll  <- matrix(.location_selection_cluster_loglik(setup, model, s, mm), nrow = latent_nodes)
      top <- apply(ll, 2L, max)
      total <- total + top + log(colSums(latent$weights * exp(sweep(ll, 2L, top))))
    }
    total
  }
  log_mvn <- function(x, m, P) {
    x <- sweep(matrix(x, ncol = 2L), 2L, m)
    0.5 * as.numeric(determinant(P)$modulus) - log(2 * pi) - 0.5 * rowSums((x %*% P) * x)
  }
  patterns <- as.matrix(expand.grid(g1 = 0:1, g2 = 0:1))
  log_gate <- log(setup$weights / sum(setup$weights))
  results <- lapply(seq_len(nrow(patterns)), function(k) {
    on <- c(if (patterns[k, 1] == 1) 1L, if (patterns[k, 2] == 1) 2:3)
    log_prior_gate <- sum(log_gate[patterns[k, ] + 1])
    if (length(on) == 0L) {
      return(list(log_mass = log_prior_gate + log_marginal(matrix(0, 1, 3)), m1 = rep(0, 3), m2 = rep(0, 3), m23 = 0))
    }
    full <- function(beta) {
      beta <- matrix(beta, ncol = length(on))
      B <- matrix(0, nrow(beta), 3)
      B[, on] <- beta
      prior <- rep(0, nrow(beta))
      if (1L %in% on) prior <- prior + stats::dnorm(B[, 1], 0.2, 0.5, log = TRUE)
      if (2L %in% on) prior <- prior + log_mvn(B[, 2:3, drop = FALSE], m2, P2)
      log_marginal(B) + prior
    }
    start <- c(if (1L %in% on) 0.2, if (2L %in% on) m2)
    fit <- stats::optim(start, function(beta) -full(beta), method = "BFGS", hessian = TRUE)
    L <- t(chol(solve(fit$hessian)))
    rule <- .coarse_reference_rule(coefficient_nodes)
    grid <- as.matrix(expand.grid(rep(list(seq_len(coefficient_nodes)), length(on))))
    x <- matrix(rule$nodes[grid], ncol = length(on))
    w <- apply(matrix(rule$weights[grid], ncol = length(on)), 1L, prod)
    beta <- sweep(x %*% t(L), 2L, fit$par, `+`)
    log_f <- full(beta) + 0.5 * rowSums(x^2)
    top <- max(log_f)
    weight <- w * exp(log_f - top)
    mass <- sum(weight)
    m1 <- mm2 <- rep(0, 3)
    m1[on] <- colSums(weight * beta) / mass
    mm2[on] <- colSums(weight * beta^2) / mass
    m23 <- if (2L %in% on) sum(weight * beta[, match(2L, on)] * beta[, match(3L, on)]) / mass else 0
    list(log_mass = log_prior_gate + top + log(mass) + sum(log(diag(L))) + 0.5 * length(on) * log(2 * pi),
         m1 = m1, m2 = mm2, m23 = m23)
  })
  log_mass <- vapply(results, `[[`, numeric(1), "log_mass")
  p  <- exp(log_mass - max(log_mass))
  p  <- p / sum(p)
  M1 <- colSums(p * t(vapply(results, `[[`, numeric(3), "m1")))
  M2 <- colSums(p * t(vapply(results, `[[`, numeric(3), "m2")))
  c(g1 = sum(p * patterns[, 1]), g2 = sum(p * patterns[, 2]), b1 = M1[[1]], b21 = M1[[2]], b22 = M1[[3]],
    b1_2 = M2[[1]], b21_2 = M2[[2]], b22_2 = M2[[3]], b21b22 = sum(p * vapply(results, `[[`, numeric(1), "m23")))
}


.location_factor_observables <- function(draws) {

  cbind(g1 = draws[, "g1"] == 2, g2 = draws[, "g2"] == 2, b1 = draws[, "b1"], b21 = draws[, "b2[1]"],
        b22 = draws[, "b2[2]"], b1_2 = draws[, "b1"]^2, b21_2 = draws[, "b2[1]"]^2, b22_2 = draws[, "b2[2]"]^2,
        b21b22 = draws[, "b2[1]"] * draws[, "b2[2]"]) * 1
}


# The sweep sequences of the trace switch (bit 2; LocationBlockReport.h): a
# list of sweeps, each with its chain, direction (reversed), entry label
# (selection), the slots in the order they ran (kind, a, b) and the
# coordinates of composite slice refreshes (group, coordinate) in that order.
.location_sequences <- function(records = .Call("RoBMA_location_block_test_sequence", PACKAGE = "RoBMA")) {

  out <- list()
  i <- 1L
  while (i <= length(records)) {
    stopifnot(records[[i]] == -1)
    n <- records[[i + 4L]]
    entries <- matrix(records[i + 4L + seq_len(3L * n)], ncol = 3L, byrow = TRUE)
    coordinate <- entries[, 1L] == 100
    out[[length(out) + 1L]] <- list(chain = records[[i + 1L]], reversed = records[[i + 2L]] == 1,
                                    selection = records[[i + 3L]] == 1,
                                    slots = entries[!coordinate, , drop = FALSE],
                                    coordinates = entries[coordinate, 2:3, drop = FALSE])
    i <- i + 5L + 3L * n
  }
  out
}


# ---------------------------------------------------------------------------- #
# The dense fixture: three dselnorm_mnorm_step blocks of three equicorrelated
# rows (unit variances and SEs, correlation rho; the children of known-V fits
# that the automatic policy gives the block), a one-sided step selection with
# weight omega below zero, a gated intercept and a gated moderator that varies
# within the blocks. Its reference writes the block density in R: the
# multivariate normal density times the observed bins' weights over
# A(m) = int phi(z) prod_i (omega + (1 - omega) Phi((m_i + sqrt(rho) z) /
# sqrt(1 - rho))) dz.
# ---------------------------------------------------------------------------- #

.location_dense_setup <- function() {

  list(y = c(.30, .55, .10, .65, .20, .45, -.10, .35, .60),
       x = c(-1, 0, 1, -.5, .5, 1, -1, -.2, .8),
       block = rep(1:3, each = 3L), rho = .4, omega = .2, prior_mean = c(.2, 0),
       rules = .coarse_rules(c(7L, 15L, 31L, 63L, 127L, 255L, 511L, 1023L)))
}


.location_dense_model <- function(setup) {

  covariance <- matrix(setup$rho, 3L, 3L)
  diag(covariance) <- 1
  design <- BayesTools::selection_qmc_design(6L, 512L, 8L, 9049L)
  child <- function(s) {
    paste0("  y", s, "[1:3] ~ dselnorm_mnorm_step(m", s, "[1:3], covariance[1:6], sei[1:3], omega[1:2], ",
           "z_lower[1:2], z_upper[1:2], bin", s, "[1:3], 1, 1, 1, qmc[1:nqmc], 512, 8, 0.005, 0, nodes[1:nq], ",
           "log_weights[1:nq], orders[1:nr])")
  }
  means <- unlist(lapply(1:3, function(s) {
    rows <- which(setup$block == s)
    paste0("  m", s, "[", seq_along(rows), "] <- b[1] + b[2] * x[", rows, "]")
  }))
  syntax <- c(
    "model {",
    "  for (j in 1:2) {",
    "    g[j] ~ dcat(gate_weights)",
    "    slab[j] ~ dnorm(prior_mean[j], 4)",
    "    b[j] <- 0 * (g[j] == 1) + slab[j] * (g[j] == 2)",
    "  }",
    means,
    vapply(1:3, child, character(1)),
    paste0("  tag_group_", 1:2, " <- robma_location_group(g[", 1:2, "], slab[", 1:2, "], 2, b[", 1:2, "], ", 1:2, ")"),
    "  tag_setup <- robma_location_setup(g[1], 2, 0)",
    "}")
  data <- list(x = setup$x, prior_mean = setup$prior_mean, gate_weights = c(1, 1),
               covariance = covariance[lower.tri(covariance, diag = TRUE)], sei = rep(1, 3L),
               omega = c(1, setup$omega), z_lower = c(0, -1e300), z_upper = c(1e300, 0),
               qmc = as.numeric(design), nqmc = length(design), nodes = setup$rules$nodes,
               log_weights = setup$rules$log_weights, orders = setup$rules$orders,
               nq = length(setup$rules$nodes), nr = length(setup$rules$orders))
  for (s in 1:3) {
    rows <- which(setup$block == s)
    data[[paste0("y", s)]]   <- setup$y[rows]
    data[[paste0("bin", s)]] <- ifelse(setup$y[rows] > 0, 1L, 2L)
  }
  list(syntax = paste(syntax, collapse = "\n"), data = data)
}


.location_dense_reference <- function(setup, nodes = 20L) {

  Sigma <- matrix(setup$rho, 3L, 3L)
  diag(Sigma) <- 1
  R <- chol(Sigma)
  log_block <- function(y, m) {
    residual <- backsolve(R, y - m, transpose = TRUE)
    weight <- ifelse(y > 0, 1, setup$omega)
    A <- stats::integrate(function(z) {
      vapply(z, function(zz) {
        stats::dnorm(zz) * prod(setup$omega + (1 - setup$omega) *
                                  stats::pnorm((m + sqrt(setup$rho) * zz) / sqrt(1 - setup$rho)))
      }, numeric(1))
    }, -Inf, Inf, rel.tol = 1e-10)$value
    -sum(log(diag(R))) - 1.5 * log(2 * pi) - 0.5 * sum(residual^2) + sum(log(weight)) - log(A)
  }
  loglik <- function(B) {
    vapply(seq_len(nrow(B)), function(k) {
      m <- B[k, 1] + B[k, 2] * setup$x
      sum(vapply(1:3, function(s) {
        rows <- which(setup$block == s)
        log_block(setup$y[rows], m[rows])
      }, numeric(1)))
    }, numeric(1))
  }
  patterns <- as.matrix(expand.grid(g1 = 0:1, g2 = 0:1))
  results <- lapply(seq_len(nrow(patterns)), function(k) {
    on <- which(patterns[k, ] == 1)
    lp <- 2 * log(0.5)
    if (length(on) == 0L) return(list(log_mass = lp + loglik(matrix(0, 1, 2)), m1 = rep(0, 2), m2 = rep(0, 2)))
    full <- function(beta) {
      beta <- matrix(beta, ncol = length(on))
      B <- matrix(0, nrow(beta), 2)
      B[, on] <- beta
      prior <- rep(0, nrow(beta))
      for (i in seq_along(on)) prior <- prior + stats::dnorm(beta[, i], setup$prior_mean[on[i]], 0.5, log = TRUE)
      loglik(B) + prior
    }
    fit <- stats::optim(setup$prior_mean[on], function(beta) -full(beta), method = "BFGS", hessian = TRUE)
    L <- t(chol(solve(fit$hessian)))
    rule <- .coarse_reference_rule(nodes)
    grid <- as.matrix(expand.grid(rep(list(seq_len(nodes)), length(on))))
    x <- matrix(rule$nodes[grid], ncol = length(on))
    w <- apply(matrix(rule$weights[grid], ncol = length(on)), 1L, prod)
    beta <- sweep(x %*% t(L), 2L, fit$par, `+`)
    log_f <- full(beta) + 0.5 * rowSums(x^2)
    top <- max(log_f)
    weight <- w * exp(log_f - top)
    mass <- sum(weight)
    m1 <- m2 <- rep(0, 2)
    m1[on] <- colSums(weight * beta) / mass
    m2[on] <- colSums(weight * beta^2) / mass
    list(log_mass = lp + top + log(mass) + sum(log(diag(L))) + 0.5 * length(on) * log(2 * pi), m1 = m1, m2 = m2)
  })
  log_mass <- vapply(results, `[[`, numeric(1), "log_mass")
  p <- exp(log_mass - max(log_mass))
  p <- p / sum(p)
  M1 <- colSums(p * t(vapply(results, `[[`, numeric(2), "m1")))
  M2 <- colSums(p * t(vapply(results, `[[`, numeric(2), "m2")))
  c(g1 = sum(p * patterns[, 1]), g2 = sum(p * patterns[, 2]), b1 = M1[[1]], b2 = M1[[2]],
    b1_2 = M2[[1]], b2_2 = M2[[2]])
}


.location_dense_observables <- function(draws) {

  cbind(g1 = draws[, "g[1]"] == 2, g2 = draws[, "g[2]"] == 2, b1 = draws[, "b[1]"], b2 = draws[, "b[2]"],
        b1_2 = draws[, "b[1]"]^2, b2_2 = draws[, "b[2]"]^2) * 1
}
