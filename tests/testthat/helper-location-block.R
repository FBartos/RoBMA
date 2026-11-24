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
    return(data.frame(type = integer(), selection = logical(), stage1 = numeric(), full = numeric(),
                      selection_part = numeric(), log_odds = numeric()))
  }
  m <- matrix(records, ncol = 6L, byrow = TRUE)
  data.frame(type = c("flip_on", "flip_off", "refresh", "study", "swap", "joint", "ellipse",
                      "family", "slope", "label_swap")[m[, 1L] + 1L],
             selection = m[, 2L] == 1, stage1 = m[, 3L], full = m[, 4L], selection_part = m[, 5L],
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
# N(0, 1) priors; the gate probabilities and the first two moments of the
# physical coefficients.
.location_selection_reference <- function(setup, model, coefficient_nodes = 16L, latent_nodes = 32L) {

  latent <- .coarse_reference_rule(latent_nodes)
  X <- cbind(1, setup$x1, setup$x2)
  # log of the likelihood integrated over the cluster effects, at coefficient
  # points (rows of B, three columns)
  log_marginal <- function(B) {
    total <- numeric(nrow(B))
    for (s in 1:3) {
      rows  <- which(setup$cluster == s)
      base  <- B %*% t(X[rows, , drop = FALSE])
      means <- base[rep(seq_len(nrow(B)), each = latent_nodes), , drop = FALSE] +
        setup$scale * latent$nodes[rep(seq_len(latent_nodes), nrow(B))]
      ll  <- matrix(.location_selection_cluster_loglik(setup, model, s, means), nrow = latent_nodes)
      top <- apply(ll, 2L, max)
      total <- total + top + log(colSums(latent$weights * exp(sweep(ll, 2L, top))))
    }
    total
  }
  patterns <- as.matrix(expand.grid(g1 = 0:1, g2 = 0:1, g3 = 0:1))
  log_gate <- log(setup$weights / sum(setup$weights))
  results <- lapply(seq_len(nrow(patterns)), function(k) {
    on <- which(patterns[k, ] == 1)
    log_prior_gate <- sum(log_gate[patterns[k, ] + 1])
    if (length(on) == 0L) {
      return(list(log_mass = log_prior_gate + log_marginal(matrix(0, 1, 3)), m1 = rep(0, 3), m2 = rep(0, 3)))
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
      log_marginal(B) + prior
    }
    fit <- stats::optim(setup$prior_mean[on], function(beta) -full(beta), method = "BFGS", hessian = TRUE)
    L <- t(chol(solve(fit$hessian)))
    rule <- .coarse_reference_rule(coefficient_nodes)
    grid <- as.matrix(expand.grid(rep(list(seq_len(coefficient_nodes)), length(on))))
    x <- matrix(rule$nodes[grid], ncol = length(on))
    w <- apply(matrix(rule$weights[grid], ncol = length(on)), 1L, prod)
    beta <- sweep(x %*% t(L), 2L, fit$par, `+`)
    # the N(0, I) weight of the rule is divided out
    log_f <- full(beta) + 0.5 * rowSums(x^2)
    top <- max(log_f)
    weight <- w * exp(log_f - top)
    mass <- sum(weight)
    m1 <- m2 <- rep(0, 3)
    m1[on] <- colSums(weight * beta) / mass
    m2[on] <- colSums(weight * beta^2) / mass
    list(log_mass = log_prior_gate + top + log(mass) + sum(log(diag(L))) + 0.5 * length(on) * log(2 * pi),
         m1 = m1, m2 = m2)
  })
  log_mass <- vapply(results, `[[`, numeric(1), "log_mass")
  p  <- exp(log_mass - max(log_mass))
  p  <- p / sum(p)
  m1 <- colSums(p * t(vapply(results, `[[`, numeric(3), "m1")))
  m2 <- colSums(p * t(vapply(results, `[[`, numeric(3), "m2")))
  c(g1 = sum(p * patterns[, 1]), g2 = sum(p * patterns[, 2]), g3 = sum(p * patterns[, 3]),
    b1 = m1[[1]], b2 = m1[[2]], b3 = m1[[3]], b1_2 = m2[[1]], b2_2 = m2[[2]], b3_2 = m2[[3]])
}


.location_selection_observables <- function(draws) {

  b <- draws[, c("b[1]", "b[2]", "b[3]")]
  cbind(g1 = draws[, "g[1]"] == 2, g2 = draws[, "g[2]"] == 2, g3 = draws[, "g[3]"] == 2,
        b1 = b[, 1], b2 = b[, 2], b3 = b[, 3], b1_2 = b[, 1]^2, b2_2 = b[, 2]^2, b3_2 = b[, 3]^2) * 1
}
