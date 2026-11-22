#!/usr/bin/env Rscript
# Independent reference for the legacy clustered Gaussian model. Requires
# installed development BayesTools and RoBMA using the same JAGS 5 runtime.
# From the RoBMA checkout, Rscript tools/verify-clustered-sampler.R [draws] [seeds]
args <- commandArgs(trailingOnly = TRUE)
n_draws <- if (length(args) >= 1L) as.integer(args[1L]) else 1500L
n_seeds <- if (length(args) >= 2L) as.integer(args[2L]) else 20L
stopifnot(n_draws >= 100L, n_seeds >= 8L)

if (.Platform$OS.type == "windows" && nzchar(Sys.getenv("JAGS_HOME"))) {
  # The agent SDK stores the vendor libraries here. An ordinary installation
  # typically finds them via the system DLL search path.
  for (dll in c("libjags-5.dll", "libjrmath-0.dll")) {
    path <- file.path(Sys.getenv("JAGS_HOME"), "x64", "bin", dll)
    if (file.exists(path)) dyn.load(path)
  }
}

data <- data.frame(
  y = c(.15, .20, -.15, .31, .17, -.21, .08, .26, -.19),
  v = seq(.04, .18, length.out = 9),
  study = rep(letters[1:3], each = 3),
  group = factor(rep(letters[1:3], 3), levels = c("c", "a", "b"))
)
fixture <- function(...) {
  RoBMA::BMA(yi = y, vi = v, cluster = study, mods = ~ group, data = data,
             measure = "SMD", ...)
}

legendre <- function(n) {
  k <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- J[cbind(k + 1L, k)] <- k / sqrt(4 * k^2 - 1)
  e <- eigen(J, symmetric = TRUE)
  list(x = (1 + e$values) / 2, w = e$vectors[1L, ]^2)
}

normal_reference <- function(y, v, group, study, fixed_sd, tau_sd, odds, n) {
  rule <- legendre(n)
  Z <- outer(study, study, "==")
  cells <- expand.grid(intercept = 1:2, group = 1:2, heterogeneity = 1:2)
  p <- length(fixed_sd)
  result <- t(vapply(seq_len(nrow(cells)), function(k) {
    cell <- cells[k, ]
    included <- c(cell$intercept == 2L, rep(cell$group == 2L, p - 1L))
    X <- cbind(1, group)[, included, drop = FALSE]
    baseline <- if (any(included)) tcrossprod(sweep(X, 2L, fixed_sd[included], "*")) else 0
    at_scales <- function(tau, rho) {
      covariance <- diag(v + tau^2 * (1 - rho)) + tau^2 * rho * Z + baseline
      U <- chol(covariance)
      z <- forwardsolve(t(U), y)
      beta <- numeric(p)
      if (any(included)) beta[included] <- fixed_sd[included]^2 * crossprod(X, backsolve(U, z))
      mass <- exp(-0.5 * sum(z^2) - sum(log(diag(U))))
      mass * c(1, beta, tau * sqrt(rho), tau * sqrt(1 - rho), rho)
    }
    integrated <- if (cell$heterogeneity == 1L) {
      at_scales(0, .5)
    } else {
      taus <- tau_sd * qnorm((1 + rule$x) / 2)
      grid <- vapply(taus, function(tau) {
        vapply(rule$x, function(rho) at_scales(tau, rho), numeric(p + 4L)) %*% rule$w
      }, numeric(p + 4L))
      grid %*% rule$w
    }
    as.numeric(integrated) * prod(mapply(function(probabilities, label) probabilities[label], odds, cell))
  }, numeric(p + 4L)))
  mass <- result[, 1L] / sum(result[, 1L])
  inclusion <- vapply(cells, function(cell) sum(mass[cell == 2L]), numeric(1))
  c(stats::setNames(mass, paste0("cell", seq_len(nrow(cells)))), inclusion,
    stats::setNames(colSums(result[, -1L]) / sum(result[, 1L]),
                    c("mu", "contrast1", "contrast2", "between", "within", "rho")))
}

# Compare the complete original log density, log Jacobian and independent
# conditional proposal density with the cancelled formula in both directions.
# Evaluate different numerical scales and positive allocations near both ends.
check_scale_ratio <- function() {
  gamma <- c(-1.3, .7, .15)
  J <- length(gamma)
  half_sd <- .5
  for (tau in c(1e-5, .08, 1, 50)) {
    for (rho in c(1e-5, .15, .8, .99999)) {
      for (multiplier in c(.57, 1.6)) {
        w <- tau * sqrt(1 - rho)
        s <- tau * sqrt(rho)
        s_new <- s * multiplier
        tau_new <- sqrt(w^2 + s_new^2)
        rho_new <- (s_new / tau_new)^2
        gamma_new <- gamma / multiplier
        Q <- sum((s * gamma)^2)
        jacobian <- function(t, b) log(2 * w) + (2 - J) * log(b) - 3 * log(t)
        proposal <- function(b) -J * log(b) - Q / (2 * b^2)
        prior <- function(t, r) dnorm(t, 0, half_sd, log = TRUE) + dbeta(r, 1.5, 2, log = TRUE)
        before <- prior(tau, rho) + sum(dnorm(gamma, log = TRUE)) +
          jacobian(tau, s) - proposal(s)
        after <- prior(tau_new, rho_new) + sum(dnorm(gamma_new, log = TRUE)) +
          jacobian(tau_new, s_new) - proposal(s_new)
        ratio <- prior(tau_new, rho_new) - prior(tau, rho) +
          2 * log(s_new / s) - 3 * log(tau_new / tau)
        reverse <- prior(tau, rho) - prior(tau_new, rho_new) +
          2 * log(s / s_new) - 3 * log(tau / tau_new)
        stopifnot(abs(after - before - ratio) < 2e-10 * (1 + abs(ratio)),
                  abs(reverse + ratio) < 2e-10 * (1 + abs(ratio)),
                  max(abs(s * gamma - s_new * gamma_new)) < 1e-12,
                  abs(w - tau_new * sqrt(1 - rho_new)) < 1e-10)
      }
    }
  }
}

options_before <- RoBMA::RoBMA.options()[c("sampler.joint", "selection.sampler")]
run <- function() {
  on.exit(do.call(RoBMA::RoBMA.options, options_before), add = TRUE)
  RoBMA::RoBMA.options(selection.sampler = "slice")
  check_scale_ratio()
  cat("Forward and reverse scale ratios passed.\n")
  obj <- fixture(only_priors = TRUE)
  priors <- c(obj[["priors"]][["mods"]], obj[["priors"]][["outcome"]])
  # Build the independent covariance reference from the fitted contrast design
  # and the prior declarations used by the public model.
  effects <- obj[["priors"]][["outcome"]][c("mu", "tau")]
  mod_prior <- obj[["priors"]][["mods"]]
  mod_prior <- mod_prior[["group"]]
  if (is.null(effects[["mu"]])) effects[["mu"]] <- priors[["intercept"]]
  get_sd <- function(p) p[[2L]][["parameters"]][["sd"]]
  probs <- function(p) {
    w <- unname(attr(p, "prior_weights", exact = TRUE))
    w / sum(w)
  }
  model <- fixture(adapt = 100L, burnin = 100L, sample = 100L, chains = 2L,
                   seed = 31L, autofit = FALSE, silent = TRUE)
  contrast <- model[["fit"]][["data"]][["mu_data_group"]]
  fixed_sd <- c(get_sd(effects[["mu"]]), rep(get_sd(mod_prior), ncol(contrast)))
  odds <- lapply(list(effects[["mu"]], mod_prior, effects[["tau"]]), probs)
  exact <- function(n) normal_reference(data$y, data$v, contrast, data$study,
                                       fixed_sd, get_sd(effects[["tau"]]), odds, n)
  reference <- exact(96L)
  integration_error <- max(abs(reference - exact(48L)))
  stopifnot(integration_error < 1e-5)
  cat("Max quadrature difference", integration_error, "\n")

  for (policy in c("off", "clustered", "clustered.scale")) {
    RoBMA::RoBMA.options(sampler.joint = policy)
    estimates <- t(vapply(seq.int(4001L, length.out = n_seeds), function(seed) {
      fit <- fixture(adapt = 100L, burnin = 100L, sample = n_draws, chains = 2L,
                     seed = seed, autofit = FALSE, silent = TRUE)
      sampler <- fit[["fit"]][["samplers"]][["Sampler"]]
      expect <- switch(policy, off = 0L, clustered = 6L, clustered.scale = 8L)
      stopifnot(sum(sampler == "RoBMA::JointCluster") == expect)
      draws <- as.matrix(fit[["fit"]][["mcmc"]])
      labels <- draws[, c("mu_intercept_indicator", "mu_group_indicator", "tau_indicator")]
      cell <- labels[, 1L] + 2 * (labels[, 2L] - 1) + 4 * (labels[, 3L] - 1)
      c(tabulate(cell, 8L) / length(cell), colMeans(labels) - 1,
        colMeans(draws[, c("mu_intercept", "mu_group[1]", "mu_group[2]")]),
        mean(draws[, "tau"] * sqrt(draws[, "rho"])),
        mean(draws[, "tau"] * sqrt(1 - draws[, "rho"])), mean(draws[, "rho"]))
    }, numeric(length(reference))))
    mcse <- apply(estimates, 2L, stats::sd) / sqrt(n_seeds)
    z <- (colMeans(estimates) - reference) / mcse
    print(data.frame(policy = policy, target = names(reference), reference,
                     observed = colMeans(estimates), mcse, z), row.names = FALSE)
    stopifnot(all(abs(z) < 5))
  }
}

run()
