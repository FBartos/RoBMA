# Fixtures of the row cache of singleton step-selection rows
# (src/samplers/StepRowCache.h): hand-written rjags models whose selection
# children are dselnorm_step_switch rows, the cache's test hook, and an exact
# likelihood of the rows in R for the quadrature references.
#
# Label fixture (.rows_setup, .rows_model): rows y_i = direction z_i se_i with
# means mu_i + pet se_i equals(b, 4), SDs sqrt(se_i^2 + tau^2), powers pw_i,
# three bins of direction y / se (z >= 1.96, 0 <= z < 1.96, z < 0; open at both
# ends) and a label b over
#   1  no selection (normal kernel);
#   2  a two-bin weight function: r ~ Beta(2, 2), weights (1, r, r);
#   3  a three-bin weight function: eta_j ~ Gamma(1, 1), p = eta / sum(eta),
#      weights (1, p2 + p3, p3);
#   4  PET (normal kernel, the means move by pet se_i);
#   5  fixed weights (1, w5, w5).
# r, eta and pet carry gate markers active at b = 2, 3 and 4. Row 2 is
# observed at the cut point 1.96, rows 2, 4 and 6 have powers other than one,
# and row 8 lies far below the cut points: under label 5 (w5 = 1e-300) its
# normalizer is below the linear threshold, so the cache evaluates it by its
# node (the log-scale path).

.rows_setup <- function(direction = 1, telescope = 1, tau = .1,
                        z = c(2.4, 1.96, .6, -.3, 1.1, 2.9, 3.5, -.9),
                        se = c(.2, .15, .3, .25, .1, .2, .12, .05),
                        mu = c(.3, .3, .3, .3, .3, .3, .3, -5),
                        power = c(1, .5, 1, 2, 1, .7, 1, 1),
                        w5 = c(1, 1e-300, 1e-300)) {

  obs_bin <- ifelse(z >= 1.96, 1, ifelse(z >= 0, 2, 3))
  list(direction = direction, telescope = telescope, z = z, se = se, y = direction * z * se,
       mu = direction * mu, sd = sqrt(se^2 + tau^2), power = power, obs_bin = obs_bin, w5 = w5,
       z_lower = c(1.96, 0, -1e300), z_upper = c(1e300, 1.96, 0))
}


# The label fixture's model with the corrected sampler options
# (.block_hand_model), initialized at label `b` and the given weights;
# `heterogeneity` adds an indicator h ~ Bernoulli(1/2) of the rows' SDs
# sqrt(sd^2 + h 0.04).
.rows_model <- function(setup, b = 1L, r = .4, eta = c(1, 2, .5), pet = .8, prior = rep(.2, 5),
                        heterogeneity = FALSE, ...) {

  sd <- if (heterogeneity) "sqrt(sd[i] * sd[i] + h * 0.04)" else "sd[i]"
  syntax <- paste0("model {\n",
    if (heterogeneity) "  h ~ dbern(0.5)\n",
    "  b ~ dcat(prior[1:5])\n",
    "  r ~ dbeta(2, 2)\n",
    "  gate_r ~ dbt_gate(r, equals(b, 2))\n",
    "  for (j in 1:3) {\n",
    "    eta[j] ~ dgamma(1, 1)\n",
    "    p[j] <- eta[j] / sum(eta)\n",
    "    gate_eta[j] ~ dbt_gate(eta[j], equals(b, 3))\n",
    "  }\n",
    "  pet ~ dnorm(0, 1)T(0,)\n",
    "  gate_pet ~ dbt_gate(pet, equals(b, 4))\n",
    "  omega[1] <- 1\n",
    "  omega[2] <- equals(b, 1) + r * equals(b, 2) + (p[2] + p[3]) * equals(b, 3) + equals(b, 4) + ",
    "w5[2] * equals(b, 5)\n",
    "  omega[3] <- equals(b, 1) + r * equals(b, 2) + p[3] * equals(b, 3) + equals(b, 4) + w5[3] * equals(b, 5)\n",
    "  mode <- 1 - equals(b, 1) - equals(b, 4)\n",
    "  for (i in 1:N) {\n",
    "    m[i] <- mu[i] + pet * se[i] * equals(b, 4)\n",
    "    y[i] ~ dselnorm_step_switch(m[i], ", sd, ", se[i], pw[i], omega[1:3], z_lower[1:3], z_upper[1:3], ",
    "obs[i], direction, mode, telescope)\n",
    "  }\n}")
  data <- list(y = setup$y, mu = setup$mu, sd = setup$sd, se = setup$se, pw = setup$power, obs = setup$obs_bin,
               N = length(setup$y), z_lower = setup$z_lower, z_upper = setup$z_upper, w5 = setup$w5,
               direction = setup$direction, telescope = setup$telescope, prior = prior,
               gate_r = 0, gate_eta = rep(0, 3), gate_pet = 0)
  inits <- list(b = b, r = r, eta = eta, pet = pet)
  if (heterogeneity) inits$h <- 1L
  .block_hand_model(syntax, data, inits = inits, ...)
}


# A small RoBMA() meta-regression on the Assink data (singleton step-selection
# rows, nine bias labels; `nested` the nested weight layout, `joint` the
# sampler.joint policy): its JAGS fit object (model syntax, data, end state,
# sampler table).
.rows_robma_fit <- function(nested, joint = "auto") {

  data("dat.assink2016", package = "metadat", envir = environment())
  dat <- dat.assink2016[dat.assink2016$study %in% c(1, 2, 4, 12, 13), ]
  previous <- RoBMA::RoBMA.options()[c("sampler.nested.weights", "selection.sampler", "sampler.joint")]
  on.exit(do.call(RoBMA::RoBMA.options, previous), add = TRUE)
  RoBMA::RoBMA.options(sampler.nested.weights = nested, selection.sampler = "slice.coarse", sampler.joint = joint)
  fit <- suppressWarnings(RoBMA::RoBMA(
    yi = dat$yi, sei = sqrt(dat$vi), mods = ~ deltype, measure = "SMD", data = dat, chains = 1L, adapt = 100L,
    burnin = 100L, sample = 100L, seed = 1L, parallel = FALSE, silent = TRUE,
    convergence_checks = RoBMA::set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)))
  fit[["fit"]]
}

# The fit's model compiled again (.block_hand_model) at its end state with the
# bias indicator at `label`; deterministic entries of the end state (any NA)
# are left out of the initial values.
.rows_recompile <- function(fit, label) {

  invisible(gc())
  state <- fit[["end_state"]][[1L]]
  state <- state[!names(state) %in% c(".RNG.state", ".RNG.name")]
  state <- state[!vapply(state, anyNA, logical(1L))]
  state[["bias_indicator"]] <- label
  .block_hand_model(paste(as.character(fit[["model"]]), collapse = "\n"), fit[["data"]], inits = state)
}


# The test hook of the row cache: flags of the samplers built from now on,
# counters (cleared on request) and the frozen-state probe of the live blocks.
.rows_test_hook <- function(flags = NULL, clear = FALSE, probe = FALSE) {

  .Call("RoBMA_step_row_cache_test", flags, clear, probe, PACKAGE = "RoBMA")
}

# Runs `code` with the row cache's test flags set, restoring the defaults.
.with_row_flags <- function(flags, code) {

  .rows_test_hook(flags)
  on.exit(.rows_test_hook(list(mutant = 0, perturbation = 0, disable = FALSE)), add = TRUE)
  force(code)
}

# The probe records as a matrix (chain, cache value, sum of the rows' graph log
# densities), one row per evaluated point.
.rows_probe <- function() {

  values <- .rows_test_hook(probe = TRUE)$probe
  matrix(values, ncol = 3L, byrow = TRUE, dimnames = list(NULL, c("chain", "cache", "graph")))
}

# The largest relative difference |a - b| / max(1, |b|) over finite pairs; the
# pairs must agree on which values are -Inf.
.rows_relative_gap <- function(a, b) {

  expect_identical(is.finite(a), is.finite(b))
  finite <- is.finite(a) & is.finite(b)
  if (!any(finite)) return(0)
  max(abs(a[finite] - b[finite]) / pmax(1, abs(b[finite])))
}

# The frozen-state gaps of a compiled model: the screened indicator's surrogate
# against its exact score of every level, and the live blocks' cache values
# against the rows' graph densities at their probe points (NA without a
# screened indicator or a live block).
.rows_frozen_gaps <- function() {

  scores <- .weights_screened_scores()
  probe <- .rows_probe()
  c(label = if (length(scores) > 0L) max(vapply(scores, function(s) .rows_relative_gap(s$surrogate, s$exact),
                                                numeric(1L))) else NA_real_,
    weights = if (nrow(probe) > 0L) .rows_relative_gap(probe[, "cache"], probe[, "graph"]) else NA_real_)
}


# The exact log likelihood of the rows of a setup at weights `omega` and kernel
# `mode` (0 normal, 1 step), from the definition of the step-selection density:
# pw [log phi(y; m, sd) + log omega_obs - log sum_b omega_b P(lower_b < direction
# Y / se < upper_b)] in step mode, pw log phi(y; m, sd) in normal mode.
.rows_children_log <- function(setup, omega, mode) {

  y  <- setup$y
  m  <- setup$mu
  sd <- setup$sd
  se <- setup$se
  gaussian <- stats::dnorm(y, m, sd, log = TRUE)
  if (mode == 0) return(sum(setup$power * gaussian))
  lower <- ifelse(setup$z_lower <= -1e299, -Inf, setup$z_lower)
  upper <- ifelse(setup$z_upper >= 1e299, Inf, setup$z_upper)
  shift <- setup$direction * m
  normalizer <- vapply(seq_along(y), function(i) {
    sum(omega * (stats::pnorm((upper * se[i] - shift[i]) / sd[i]) - stats::pnorm((lower * se[i] - shift[i]) / sd[i])))
  }, numeric(1L))
  sum(setup$power * (gaussian + log(omega[setup$obs_bin]) - log(normalizer)))
}


# The weights fixture of helper-selection-corrected-weights.R with singleton
# step-selection rows: the rows of `x` (by row) at the fixed mean, SD 1 and
# selection SE 1, labels 1 (normal kernel), 2 (r) and 3 (eta) as there.
.rows_weights_setup <- function(branches = 1:3, x, mean = .25, prior = NULL) {

  z <- as.vector(t(x))
  if (is.null(prior)) prior <- rep(1, length(branches)) / length(branches)
  list(branches = branches, direction = 1, telescope = 1, z = z, y = z, se = rep(1, length(z)),
       mu = rep(mean, length(z)), sd = rep(1, length(z)), power = rep(1, length(z)),
       obs_bin = ifelse(z >= 1.96, 1, ifelse(z >= 0, 2, 3)), z_lower = c(1.96, 0, -1e300),
       z_upper = c(1e300, 1.96, 0), alpha = c(1, 1, 1), beta = c(1, 1), prior = prior)
}

.rows_weights_model <- function(setup, ...) {

  branches <- setup$branches
  label <- function(branch) match(branch, branches)
  composition <- function(bin) {
    terms <- c(
      if (1L %in% branches) paste0("equals(b, ", label(1L), ")"),
      if (2L %in% branches) paste0(c("1", "r", "r")[bin], " * equals(b, ", label(2L), ")"),
      if (3L %in% branches) paste0(c("1", "(p[2] + p[3])", "p[3]")[bin], " * equals(b, ", label(3L), ")"))
    paste0("  omega[", bin, "] <- ", paste(terms, collapse = " + "))
  }
  mode <- paste0("1", if (1L %in% branches) paste0(" - equals(b, ", label(1L), ")"))
  syntax <- paste0("model {\n",
    "  b ~ dcat(prior[1:", length(branches), "])\n",
    if (2L %in% branches) paste0("  r ~ dbeta(beta[1], beta[2])\n",
                                 "  gate_r ~ dbt_gate(r, equals(b, ", label(2L), "))\n"),
    if (3L %in% branches) paste0(
      "  for (j in 1:3) {\n    eta[j] ~ dgamma(alpha[j], 1)\n    p[j] <- eta[j] / sum(eta)\n",
      "    gate_eta[j] ~ dbt_gate(eta[j], equals(b, ", label(3L), "))\n  }\n"),
    paste(vapply(1:3, composition, ""), collapse = "\n"), "\n",
    "  mode <- ", mode, "\n",
    "  for (i in 1:N) {\n",
    "    y[i] ~ dselnorm_step_switch(mu[i], sd[i], se[i], 1, omega[1:3], z_lower[1:3], z_upper[1:3], ",
    "obs[i], 1, mode, 1)\n",
    "  }\n}")
  data <- list(y = setup$y, mu = setup$mu, sd = setup$sd, se = setup$se, obs = setup$obs_bin, N = length(setup$y),
               z_lower = setup$z_lower, z_upper = setup$z_upper, prior = setup$prior)
  inits <- list(b = 1L)
  if (2L %in% branches) {
    data <- c(data, list(beta = setup$beta, gate_r = 0))
    inits$r <- .5
  }
  if (3L %in% branches) {
    data <- c(data, list(alpha = setup$alpha, gate_eta = rep(0, 3)))
    inits$eta <- c(1, 1, 1)
  }
  .block_hand_model(syntax, data, inits = inits, ...)
}

# z-scores of a chain of the step-row weights fixture against the reference.
.rows_weights_z <- function(setup, reference, draws, flags = list(), row_flags = list(), fixed_order = FALSE) {

  monitor <- c("b", if (2L %in% setup$branches) "r", if (3L %in% setup$branches) "eta")
  .corrected_test_hook(clear = TRUE)
  chain <- .with_row_flags(row_flags, .with_corrected_flags(flags, .block_chain(
    .rows_weights_model(setup, fixed_order = fixed_order), monitor, draws, adapt = 500L, burnin = 500L)))
  observables <- .weights_observables(chain$draws, setup)
  targets <- c(reference$estimate, if (3L %in% setup$branches) .weights_total_reference(setup))
  allowance <- c(reference$allowance, if (3L %in% setup$branches) c(T = 0, T2 = 0))
  comparison <- .weights_compare(observables, targets[colnames(observables)], allowance[colnames(observables)])
  comparison <- comparison[is.finite(comparison$z), ]
  # the corrected samplers' proposals and accepted corrections (adaptation included)
  attr(comparison, "kinds") <- .corrected_test_hook(clear = TRUE)$kinds
  comparison
}
