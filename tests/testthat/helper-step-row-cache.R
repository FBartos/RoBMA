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

# The fit's model compiled again (.block_hand_model, one at a time:
# .fresh_model) at its end state with the bias indicator at `label`;
# deterministic entries of the end state (any NA) are left out of the initial
# values.
.rows_recompile <- function(fit, label) {

  state <- fit[["end_state"]][[1L]]
  state <- state[!names(state) %in% c(".RNG.state", ".RNG.name")]
  state <- state[!vapply(state, anyNA, logical(1L))]
  state[["bias_indicator"]] <- label
  .fresh_model(.block_hand_model(paste(as.character(fit[["model"]]), collapse = "\n"), fit[["data"]], inits = state))
}


# The public RoBMA() fixture of a narrow central bin: one estimate, fixed
# effect and heterogeneity, a fixed weight function (1, 1e60, 1) on the
# one-sided cut points 0.5 - eps / 2 and 0.5 (a central bin of width 2.8e-16
# on the z scale) against a fixed PET slope, without joint samplers (the
# screened indicator owns the bias label). The branch prior makes both labels
# about equally probable. The linear interval mass of the central bin rounds
# to zero, while the other bins keep the cached normalizer above the linear
# threshold. Returns the backend fit (model syntax, data).
.rows_narrow_bin_fit <- function() {

  steps <- c(.5 - .Machine$double.eps / 2, .5)
  cut   <- stats::qnorm(steps[1L], lower.tail = FALSE)
  tau   <- sqrt(15)
  y     <- cut / 2
  pet   <- .2
  relative_mass <- stats::dnorm(y, y + pet, sqrt(1 + tau^2)) * cut
  bias <- list(
    BayesTools::prior_weightfunction("one-sided", steps, BayesTools::wf_fixed(c(1, 1e60, 1)),
                                     prior_weights = relative_mass),
    BayesTools::prior_PET("point", list(location = pet), prior_weights = 1))
  # the 1e60 weight disables probability telescoping (a warning of the fixture)
  suppressWarnings(.rows_screened_fit(y, effect = y, heterogeneity = tau, bias = bias))
}

# The public RoBMA() fixture of two tied labels: one estimate at zero, fixed
# effect and heterogeneity zero, a flat one-sided weight function (1, 1)
# against no bias, equal prior weights: both labels have the same likelihood,
# so each has probability 1/2. Returns the backend fit (model syntax, data).
.rows_tie_fit <- function() {

  bias <- list(BayesTools::prior_weightfunction("one-sided", .5, BayesTools::wf_fixed(c(1, 1))),
               BayesTools::prior_none())
  .rows_screened_fit(0, effect = 0, heterogeneity = 0, bias = bias)
}

# A one-estimate RoBMA() fit with fixed effect and heterogeneity and the bias
# priors `bias`, without joint samplers (the screened indicator owns the bias
# label). RoBMA() is asked for the object only and .fit() draws, since the
# summary of the fixed parameters is not available.
.rows_screened_fit <- function(y, effect, heterogeneity, bias) {

  fixed <- function(location) BayesTools::prior("point", list(location = location))
  previous <- RoBMA::RoBMA.options()[c("sampler.joint", "selection.sampler")]
  on.exit(do.call(RoBMA::RoBMA.options, previous), add = TRUE)
  RoBMA::RoBMA.options(sampler.joint = "off", selection.sampler = "slice.coarse")
  object <- RoBMA::RoBMA(
    yi = y, sei = 1, measure = "SMD", effect_direction = "positive", prior_effect = fixed(effect),
    prior_effect_null = FALSE, prior_heterogeneity = fixed(heterogeneity), prior_heterogeneity_null = FALSE,
    prior_bias = bias, prior_bias_null = FALSE, seed = 15L, sample = 100L, burnin = 50L, adapt = 50L,
    chains = 1L, autofit = FALSE, silent = TRUE, only_priors = TRUE,
    convergence_checks = RoBMA::set_convergence_checks(min_ESS = NULL, max_Rhat = NULL,
                                                       allow_not_assessable = TRUE))
  .fit(object)
}

# The screened indicator's frozen scores of a model compiled at entry level
# `entry` of node `node` (the scores of the most recent live sampler).
.rows_entry_scores <- function(syntax, data, node, entry, inits = list()) {

  inits[[node]] <- entry
  hand <- .fresh_model(.block_hand_model(syntax, data, inits = inits))
  scores <- .weights_screened_scores()
  scores[[length(scores)]]
}

# The stationary probability of level 2 of a two-level screened update whose
# moves from level 1 use the surrogate scores `from1` and those from level 2
# `from2` (each the scores of both levels), with exact scores `exact`: the
# proposal's probability of the other level (Liu's step: the other level,
# accepted with min(1, p_other / p_current); the plain draw: p_other) times
# the exact correction. With one surrogate for both ends it is the exact
# conditional probability for either proposal.
.rows_two_state <- function(from1, from2, exact, proposal = c("liu", "plain")) {

  proposal <- match.arg(proposal)
  move <- function(scores, from, to) {
    d <- scores[to] - scores[from]
    first <- if (proposal == "liu") exp(min(d, 0)) else stats::plogis(d)
    first * exp(min((exact[to] - scores[to]) - (exact[from] - scores[from]), 0))
  }
  up   <- move(from1, 1L, 2L)
  down <- move(from2, 2L, 1L)
  up / (up + down)
}


# The nested weight layout with a PET label on singleton step-selection rows
# (the row cache): 16 rows at a fixed mean and SD = SE, four bins (z >= 2.326,
# [1.96, 2.326), [0, 1.96), z < 0) and a label b over
#   1  no selection (normal kernel);
#   2  the coarse weight function (1, 1, a, a v), a ~ Beta(2, 1);
#   3  the fine weight function (1, a + (1 - a) u, a, a v), a ~ Beta(2, 2);
#   4  PET (normal kernel, the means move by pet se_i), pet ~ N+(0, 1);
# with the label-dependent prior of the shared coordinate a (BayesTools'
# nested layout), gates {a, v} open under labels 2 and 3, {u} under 3 and
# {pet} under 4, sampled by the gated block, the screened indicator and the
# gated slice sampler (pet). The prior balances the labels' posterior
# probabilities (each about 1/4).
.rows_nested_setup <- function() {

  se <- rep(c(.35, .6, .9, 1.3), 4)
  z  <- c(2.5, 2.1, 1.0, -.4, 2.4, 1.7, .3, -.9, 2.8, 2.0, .7, .1, 2.2, 1.2, .5, -.3)
  setup <- list(se = se, z = z, y = z * se, mu = .1, sd = se, z_lower = c(2.326, 1.96, 0, -1e300),
                z_upper = c(1e300, 2.326, 1.96, 0),
                obs = ifelse(z >= 2.326, 1, ifelse(z >= 1.96, 2, ifelse(z >= 0, 3, 4))))
  base <- .rows_nested_reference_at(setup, rep(.25, 4), 24L)
  prior <- .25 / base[c("I1", "I2", "I3", "I4")]
  setup$prior <- unname(prior / sum(prior))
  setup
}

.rows_nested_model <- function(setup, inits = list(b = 1L, a = .5, v = .5, u = .5, pet = .2), ...) {

  syntax <- paste0("model {\n",
    "  b ~ dcat(prior[1:4])\n",
    "  a ~ dbeta(1 + equals(b, 2) + equals(b, 3), 1 + equals(b, 3))\n",
    "  v ~ dbeta(1, 1)\n  u ~ dbeta(1, 1)\n  pet ~ dnorm(0, 1)T(0,)\n",
    "  gate_a ~ dbt_gate(a, equals(b, 2) + equals(b, 3))\n",
    "  gate_v ~ dbt_gate(v, equals(b, 2) + equals(b, 3))\n",
    "  gate_u ~ dbt_gate(u, equals(b, 3))\n",
    "  gate_pet ~ dbt_gate(pet, equals(b, 4))\n",
    "  omega[1] <- 1\n",
    "  omega[2] <- equals(b, 1) + equals(b, 4) + equals(b, 2) + (a + (1 - a) * u) * equals(b, 3)\n",
    "  omega[3] <- equals(b, 1) + equals(b, 4) + a * equals(b, 2) + a * equals(b, 3)\n",
    "  omega[4] <- equals(b, 1) + equals(b, 4) + a * v * equals(b, 2) + a * v * equals(b, 3)\n",
    "  mode <- equals(b, 2) + equals(b, 3)\n",
    "  for (i in 1:N) {\n",
    "    m[i] <- mu + pet * se[i] * equals(b, 4)\n",
    "    y[i] ~ dselnorm_step_switch(m[i], sd[i], se[i], 1, omega[1:4], z_lower[1:4], z_upper[1:4], obs[i], 1, ",
    "mode, 1)\n",
    "  }\n}")
  data <- list(y = setup$y, mu = setup$mu, sd = setup$sd, se = setup$se, obs = setup$obs, N = length(setup$y),
               z_lower = setup$z_lower, z_upper = setup$z_upper, prior = setup$prior, gate_a = 0, gate_v = 0,
               gate_u = 0, gate_pet = 0)
  .block_hand_model(syntax, data, inits = inits, ...)
}

# The label probabilities and moments of the nested fixture at branch prior
# `prior` by Gauss-Legendre rules of order n over a, v, u on (0, 1) and pet on
# (0, 8), from the rows' likelihood written in R.
.rows_nested_reference_at <- function(setup, prior, n) {

  y <- setup$y
  lower <- ifelse(setup$z_lower <= -1e299, -Inf, setup$z_lower)
  upper <- ifelse(setup$z_upper >= 1e299, Inf, setup$z_upper)
  bins <- t(vapply(seq_along(y), function(i) {
    stats::pnorm((upper * setup$se[i] - setup$mu) / setup$sd[i]) -
      stats::pnorm((lower * setup$se[i] - setup$mu) / setup$sd[i])
  }, numeric(4L)))
  gaussian <- sum(stats::dnorm(y, setup$mu, setup$sd, log = TRUE))
  # the log likelihood of the step labels at weight matrices W (points x 4),
  # relative to the normal kernel
  step_log <- function(W) {
    as.numeric(log(W[, setup$obs, drop = FALSE]) %*% rep(1, length(y))) - rowSums(log(W %*% t(bins)))
  }
  rule <- .weights_legendre(n)
  Z <- numeric(4L)
  m <- list()
  Z[1L] <- prior[1L]
  grid <- expand.grid(i = seq_len(n), j = seq_len(n))
  a <- rule$nodes[grid$i]
  v <- rule$nodes[grid$j]
  w <- rule$weights[grid$i] * rule$weights[grid$j] * stats::dbeta(a, 2, 1) * exp(step_log(cbind(1, 1, a, a * v)))
  Z[2L] <- prior[2L] * sum(w)
  m$I2_a  <- prior[2L] * sum(w * a)
  m$I2_v  <- prior[2L] * sum(w * v)
  m$I2_a2 <- prior[2L] * sum(w * a^2)
  grid <- expand.grid(i = seq_len(n), j = seq_len(n), k = seq_len(n))
  a <- rule$nodes[grid$i]
  v <- rule$nodes[grid$j]
  u <- rule$nodes[grid$k]
  w <- rule$weights[grid$i] * rule$weights[grid$j] * rule$weights[grid$k] * stats::dbeta(a, 2, 2) *
    exp(step_log(cbind(1, a + (1 - a) * u, a, a * v)))
  Z[3L] <- prior[3L] * sum(w)
  m$I3_a  <- prior[3L] * sum(w * a)
  m$I3_v  <- prior[3L] * sum(w * v)
  m$I3_u  <- prior[3L] * sum(w * u)
  m$I3_a2 <- prior[3L] * sum(w * a^2)
  pet <- 8 * rule$nodes
  like <- vapply(pet, function(slope) {
    exp(sum(stats::dnorm(y, setup$mu + slope * setup$se, setup$sd, log = TRUE)) - gaussian)
  }, numeric(1L))
  w <- 8 * rule$weights * 2 * stats::dnorm(pet) * like
  Z[4L] <- prior[4L] * sum(w)
  m$I4_pet  <- prior[4L] * sum(w * pet)
  m$I4_pet2 <- prior[4L] * sum(w * pet^2)
  c(I1 = Z[1L], I2 = Z[2L], I3 = Z[3L], I4 = Z[4L], unlist(m)) / sum(Z)
}

# The reference of the fixture's prior (order 32) with the change from order
# 16 as its allowance.
.rows_nested_reference <- function(setup) {

  fine   <- .rows_nested_reference_at(setup, setup$prior, 32L)
  coarse <- .rows_nested_reference_at(setup, setup$prior, 16L)
  list(estimate = fine, allowance = abs(fine - coarse))
}

.rows_nested_observables <- function(draws) {

  b <- draws[, "b"]
  I <- lapply(1:4, function(k) as.numeric(b == k))
  cbind(I1 = I[[1L]], I2 = I[[2L]], I3 = I[[3L]], I4 = I[[4L]],
        I2_a = I[[2L]] * draws[, "a"], I2_v = I[[2L]] * draws[, "v"], I2_a2 = I[[2L]] * draws[, "a"]^2,
        I3_a = I[[3L]] * draws[, "a"], I3_v = I[[3L]] * draws[, "v"], I3_u = I[[3L]] * draws[, "u"],
        I3_a2 = I[[3L]] * draws[, "a"]^2, I4_pet = I[[4L]] * draws[, "pet"], I4_pet2 = I[[4L]] * draws[, "pet"]^2)
}

# z-scores of a chain of the nested fixture against the reference.
.rows_nested_z <- function(setup, reference, draws, flags = list(), row_flags = list()) {

  .corrected_test_hook(clear = TRUE)
  chain <- .with_row_flags(row_flags, .with_corrected_flags(flags, .block_chain(
    .rows_nested_model(setup), c("b", "a", "v", "u", "pet"), draws, adapt = 500L, burnin = 500L)))
  observables <- .rows_nested_observables(chain$draws)
  comparison <- .weights_compare(observables, reference$estimate, reference$allowance)
  attr(comparison, "kinds") <- .corrected_test_hook(clear = TRUE)$kinds
  attr(comparison, "samplers") <- chain$samplers
  comparison
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
  # a single label's indicator is a structural constant
  comparison <- .weights_compare(observables, targets, allowance, constants = .weights_constants(setup))
  # the corrected samplers' proposals and accepted corrections (adaptation included)
  attr(comparison, "kinds") <- .corrected_test_hook(clear = TRUE)$kinds
  comparison
}
