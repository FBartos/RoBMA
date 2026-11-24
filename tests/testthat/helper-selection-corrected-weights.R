# Fixtures and quadrature references for the corrected samplers of the
# publication-bias side of model-averaged fits: the gated weight block
# (RoBMA::GatedCorrectedBlock), the screened indicator (RoBMA::ScreenedFinite)
# and the gamma-sum coordinates of the corrected block samplers, in
# hand-written rjags models with rank-one cluster selection children.
#
# Weights fixture: K cluster children of 4 rows each at a fixed mean, observed
# in three p-value bins (z >= 1.96, 0 <= z < 1.96, z < 0), and a bias indicator
# b over the branches
#   1  no selection: weights (1, 1, 1), normal kernel;
#   2  a two-bin cumulative weight function: r ~ Beta(ar, br), weights (1, r, r);
#   3  a three-bin cumulative weight function: eta_j ~ Gamma(alpha_j, 1),
#      p = eta / sum(eta), weights (1, p2 + p3, p3);
# r and eta carry gate markers active at b = 2 and b = 3. The two-branch
# ensemble has branches 1 and 3, the three-branch ensemble all three. The
# posterior label probabilities and weight moments are integrals of the
# package's exact likelihood of the children (the probe of the distribution
# class the sampler evaluates, with the model's schedule and controls) over
# r and over the simplex of p.
#
# Heterogeneity fixture: an inclusion indicator h ~ Bernoulli(1/2), an
# allocation SD s ~ half-normal and Dirichlet allocation weights from two
# gamma members e (a declared block of s, e1, e2), and cluster children whose
# residual SD and loading grow with h s^2 w1 and h s^2 w2 under fixed selection
# weights. The reference integrates over s and w1.

# The crude schedule of .coarse_exactness_setup("cluster"): rule 3 has three
# nodes, so the surrogate (max_rules = 3) is far from the exact normalizer.
.weights_orders <- c(1L, 2L, 3L, 7L, 15L, 31L, 63L, 127L, 255L, 511L, 1023L)

.weights_setup <- function(branches = 1:3, children = 3L, mean = .25, alpha = c(1, 1, 1),
                           beta = c(1, 1), prior = NULL, x = NULL, correlation = .8) {

  # child c (a row of x): by default rows observed in bins 1, 2, 2, 3 at
  # standardized values that differ by child
  if (is.null(x)) {
    x <- rbind(c(2.4, 1.3, .6, -.3), c(2.1, 1.7, .9, -.8), c(2.9, 1.1, .2, -.5),
               c(2.2, 1.5, .4, -.2))[seq_len(children), , drop = FALSE]
  }
  children <- nrow(x)
  rows <- ncol(x)
  obs_bin <- ifelse(x >= 1.96, 1, ifelse(x >= 0, 2, 3))
  if (is.null(prior)) prior <- rep(1, length(branches)) / length(branches)
  list(branches = branches, children = children, rows = rows, mean = mean, x = x, obs_bin = obs_bin,
       loading = rep(sqrt(correlation), rows), residual_sd = rep(sqrt(1 - correlation), rows),
       z_lower = c(1.96, 0, -1e300), z_upper = c(1e300, 1.96, 0),
       alpha = alpha, beta = beta, prior = prior, rules = .coarse_rules(.weights_orders))
}


# The rjags model of the weights fixture with the corrected sampler options
# (.block_hand_model); the label is `b` with values 1..length(branches) mapped
# to the setup's branches. Variants for the ownership tests: `rates` of the eta
# gammas (syntax), the kernel `mode` (syntax; by default the normal kernel for
# branch 1 and the selection kernel otherwise), `means` (syntax of row j's mean;
# by default the fixed mean as data) and `extra` syntax with `extra_data`.
.weights_model <- function(setup, inits = NULL, rates = c("1", "1", "1"), mode = NULL, means = NULL,
                           extra = "", extra_data = list(), ...) {

  branches <- setup$branches
  label <- function(branch) match(branch, branches)
  rows <- setup$rows
  if (is.null(mode)) mode <- paste0("1", if (1L %in% branches) paste0(" - equals(b, ", label(1L), ")"))
  composition <- function(bin) {
    terms <- c(
      if (1L %in% branches) paste0("equals(b, ", label(1L), ")"),
      if (2L %in% branches) paste0(c("1", "r", "r")[bin], " * equals(b, ", label(2L), ")"),
      if (3L %in% branches) paste0(c("1", "(p[2] + p[3])", "p[3]")[bin], " * equals(b, ", label(3L), ")"))
    paste0("  omega[", bin, "] <- ", paste(terms, collapse = " + "))
  }
  syntax <- paste0("model {\n",
    "  b ~ dcat(prior[1:", length(branches), "])\n",
    if (2L %in% branches) paste0("  r ~ dbeta(beta[1], beta[2])\n",
                                 "  gate_r ~ dbt_gate(r, equals(b, ", label(2L), "))\n"),
    if (3L %in% branches) paste0(
      paste0("  eta[", 1:3, "] ~ dgamma(alpha[", 1:3, "], ", rates, ")\n", collapse = ""),
      "  for (j in 1:3) {\n    p[j] <- eta[j] / sum(eta)\n",
      "    gate_eta[j] ~ dbt_gate(eta[j], equals(b, ", label(3L), "))\n  }\n"),
    paste(vapply(1:3, composition, ""), collapse = "\n"), "\n",
    "  mode <- ", mode, "\n",
    if (!is.null(means)) paste0("  for (j in 1:", rows, ") { means[j] <- ", means, " }\n"),
    extra, "\n",
    "  for (c in 1:K) {\n",
    "    y[c, 1:", rows, "] ~ dselnorm_cluster_step(means[1:", rows, "], residual_sd[1:", rows, "], loading[1:", rows,
    "], sei[1:", rows, "], omega[1:3], z_lower[1:3], z_upper[1:3], obs_bin[c, 1:", rows, "], 1, 1, mode, ",
    "nodes[1:nq], log_weights[1:nq], orders[1:nr], qmc[1:nqmc], 64, 512, 8, 0.005, 0)\n",
    "  }\n}")
  design <- BayesTools::selection_qmc_design(2L, 512L, 8L, 9049L)
  data <- list(y = setup$x, means = rep(setup$mean, rows), residual_sd = setup$residual_sd,
               loading = setup$loading, sei = rep(1, rows), z_lower = setup$z_lower, z_upper = setup$z_upper,
               obs_bin = setup$obs_bin, K = setup$children, prior = setup$prior,
               nodes = setup$rules$nodes, log_weights = setup$rules$log_weights, orders = setup$rules$orders,
               nq = length(setup$rules$nodes), nr = length(setup$rules$orders),
               qmc = as.numeric(design), nqmc = length(design))
  if (2L %in% branches) data <- c(data, list(beta = setup$beta, gate_r = 0))
  if (3L %in% branches) data <- c(data, list(alpha = setup$alpha, gate_eta = rep(0, 3)))
  if (!is.null(means)) data[["means"]] <- NULL
  data <- c(data, extra_data)
  if (is.null(inits)) {
    inits <- list(b = 1L)
    if (2L %in% branches) inits$r <- .5
    if (3L %in% branches) inits$eta <- c(1, 1, 1)
  }
  .block_hand_model(syntax, data, inits = inits, ...)
}


# log of the product over the children of the exact density (and of the
# surrogate) at weights `omega` and kernel `mode`, through the distribution
# class the sampler evaluates.
.weights_children_log <- function(setup, omega, mode, what = "exact") {

  total <- 0
  for (child in seq_len(setup$children)) {
    block <- .coarse_block("cluster", mean = rep(setup$mean, setup$rows), residual_sd = setup$residual_sd,
                           loading = setup$loading, sei = rep(1, setup$rows), omega = omega,
                           z_lower = setup$z_lower, z_upper = setup$z_upper, obs_bin = setup$obs_bin[child, ],
                           x = setup$x[child, ], rules = setup$rules, mode = mode, max_points = 512,
                           scrambles = 8, initial_points = 64, seed = 9049L)
    probe <- .coarse_probe(block, max_rules = 3L)
    value <- probe[[what]]
    stopifnot(is.finite(value))
    total <- total + value
  }
  total
}


# Gauss-Legendre rule of n nodes on (0, 1) by Golub-Welsch.
.weights_legendre <- function(n) {

  k <- seq_len(n - 1L)
  off <- k / sqrt(4 * k^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(k, k + 1L)] <- off
  jacobi[cbind(k + 1L, k)] <- off
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(nodes = (decomposition$values + 1) / 2, weights = decomposition$vectors[1L, ]^2)
}


# Posterior label probabilities and weight moments of the weights fixture by
# quadrature of the exact likelihood: E[I_k] and, for branch 2, E[I_2 r] and
# E[I_2 r^2]; for branch 3, E[I_3 w2], E[I_3 w3], E[I_3 w2^2] and E[I_3 w3^2]
# (w2 = p2 + p3, w3 = p3), and with `setup$log_moments` E[I_3 log p_k] (the
# Dirichlet's sufficient statistics; log p2 is singular at t = 1, so its rule
# converges slowly: 1.7e-4 from a 192-node rule at 64 nodes, within the
# allowance, in the open step-row fixture). The exact likelihood is a numerical integral whose
# accepted rung changes with the weights (relative steps within its tolerance),
# so the references use fixed Gauss-Legendre rules (over r; over w2 and
# t = w3 / w2 with Jacobian w2) and allow the change from a rule of half the
# order.
.weights_reference <- function(setup, order = 64L, children_log = .weights_children_log) {

  branches <- setup$branches
  prior <- setup$prior
  # a common shift keeps the integrands representable
  shift <- children_log(setup, c(1, 1, 1), 0)
  one_rule <- function(n) {
    rule <- .weights_legendre(n)
    out <- list()
    for (i in seq_along(branches)) {
      branch <- branches[i]
      if (branch == 1L) {
        out[["I1"]] <- prior[i]
      } else if (branch == 2L) {
        r <- rule$nodes
        weight <- rule$weights * stats::dbeta(r, setup$beta[1], setup$beta[2]) *
          vapply(r, function(value) exp(children_log(setup, c(1, value, value), 1) - shift), numeric(1L))
        out[["I2"]] <- prior[i] * sum(weight)
        out[["I2_r"]] <- prior[i] * sum(weight * r)
        out[["I2_r2"]] <- prior[i] * sum(weight * r^2)
      } else {
        alpha <- setup$alpha
        grid <- expand.grid(i = seq_len(n), j = seq_len(n))
        w2 <- rule$nodes[grid$i]
        w3 <- w2 * rule$nodes[grid$j]
        p <- cbind(1 - w2, w2 - w3, w3)
        density <- exp(lgamma(sum(alpha)) - sum(lgamma(alpha)) + as.numeric(log(p) %*% (alpha - 1)))
        like <- vapply(seq_along(w2), function(k) {
          exp(children_log(setup, c(1, w2[k], w3[k]), 1) - shift)
        }, numeric(1L))
        weight <- rule$weights[grid$i] * rule$weights[grid$j] * w2 * density * like
        out[["I3"]] <- prior[i] * sum(weight)
        out[["I3_w2"]] <- prior[i] * sum(weight * w2)
        out[["I3_w3"]] <- prior[i] * sum(weight * w3)
        out[["I3_w22"]] <- prior[i] * sum(weight * w2^2)
        out[["I3_w32"]] <- prior[i] * sum(weight * w3^2)
        if (isTRUE(setup$log_moments)) {
          for (k in 1:3) out[[paste0("I3_logp", k)]] <- prior[i] * sum(weight * log(p[, k]))
        }
      }
    }
    values <- unlist(out)
    values / sum(values[intersect(c("I1", "I2", "I3"), names(values))])
  }
  fine <- one_rule(order)
  coarse <- one_rule(order %/% 2L)
  list(estimate = fine, allowance = abs(fine - coarse))
}


# The observables of the weights fixture from the draws of b, r, eta: means
# estimate the reference's entries (the log proportions with
# `setup$log_moments`); T = sum(eta) has the Gamma(sum alpha, 1) law under
# every label (its mean and second moment).
.weights_observables <- function(draws, setup) {

  branches <- setup$branches
  b <- draws[, "b"]
  out <- list()
  for (i in seq_along(branches)) {
    indicator <- as.numeric(b == i)
    branch <- branches[i]
    out[[paste0("I", branch)]] <- indicator
    if (branch == 2L) {
      out[["I2_r"]] <- indicator * draws[, "r"]
      out[["I2_r2"]] <- indicator * draws[, "r"]^2
    }
    if (branch == 3L) {
      eta <- draws[, paste0("eta[", 1:3, "]")]
      total <- rowSums(eta)
      w2 <- (eta[, 2] + eta[, 3]) / total
      w3 <- eta[, 3] / total
      out[["I3_w2"]] <- indicator * w2
      out[["I3_w3"]] <- indicator * w3
      out[["I3_w22"]] <- indicator * w2^2
      out[["I3_w32"]] <- indicator * w3^2
      if (isTRUE(setup$log_moments)) {
        for (k in 1:3) out[[paste0("I3_logp", k)]] <- ifelse(indicator == 1, log(eta[, k] / total), 0)
      }
      out[["T"]] <- total
      out[["T2"]] <- total^2
    }
  }
  do.call(cbind, out)
}


# Reference values of the gamma total (its law is the prior's under every label).
.weights_total_reference <- function(setup) {

  a <- sum(setup$alpha)
  c(T = a, T2 = a * (a + 1))
}


# Batch-means comparison of a chain's observables with a reference: the
# estimate, its standard error from `batches` batch means and the z-score
# (estimate - reference) / sqrt(se^2 + allowance^2).
.weights_compare <- function(observables, reference, allowance = 0 * reference, batches = 50L) {

  n <- nrow(observables)
  size <- n %/% batches
  index <- rep(seq_len(batches), each = size)
  kept <- observables[seq_len(batches * size), , drop = FALSE]
  means <- apply(kept, 2L, function(column) tapply(column, index, mean))
  estimate <- colMeans(kept)
  se <- apply(means, 2L, stats::sd) / sqrt(batches)
  names <- names(reference)
  data.frame(quantity = names, estimate = estimate[names], reference = unname(reference),
             se = se[names], z = (estimate[names] - reference) / sqrt(se[names]^2 + allowance[names]^2),
             row.names = NULL)
}


# The heterogeneity fixture: h, s ~ N+(0, s_sd^2), e1, e2 ~ Gamma(1, 1) in one
# declared block; K cluster children with residual SD sqrt(s0^2 + h s^2 w1)
# and loading sqrt(l0^2 + h s^2 w2) at a fixed mean, fixed weights (1, omega2)
# under the selection kernel.
.heterogeneity_setup <- function(children = 3L, mean = .2, s_sd = .5, omega2 = .3, s0 = .7, l0 = .6,
                                 x = NULL) {

  if (is.null(x)) {
    x <- rbind(c(1.4, .9, .3, -.2), c(1.9, .4, .7, -.6), c(2.6, 1.2, .1, -1.0),
               c(1.1, .8, .5, -.4))[seq_len(children), , drop = FALSE]
  }
  children <- nrow(x)
  rows <- ncol(x)
  obs_bin <- matrix(as.numeric(x < 0) + 1, children, rows)
  list(children = children, rows = rows, mean = mean, x = x, obs_bin = obs_bin, s_sd = s_sd,
       s0 = s0, l0 = l0, omega = c(1, omega2), z_lower = c(0, -1e300), z_upper = c(1e300, 0),
       rules = .coarse_rules(.weights_orders))
}

.heterogeneity_model <- function(setup, inits = list(h = 1L, s = .3, e = c(1, 1)), ...) {

  rows <- setup$rows
  markers <- .block_markers(c("s", "e[1]", "e[2]"))
  syntax <- paste0("model {\n",
    "  h ~ dbern(.5)\n  s ~ dnorm(0, 1 / (s_sd * s_sd))T(0,)\n",
    "  for (j in 1:2) {\n    e[j] ~ dgamma(1, 1)\n    w[j] <- e[j] / sum(e)\n  }\n",
    "  for (j in 1:", rows, ") {\n",
    "    residual_sd[j] <- sqrt(s0 * s0 + h * s * s * w[1])\n",
    "    loading[j] <- sqrt(l0 * l0 + h * s * s * w[2])\n  }\n",
    "  for (c in 1:K) {\n",
    "    y[c, 1:", rows, "] ~ dselnorm_cluster_step(means[1:", rows, "], residual_sd[1:", rows, "], loading[1:", rows,
    "], sei[1:", rows, "], omega[1:2], z_lower[1:2], z_upper[1:2], obs_bin[c, 1:", rows, "], 1, 1, 1, ",
    "nodes[1:nq], log_weights[1:nq], orders[1:nr], qmc[1:nqmc], 64, 512, 8, 0.005, 0)\n",
    "  }\n", markers$syntax, "\n}")
  design <- BayesTools::selection_qmc_design(2L, 512L, 8L, 9049L)
  data <- c(list(y = setup$x, means = rep(setup$mean, rows), sei = rep(1, rows), omega = setup$omega,
                 z_lower = setup$z_lower, z_upper = setup$z_upper, obs_bin = setup$obs_bin, K = setup$children,
                 s_sd = setup$s_sd, s0 = setup$s0, l0 = setup$l0, nodes = setup$rules$nodes,
                 log_weights = setup$rules$log_weights, orders = setup$rules$orders,
                 nq = length(setup$rules$nodes), nr = length(setup$rules$orders),
                 qmc = as.numeric(design), nqmc = length(design)), markers$data)
  .block_hand_model(syntax, data, inits = inits, ...)
}

.heterogeneity_children_log <- function(setup, variance1, variance2) {

  total <- 0
  for (child in seq_len(setup$children)) {
    block <- .coarse_block("cluster", mean = rep(setup$mean, setup$rows),
                           residual_sd = rep(sqrt(setup$s0^2 + variance1), setup$rows),
                           loading = rep(sqrt(setup$l0^2 + variance2), setup$rows), sei = rep(1, setup$rows),
                           omega = setup$omega, z_lower = setup$z_lower, z_upper = setup$z_upper,
                           obs_bin = setup$obs_bin[child, ], x = setup$x[child, ], rules = setup$rules, mode = 1,
                           max_points = 512, scrambles = 8, initial_points = 64, seed = 9049L)
    value <- .coarse_probe(block, max_rules = 3L)[["exact"]]
    stopifnot(is.finite(value))
    total <- total + value
  }
  total
}

# E[h], E[h s], E[h s^2], E[h w1], E[h w1^2], E[s] and E[w1] by fixed
# Gauss-Legendre rules over s on (0, 12 s_sd) (half-normal prior) and w1 on
# (0, 1) (the uniform Dirichlet(1, 1) law of e1 / (e1 + e2)), allowing the
# change from rules of half the order (see .weights_reference).
.heterogeneity_reference <- function(setup, order = 64L) {

  shift <- .heterogeneity_children_log(setup, 0, 0)
  upper <- 12 * setup$s_sd
  one_rule <- function(n) {
    rule <- .weights_legendre(n)
    grid <- expand.grid(i = seq_len(n), j = seq_len(n))
    s <- upper * rule$nodes[grid$i]
    w1 <- rule$nodes[grid$j]
    like <- vapply(seq_along(s), function(k) {
      exp(.heterogeneity_children_log(setup, s[k]^2 * w1[k], s[k]^2 * (1 - w1[k])) - shift)
    }, numeric(1L))
    weight <- upper * rule$weights[grid$i] * rule$weights[grid$j] * 2 * stats::dnorm(s, 0, setup$s_sd) * like
    # the h = 0 state: prior law of s and w1, likelihood one after the shift
    z1 <- .5 * sum(weight)
    total <- .5 + z1
    values <- c(h = z1, h_s = .5 * sum(weight * s), h_s2 = .5 * sum(weight * s^2),
                h_w1 = .5 * sum(weight * w1), h_w12 = .5 * sum(weight * w1^2)) / total
    prior_mean_s <- setup$s_sd * sqrt(2 / pi)
    c(values, s = values[["h_s"]] + (1 - values[["h"]]) * prior_mean_s,
      w1 = values[["h_w1"]] + (1 - values[["h"]]) * .5)
  }
  fine <- one_rule(order)
  coarse <- one_rule(order %/% 2L)
  list(estimate = fine, allowance = abs(fine - coarse))
}


.heterogeneity_observables <- function(draws) {

  h <- draws[, "h"]
  s <- draws[, "s"]
  e <- draws[, c("e[1]", "e[2]")]
  w1 <- e[, 1] / rowSums(e)
  cbind(h = h, h_s = h * s, h_s2 = h * s^2, h_w1 = h * w1, h_w12 = h * w1^2, s = s, w1 = w1,
        T = rowSums(e), T2 = rowSums(e)^2)
}


# The test hook of the corrected samplers (flags of the factories built from
# now on, per-kind counters, memo counters, live screened scores).
.corrected_test_hook <- function(flags = NULL, clear = FALSE, scores = FALSE) {

  .Call("RoBMA_corrected_sampler_test", flags, clear, scores, PACKAGE = "RoBMA")
}

# Runs `code` with the test flags set, restoring the defaults afterwards.
.with_corrected_flags <- function(flags, code) {

  .corrected_test_hook(flags)
  on.exit(.corrected_test_hook(list(disable_new = FALSE, drop_sum_jacobian = FALSE, finite_proposal = 2,
                                    force_fallback = FALSE, freeze_finite = FALSE, disable_gated = FALSE,
                                    disable_finite = FALSE, disable_sum = FALSE,
                                    memo_slots = 2)), add = TRUE)
  force(code)
}

# The widths of the live block kernels (.block_test_hook()$widths) as a list of
# (chain, d, widths).
.weights_block_widths <- function(values = .block_test_hook()$widths) {

  out <- list()
  i <- 1L
  while (i <= length(values)) {
    d <- values[i + 1L]
    out[[length(out) + 1L]] <- list(chain = values[i], d = d, widths = values[i + 1L + seq_len(d)])
    i <- i + 2L + d
  }
  out
}

# The scores of the live screened samplers (.corrected_test_hook()$screened_scores)
# as a list of (chain, surrogate, exact).
.weights_screened_scores <- function(values = .corrected_test_hook(scores = TRUE)$screened_scores) {

  out <- list()
  i <- 1L
  while (i <= length(values)) {
    levels <- values[i + 1L]
    out[[length(out) + 1L]] <- list(chain = values[i], surrogate = values[i + 1L + seq_len(levels)],
                                    exact = values[i + 1L + levels + seq_len(levels)])
    i <- i + 2L + 2L * levels
  }
  out
}
