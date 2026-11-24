# Fixtures and independent references for the heterogeneity scale move of the
# location block sampler (src/samplers/LocationScale.cc) in hand-written rjags
# models that declare the location block's and the scale's terminal tags.


# The native switches frozen into scale owners built afterwards: the mutant
# (0 none, 1 exponent c^(1 - K), 2 no a'^-(2L - 1), 3 the last latent left
# out, 4 the gamma total not preserved), the decline switch and the timing
# switch; negative values keep a setting. Returns the previous settings.
.scale_config <- function(mutant = -1L, decline = -1L, timing = -1L) {

  .Call("RoBMA_location_scale_test_config", as.integer(c(mutant, decline, timing)), PACKAGE = "RoBMA")
}


# The test hook (r-location-scale.cc): "last", "state", "set", "transform",
# "density" and "slice" for the live owner with serial number `index` and
# chain `chain` (0-based).
.scale_hook <- function(what, index = 0L, chain = 0L, values = numeric()) {

  .Call("RoBMA_location_scale_test", what, as.integer(index), as.integer(chain), as.double(values),
        PACKAGE = "RoBMA")
}


# The most recently built live scale owner.
.scale_last_owner <- function() .scale_hook("last")


# ---------------------------------------------------------------------------- #
# The allocation layout: a gated intercept, the total SD a (half-normal,
# optionally truncated above), an inclusion indicator, L gamma sources of the
# weights (w_j = eta_j / sum(eta)), the sampled component s whose SD scales
# the cluster latents in the means and the other components in the rows'
# variances (Gaussian children).
# ---------------------------------------------------------------------------- #

.scale_allocation_setup <- function(shapes = c(1, 1), s = 2L, a_sd = 0.5, a_upper = Inf, include = 0.5) {

  cluster <- rep(1:4, times = c(4L, 3L, 4L, 3L))
  list(
    y        = c(0.301, 0.224, 0.350, 0.263, -0.068, 0.026, -0.029, 0.152, 0.185, 0.114, 0.218, 0.455, 0.521,
                 0.422),
    variance = rep(c(0.02, 0.03), 7L),
    cluster  = cluster,
    shapes   = shapes,
    s        = s,
    a_sd     = a_sd,
    a_upper  = a_upper,
    include  = include,
    prior_mean = 0,
    prior_precision = 1,
    weights  = c(1, 1)
  )
}


# `extra`: an additional dependence of the rows' variances on the sampled SD
# (the invariance check's mutant); `tag`: declare the scale.
.scale_allocation_model <- function(setup, extra = FALSE, tag = TRUE) {

  L <- length(setup$shapes)
  K <- max(setup$cluster)
  others <- setdiff(seq_len(L), setup$s)
  truncation <- if (is.finite(setup$a_upper)) "T(0, a_upper)" else "T(0,)"
  integrated <- paste0("pow(sd_component[", others, "], 2)", collapse = " + ")
  if (extra) integrated <- paste(integrated, "+ 0.5 * pow(sigma, 2)")
  syntax <- c(
    "model {",
    "  g ~ dcat(gate_weights)",
    "  slab ~ dnorm(prior_mean, prior_precision)",
    "  b <- 0 * (g == 1) + slab * (g == 2)",
    paste0("  a ~ dnorm(0, a_precision) ", truncation),
    "  I ~ dbern(include)",
    "  for (j in 1:L) { eta[j] ~ dgamma(shapes[j], 1) }",
    "  for (j in 1:L) { w[j] <- eta[j] / sum(eta[1:L]) }",
    "  for (j in 1:L) { sd_component[j] <- a * I * sqrt(w[j]) }",
    paste0("  sigma <- a * I * sqrt(w[", setup$s, "])"),
    paste0("  integrated <- ", integrated),
    "  for (k in 1:K) { Z[k] ~ dnorm(0, 1) }",
    "  for (i in 1:n) {",
    "    mu[i] <- b + sigma * Z[cluster[i]]",
    "    y[i] ~ dnorm(mu[i], 1 / (variance[i] + integrated))",
    "  }",
    "  tag_group <- robma_location_group(g, slab, 2, b, 1)",
    paste0("  tag_latent_", 1:K, " <- robma_location_latent(Z[", 1:K, "], ", 1:K, ")"),
    paste0("  tag_setup <- robma_location_setup(g, 1, ", K, ")"),
    if (tag) paste0("  tag_scale <- robma_location_scale(a, a, sigma, c(", paste0("w[", 1:L, "]", collapse = ", "),
                    "), 0, I, 1, 0, c(", L, ", ", setup$s, ", ", K, "))"),
    "}")
  data <- list(y = setup$y, variance = setup$variance, cluster = setup$cluster, n = length(setup$y), K = K, L = L,
               shapes = setup$shapes, a_precision = 1 / setup$a_sd^2, include = setup$include,
               prior_mean = setup$prior_mean, prior_precision = setup$prior_precision, gate_weights = setup$weights)
  if (is.finite(setup$a_upper)) data$a_upper <- setup$a_upper
  list(syntax = paste(syntax, collapse = "\n"), data = data, layout = "allocation", L = L, K = K)
}


# The clustered layout of RoBMA(): tau and rho behind their gates (tau's
# spike at zero), tau_between scales the cluster effects gamma[k] in the
# means and tau_within enters the rows' variances.
.scale_clustered_model <- function(setup, rho_shapes = c(2, 3)) {

  K <- max(setup$cluster)
  syntax <- c(
    "model {",
    "  g ~ dcat(gate_weights)",
    "  slab ~ dnorm(prior_mean, prior_precision)",
    "  b <- 0 * (g == 1) + slab * (g == 2)",
    "  tau_indicator ~ dcat(tau_weights)",
    "  tau_component_1 <- 0",
    "  tau_component_2 ~ dnorm(0, a_precision) T(0,)",
    "  tau <- tau_component_1 * (tau_indicator == 1) + tau_component_2 * (tau_indicator == 2)",
    "  rho_indicator ~ dcat(rho_weights)",
    "  rho_component_1 ~ dbeta(rho_shapes[1], rho_shapes[2])",
    "  rho <- rho_component_1 * (rho_indicator == 1)",
    "  tau_within <- tau * sqrt(1 - rho)",
    "  tau_between <- tau * sqrt(rho)",
    "  for (k in 1:K) { gamma[k] ~ dnorm(0, 1) }",
    "  for (i in 1:n) {",
    "    mu[i] <- b + gamma[cluster[i]] * tau_between",
    "    y[i] ~ dnorm(mu[i], 1 / (variance[i] + pow(tau_within, 2)))",
    "  }",
    "  tag_group <- robma_location_group(g, slab, 2, b, 1)",
    paste0("  tag_latent_", 1:K, " <- robma_location_latent(gamma[", 1:K, "], ", 1:K, ")"),
    paste0("  tag_setup <- robma_location_setup(g, 1, ", K, ")"),
    paste0("  tag_scale <- robma_location_scale(tau_component_2, tau, tau_between, rho, rho_component_1, ",
           "c(tau_indicator, rho_indicator), c(2, 1), 0, c(2, 1, ", K, "))"),
    "}")
  data <- list(y = setup$y, variance = setup$variance, cluster = setup$cluster, n = length(setup$y), K = K,
               a_precision = 1 / setup$a_sd^2, rho_shapes = rho_shapes, tau_weights = c(1, 1), rho_weights = 1,
               prior_mean = setup$prior_mean, prior_precision = setup$prior_precision, gate_weights = setup$weights)
  list(syntax = paste(syntax, collapse = "\n"), data = data, layout = "clustered", L = 2L, K = K)
}


# Compiles a hand model (one chain, the scale switches frozen into its
# sampler) without running it; returns the rjags model, the sampler of each
# node, the owner index and the factory decisions.
.scale_compile <- function(model, inits = list(), mutant = 0L, decline = 0L, seed = 2711L) {

  BayesTools::BayesTools_load_JAGS_module(quiet = TRUE)
  previous <- .scale_config(mutant, decline, 0L)
  on.exit(.scale_config(previous[["mutant"]], previous[["decline"]], previous[["timing"]]), add = TRUE)
  RoBMA:::.location_scale_report_read()
  connection <- textConnection(model$syntax)
  on.exit(close(connection), add = TRUE)
  jags <- rjags::jags.model(connection, data = model$data, n.chains = 1L, n.adapt = 0L, quiet = TRUE,
                            inits = c(inits, list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed)))
  samplers <- rjags::list.samplers(jags)
  assigned <- unlist(lapply(seq_along(samplers), function(i) {
    stats::setNames(rep(names(samplers)[i], length(samplers[[i]])), samplers[[i]])
  }))
  report <- RoBMA:::.location_scale_report_read()
  list(jags = jags, samplers = assigned, decisions = as.data.frame(report$decisions),
       owner = if (any(report$decisions$status == "owned")) .scale_last_owner() else NA_integer_)
}


# Runs a compiled model: adapt + burnin, then draws of the monitored nodes.
.scale_run <- function(compiled, monitor, draws, burnin = 500L, adapt = 200L) {

  jags <- compiled$jags
  if (adapt > 0L) rjags::adapt(jags, adapt, end.adaptation = FALSE)
  rjags::adapt(jags, 0L, end.adaptation = TRUE)
  if (burnin > 0L) stats::update(jags, n.iter = burnin, progress.bar = "none")
  RoBMA:::.location_scale_report_read()
  values <- rjags::coda.samples(jags, monitor, n.iter = draws, progress.bar = "none")
  list(draws = as.matrix(values[[1L]]), report = RoBMA:::.location_scale_report_read())
}


# ---------------------------------------------------------------------------- #
# Independent references
# ---------------------------------------------------------------------------- #

# T_c of the move's coordinates (latents, a, weight sources) in R.
.scale_transform <- function(x, u, layout, L, K, s) {

  z <- x[seq_len(K)]
  a <- x[K + 1L]
  sources <- x[-seq_len(K + 1L)]
  c2 <- exp(2 * u)
  ws <- if (layout == "allocation") sources[s] / sum(sources) else sources[1L]
  d <- 1 - ws + c2 * ws
  moved <- if (layout == "allocation") {
    sources * ifelse(seq_along(sources) == s, c2, 1) / d
  } else {
    c2 * ws / d
  }
  c(z / exp(u), a * sqrt(d), moved)
}


# log pi of the move's coordinates (the sources' and latents' own prior
# densities; the children are unchanged along the orbit).
.scale_log_prior <- function(x, model, K) {

  z <- x[seq_len(K)]
  a <- x[K + 1L]
  sources <- x[-seq_len(K + 1L)]
  out <- sum(stats::dnorm(z, log = TRUE)) + stats::dnorm(a, 0, 1 / sqrt(model$data$a_precision), log = TRUE)
  if (model$layout == "allocation") {
    out + sum(stats::dgamma(sources, model$data$shapes, 1, log = TRUE))
  } else {
    out + stats::dbeta(sources, model$data$rho_shapes[1L], model$data$rho_shapes[2L], log = TRUE)
  }
}


# log |det J| of x -> T_u(x) by Richardson-extrapolated central differences.
.scale_log_jacobian <- function(x, u, layout, L, K, s) {

  p <- length(x)
  J <- matrix(0, p, p)
  for (i in seq_len(p)) {
    h <- 1e-3 * max(abs(x[i]), 1e-3)
    diff <- function(h) {
      up <- x; up[i] <- up[i] + h
      down <- x; down[i] <- down[i] - h
      (.scale_transform(up, u, layout, L, K, s) - .scale_transform(down, u, layout, L, K, s)) / (2 * h)
    }
    J[, i] <- (4 * diff(h / 2) - diff(h)) / 3
  }
  as.numeric(determinant(J, logarithm = TRUE)$modulus)
}


# log f(u) by Liu and Sabatti: log pi(T_u x) + log |det J_u(x)| (Haar measure
# du), relative to u = 0.
.scale_direct_density <- function(x, u, model, s) {

  K <- model$K
  value <- function(v) {
    .scale_log_prior(.scale_transform(x, v, model$layout, model$L, K, s), model, K) +
      .scale_log_jacobian(x, v, model$layout, model$L, K, s)
  }
  vapply(u, value, numeric(1)) - value(0)
}


# ---------------------------------------------------------------------------- #
# Quadrature reference of the allocation fixture: the intercept and the
# latents are Gaussian and integrated in closed form; for each gate and
# inclusion state the posterior over (a, w_s) is a two-dimensional
# Gauss-Legendre product rule.
# ---------------------------------------------------------------------------- #

.scale_gauss_legendre <- function(n, lower, upper) {

  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  b <- i / sqrt(4 * i^2 - 1)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(nodes = (upper - lower) / 2 * e$values + (upper + lower) / 2,
       weights = (upper - lower) * e$vectors[1L, ]^2)
}


.scale_allocation_reference <- function(setup, a_nodes = 120L, w_nodes = 120L, a_max = 4) {

  stopifnot(length(setup$shapes) == 2L, setup$s == 2L, !is.finite(setup$a_upper))
  n <- length(setup$y)
  H <- outer(setup$cluster, seq_len(max(setup$cluster)), `==`) * 1
  ga <- .scale_gauss_legendre(a_nodes, 0, a_max)
  gw <- .scale_gauss_legendre(w_nodes, 0, 1)
  log_prior_a <- stats::dnorm(ga$nodes, 0, setup$a_sd, log = TRUE) + log(2)
  log_prior_w <- stats::dbeta(gw$nodes, setup$shapes[2L], setup$shapes[1L], log = TRUE)
  log_gate <- log(setup$weights / sum(setup$weights))
  log_ml <- function(Sigma, gate) {
    if (gate == 2L) Sigma <- Sigma + 1 / setup$prior_precision
    mean <- if (gate == 2L) setup$prior_mean else 0
    R <- chol(Sigma)
    -sum(log(diag(R))) - 0.5 * n * log(2 * pi) - 0.5 * sum(backsolve(R, setup$y - mean, transpose = TRUE)^2)
  }
  rows <- list()
  for (gate in 1:2) {
    # heterogeneity excluded: the prior of (a, w_s)
    excluded <- log_ml(diag(setup$variance), gate)
    rows[[length(rows) + 1L]] <- list(log_mass = log_gate[gate] + log(1 - setup$include) + excluded,
                                      gate = gate, I = 0, grid = NULL)
    grid <- expand.grid(i = seq_len(a_nodes), j = seq_len(w_nodes))
    log_w <- numeric(nrow(grid))
    for (k in seq_len(nrow(grid))) {
      a  <- ga$nodes[grid$i[k]]
      ws <- gw$nodes[grid$j[k]]
      Sigma <- diag(setup$variance + a^2 * (1 - ws)) + a^2 * ws * H %*% t(H)
      log_w[k] <- log_prior_a[grid$i[k]] + log_prior_w[grid$j[k]] + log_ml(Sigma, gate) +
        log(ga$weights[grid$i[k]]) + log(gw$weights[grid$j[k]])
    }
    top <- max(log_w)
    weight <- exp(log_w - top)
    rows[[length(rows) + 1L]] <- list(log_mass = log_gate[gate] + log(setup$include) + top + log(sum(weight)),
                                      gate = gate, I = 1,
                                      grid = list(a = ga$nodes[grid$i], ws = gw$nodes[grid$j], weight = weight / sum(weight)))
  }
  log_mass <- vapply(rows, `[[`, numeric(1), "log_mass")
  p <- exp(log_mass - max(log_mass))
  p <- p / sum(p)
  # prior moments for the excluded states
  prior_a  <- c(setup$a_sd * sqrt(2 / pi), setup$a_sd^2)
  alpha <- setup$shapes[2L]
  beta  <- setup$shapes[1L]
  prior_w  <- c(alpha / (alpha + beta), alpha * (alpha + 1) / ((alpha + beta) * (alpha + beta + 1)))
  moment <- function(f_included, excluded) {
    sum(vapply(seq_along(rows), function(r) {
      row <- rows[[r]]
      p[r] * if (is.null(row$grid)) excluded else sum(row$grid$weight * f_included(row$grid$a, row$grid$ws))
    }, numeric(1)))
  }
  c(g     = sum(p[vapply(rows, `[[`, numeric(1), "gate") == 2]),
    I     = sum(p[vapply(rows, `[[`, numeric(1), "I") == 1]),
    a     = moment(function(a, w) a, prior_a[1L]),
    a2    = moment(function(a, w) a^2, prior_a[2L]),
    ws    = moment(function(a, w) w, prior_w[1L]),
    ws2   = moment(function(a, w) w^2, prior_w[2L]),
    sigma = moment(function(a, w) a * sqrt(w), 0),
    sigma2 = moment(function(a, w) a^2 * w, 0))
}


.scale_allocation_observables <- function(draws, s = 2L) {

  ws <- draws[, paste0("w[", s, "]")]
  cbind(g = draws[, "g"] == 2, I = draws[, "I"], a = draws[, "a"], a2 = draws[, "a"]^2, ws = ws, ws2 = ws^2,
        sigma = draws[, "sigma"], sigma2 = draws[, "sigma"]^2) * 1
}


# ---------------------------------------------------------------------------- #
# The selection fixture of the location block (helper-location-block.R) with
# its fixed cluster scale replaced by a random SD: sigma = a I with a
# half-normal a and an inclusion indicator I (the scale's L = 1 layout), and
# the gated intercept as the only location term. Rank-one cluster selection
# children on the crude coarse schedule; the source a is updated by the
# corrected scalar sampler embedded in the scale owner.
# ---------------------------------------------------------------------------- #

.scale_selection_model <- function(setup, a_sd = 0.3, include = 0.5) {

  clusters <- split(seq_along(setup$y), setup$cluster)
  child <- function(s) {
    n <- length(clusters[[s]])
    paste0("  y", s, "[1:", n, "] ~ dselnorm_cluster_step(m", s, "[1:", n, "], rsd", s, "[1:", n, "], ld", s,
           "[1:", n, "], se", s, "[1:", n, "], omega[1:2], z_lower[1:2], z_upper[1:2], bin", s, "[1:", n,
           "], 1, 1, 1, nodes[1:nq], log_weights[1:nq], orders[1:nr], qmc[1:nqmc], 64, 512, 8, 0.005, 0)")
  }
  means <- unlist(lapply(seq_along(clusters), function(s) {
    paste0("  m", s, "[", seq_along(clusters[[s]]), "] <- b + sigma * Z[", s, "]")
  }))
  syntax <- c(
    "model {",
    "  g ~ dcat(gate_weights)",
    "  slab ~ dnorm(prior_mean, prior_precision)",
    "  b <- 0 * (g == 1) + slab * (g == 2)",
    "  a ~ dnorm(0, a_precision) T(0,)",
    "  I ~ dbern(include)",
    "  sigma <- a * I",
    "  for (s in 1:3) { Z[s] ~ dnorm(0, 1) }",
    means,
    vapply(seq_along(clusters), child, character(1)),
    "  tag_group <- robma_location_group(g, slab, 2, b, 1)",
    paste0("  tag_latent_", 1:3, " <- robma_location_latent(Z[", 1:3, "], ", 1:3, ")"),
    "  tag_setup <- robma_location_setup(g, 1, 3)",
    "  tag_scale <- robma_location_scale(a, a, sigma, 0, 0, I, 1, 0, c(1, 1, 3))",
    "}")
  design <- BayesTools::selection_qmc_design(2L, 512L, 8L, 9049L)
  data <- list(prior_mean = setup$prior_mean[[1L]], prior_precision = setup$prior_precision[[1L]],
               gate_weights = setup$weights, a_precision = 1 / a_sd^2, include = include, omega = setup$omega,
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
  list(syntax = paste(syntax, collapse = "\n"), data = data, qmc = design, a_sd = a_sd, include = include)
}


# The reference: per gate and inclusion state, the intercept on a
# Gauss-Hermite rule adapted to its conditional mode and curvature at each SD
# node, the SD on a Gauss-Legendre rule against its half-normal prior and the
# cluster effects by Gauss-Hermite rules against their N(0, 1) priors, with
# the exact cluster likelihood of the package. Moments of g, I, b, a and
# sigma = a I.
.scale_selection_reference <- function(setup, model, a_nodes = 40L, b_nodes = 20L, latent_nodes = 40L,
                                       a_max = 1.5) {

  latent <- .coarse_reference_rule(latent_nodes)
  rule <- .coarse_reference_rule(b_nodes)
  ga <- .scale_gauss_legendre(a_nodes, 0, a_max)
  # the log likelihood at intercepts b (vector) and SD sigma, each cluster
  # effect integrated by a Gauss-Hermite rule centred and scaled by its
  # Gaussian conditional without selection (cluster covariance diag(rsd^2) +
  # ld ld'), the N(0, 1) prior and the rule's weight function divided out
  log_lik <- function(b, sigma) {
    total <- numeric(length(b))
    for (s in 1:3) {
      rows   <- which(setup$cluster == s)
      ones   <- solve(diag(setup$residual_sd[rows]^2, length(rows)) + tcrossprod(setup$loading[rows]),
                      rep(1, length(rows)))
      lambda <- 1 + sigma^2 * sum(ones)
      center <- sigma * (sum(ones * setup$y[rows]) - b * sum(ones)) / lambda
      scale  <- 1 / sqrt(lambda)
      z      <- outer(latent$nodes, rep(scale, length(b))) + rep(center, each = latent_nodes)
      point  <- rep(b, each = latent_nodes) + sigma * as.vector(z)
      means  <- matrix(point, nrow = length(point), ncol = length(rows))
      ll <- matrix(.location_selection_cluster_loglik(setup, model, s, means), nrow = latent_nodes) +
        stats::dnorm(z, log = TRUE) - stats::dnorm(latent$nodes, log = TRUE) + log(scale)
      top <- apply(ll, 2L, max)
      total <- total + top + log(colSums(latent$weights * exp(sweep(ll, 2L, top))))
    }
    total
  }
  log_prior_b <- function(b) {
    stats::dnorm(b, setup$prior_mean[[1L]], 1 / sqrt(setup$prior_precision[[1L]]), log = TRUE)
  }
  # the intercept integrated at one SD: the log mass and the first two moments
  integrate_b <- function(sigma) {
    f <- function(b) log_lik(b, sigma) + log_prior_b(b)
    mode <- stats::optimize(f, c(-2, 2), maximum = TRUE, tol = 1e-8)$maximum
    h <- 1e-2
    curvature <- -(f(mode + h) - 2 * f(mode) + f(mode - h)) / h^2
    # a flat or irregular conditional (negligible mass) keeps the prior scale
    sd <- if (is.finite(curvature) && curvature > 0) 1 / sqrt(curvature) else 1 / sqrt(setup$prior_precision[[1L]])
    b <- mode + sd * rule$nodes
    log_f <- log_lik(b, sigma) + log_prior_b(b) + 0.5 * rule$nodes^2
    top <- max(log_f)
    weight <- rule$weights * exp(log_f - top)
    mass <- sum(weight)
    list(log_mass = top + log(mass) + log(sd) + 0.5 * log(2 * pi), m1 = sum(weight * b) / mass,
         m2 = sum(weight * b^2) / mass)
  }
  log_gate <- log(setup$weights / sum(setup$weights))
  log_prior_a <- stats::dnorm(ga$nodes, 0, model$a_sd, log = TRUE) + log(2)
  states <- list()
  for (I in 0:1) {
    sigmas <- if (I == 1L) ga$nodes else 0
    for (gate in 1:2) {
      parts <- lapply(sigmas, function(sigma) {
        if (gate == 1L) list(log_mass = log_lik(0, sigma), m1 = 0, m2 = 0) else integrate_b(sigma)
      })
      log_mass <- vapply(parts, `[[`, numeric(1), "log_mass")
      if (I == 1L) log_mass <- log_mass + log_prior_a + log(ga$weights)
      top <- max(log_mass)
      w <- exp(log_mass - top)
      a_moments <- if (I == 1L) c(sum(w * ga$nodes), sum(w * ga$nodes^2)) / sum(w) else
        c(model$a_sd * sqrt(2 / pi), model$a_sd^2)
      states[[length(states) + 1L]] <- list(
        log_mass = log_gate[gate] + log(if (I == 1L) model$include else 1 - model$include) + top + log(sum(w)),
        gate = gate, I = I, a = a_moments, sigma = if (I == 1L) a_moments else c(0, 0),
        b = c(sum(w * vapply(parts, `[[`, numeric(1), "m1")), sum(w * vapply(parts, `[[`, numeric(1), "m2"))) /
          sum(w))
    }
  }
  log_mass <- vapply(states, `[[`, numeric(1), "log_mass")
  p <- exp(log_mass - max(log_mass))
  p <- p / sum(p)
  pick <- function(name, k) sum(p * vapply(states, function(state) state[[name]][[k]], numeric(1)))
  c(g = sum(p[vapply(states, `[[`, numeric(1), "gate") == 2]), I = sum(p[vapply(states, `[[`, numeric(1), "I") == 1]),
    b = pick("b", 1L), b2 = pick("b", 2L), a = pick("a", 1L), a2 = pick("a", 2L), sigma = pick("sigma", 1L),
    sigma2 = pick("sigma", 2L))
}


.scale_selection_observables <- function(draws) {

  cbind(g = draws[, "g"] == 2, I = draws[, "I"], b = draws[, "b"], b2 = draws[, "b"]^2, a = draws[, "a"],
        a2 = draws[, "a"]^2, sigma = draws[, "sigma"], sigma2 = draws[, "sigma"]^2) * 1
}
