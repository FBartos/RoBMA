# The bias-label fixture of the location block sampler (RoBMA::LocationBlock
# with an owned bias indicator, src/samplers/LocationBlockBias.h) and its
# quadrature reference: a hand-written rjags model with three correlated
# clusters (rank-one cluster selection blocks, a crude coarse schedule), a
# gated intercept with a nonzero slab mean, a gated between-cluster
# moderator, standardized cluster effects with a fixed scale and a bias
# mixture of four labels (none, one-sided step selection, PET with a
# half-Cauchy slope of scale 0.5, PEESE with a Gamma(2, 2) slope). The nested
# variant gives the step weight a beta prior whose parameters depend on the
# label (the nested weight node of BayesTools' nested layout).


.bias_setup <- function(nested = FALSE, rho = .85, second_selection = FALSE) {

  se <- c(.20, .22, .25, .21, .18, .24, .19, .23, .20, .22)
  list(
    y           = c(.25, .32, .18, .40, .45, .30, .52, .12, .28, .22),
    se          = se,
    residual_sd = se * sqrt(1 - rho),
    loading     = se * sqrt(rho),
    cluster     = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L),
    x1          = c(0, 0, 0, 0, 1, 1, 1, .5, .5, .5),
    scale       = .15,
    prior_mean  = c(.2, 0),
    prior_precision = c(4, 4),
    gate_weights = c(1, 1),
    bias_weights = c(1, 2, 1, 1, 1.5),  # none, selection, PET, PEESE (, the second selection label)
    omega       = .3,               # the step weight below zero (ordinary layout)
    nested      = nested,
    nested_shape = c(3, 2),         # the step weight's beta prior at the selection label
    pet_scale   = .5,
    peese_shape = c(2, 2),
    # variants of the frozen-state tests: the PET prior's JAGS text, a model
    # without the PEESE label (three labels), an ungated intercept
    pet_prior   = "dt(0, pet_precision, 1) T(0, )",
    peese       = TRUE,
    gated_intercept = TRUE,
    # a second step-selection label (label 5) with its own weight: labels that
    # change at identical means, two selection and three Gaussian labels
    second_selection = second_selection,
    omega_second = .8,
    rules       = .coarse_rules(c(1L, 2L, 3L, 7L, 15L, 31L, 63L, 127L, 255L, 511L, 1023L))
  )
}


.bias_model <- function(setup) {

  clusters <- split(seq_along(setup$y), setup$cluster)
  intercept <- if (setup$gated_intercept) "b[1]" else "b0"
  labels <- if (setup$second_selection) 5L else if (setup$peese) 4L else 3L
  codes  <- c(0, 1, 2, 3, 1)[seq_len(labels)]
  child <- function(s) {
    n <- length(clusters[[s]])
    paste0("  y", s, "[1:", n, "] ~ dselnorm_cluster_step(m", s, "[1:", n, "], rsd", s, "[1:", n, "], ld", s,
           "[1:", n, "], se", s, "[1:", n, "], omega[1:2], z_lower[1:2], z_upper[1:2], bin", s, "[1:", n,
           "], 1, 1, kernel, nodes[1:nq], log_weights[1:nq], orders[1:nr], qmc[1:nqmc], 64, 512, 8, 0.005, 0)")
  }
  means <- unlist(lapply(seq_along(clusters), function(s) {
    rows <- clusters[[s]]
    paste0("  m", s, "[", seq_along(rows), "] <- ", intercept, " + b[2] * x1[", rows, "] + scale * Z[", s,
           "] + PET * pet_col[", rows, "]", if (setup$peese) paste0(" + PEESE * peese_col[", rows, "]") else "")
  }))
  weight <- if (setup$nested) {
    c("  a ~ dbeta(1 + (shape[1] - 1) * equals(bias_indicator, 2), 1 + (shape[2] - 1) * equals(bias_indicator, 2))",
      "  omega[2] <- a")
  } else if (setup$second_selection) {
    "  omega[2] <- omega_fixed * equals(bias_indicator, 2) + omega_second * equals(bias_indicator, 5)"
  } else {
    "  omega[2] <- omega_fixed"
  }
  syntax <- c(
    "model {",
    paste0("  for (j in ", if (setup$gated_intercept) 1L else 2L, ":2) {"),
    "    g[j] ~ dcat(gate_weights)",
    "    slab[j] ~ dnorm(prior_mean[j], prior_precision[j])",
    "    b[j] <- 0 * (g[j] == 1) + slab[j] * (g[j] == 2)",
    "  }",
    if (!setup$gated_intercept) "  b0 ~ dnorm(prior_mean[1], prior_precision[1])",
    "  for (s in 1:3) { Z[s] ~ dnorm(0, 1) }",
    "  bias_indicator ~ dcat(bias_weights)",
    paste0("  PET_1 ~ ", setup$pet_prior),
    "  PET <- PET_1 * equals(bias_indicator, 3)",
    if (setup$peese) c("  PEESE_1 ~ dgamma(peese_shape[1], peese_shape[2])",
                        "  PEESE <- PEESE_1 * equals(bias_indicator, 4)"),
    if (setup$second_selection) "  kernel <- equals(bias_indicator, 2) + equals(bias_indicator, 5)" else
      "  kernel <- equals(bias_indicator, 2)",
    "  omega[1] <- 1",
    weight,
    means,
    vapply(seq_along(clusters), child, character(1)),
    if (setup$gated_intercept) "  tag_group_1 <- robma_location_group(g[1], slab[1], 2, b[1], 1)" else
      "  tag_group_1 <- robma_location_group(b0, b0, 0, b0, 1)",
    "  tag_group_2 <- robma_location_group(g[2], slab[2], 2, b[2], 2)",
    paste0("  tag_latent_", 1:3, " <- robma_location_latent(Z[", 1:3, "], ", 1:3, ")"),
    paste0("  tag_setup <- robma_location_setup(", if (setup$gated_intercept) "g[1]" else "b0", ", 2, 3)"),
    paste0("  tag_bias <- robma_location_bias(bias_indicator, codes[1:", labels, "], PET_1, ",
           if (setup$peese) "PEESE_1" else "0", ", 1)"),
    "}")
  design <- BayesTools::selection_qmc_design(2L, 512L, 8L, 9049L)
  data <- list(x1 = setup$x1, pet_col = setup$se,
               peese_col = if (is.null(setup$peese_column)) setup$se^2 else setup$peese_column, scale = setup$scale,
               prior_mean = setup$prior_mean, prior_precision = setup$prior_precision,
               gate_weights = setup$gate_weights, bias_weights = setup$bias_weights[seq_len(labels)],
               codes = codes, pet_precision = 1 / setup$pet_scale^2,
               z_lower = c(0, -1e300), z_upper = c(1e300, 0), qmc = as.numeric(design), nqmc = length(design),
               nodes = setup$rules$nodes, log_weights = setup$rules$log_weights, orders = setup$rules$orders,
               nq = length(setup$rules$nodes), nr = length(setup$rules$orders))
  if (setup$nested) data$shape <- setup$nested_shape else data$omega_fixed <- setup$omega
  if (setup$second_selection) data$omega_second <- setup$omega_second
  if (setup$peese) data$peese_shape <- setup$peese_shape else data$peese_col <- NULL
  if (!grepl("pet_precision", setup$pet_prior, fixed = TRUE)) data$pet_precision <- NULL
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


# Gauss-Legendre rule on (0, 1) (Golub-Welsch).
.bias_legendre <- function(n) {

  k <- seq_len(n - 1L)
  beta <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- beta
  J[cbind(k + 1L, k)] <- beta
  e <- eigen(J, symmetric = TRUE)
  list(nodes = (e$values + 1) / 2, weights = e$vectors[1L, ]^2)
}


# The normalized slope priors (PDF on [0, Inf)).
.bias_slope_density <- function(setup, which, x) {

  if (which == "PET") return(2 * stats::dcauchy(x, 0, setup$pet_scale))
  stats::dgamma(x, setup$peese_shape[1], setup$peese_shape[2])
}


# The reference: per label and gate pattern the log mass and the first two
# moments of the physical intercept, and per slope label the first two
# moments of the slope. Gaussian labels in closed form (the cluster effects
# and active coefficients integrated, the slope by adaptive quadrature against
# its prior); the selection label per gate pattern by the Laplace-adapted
# Gauss-Hermite grid of the selection fixture with the package's exact
# cluster likelihood (and, nested, a Gauss-Legendre rule in the step weight).
.bias_reference <- function(setup, model, coefficient_nodes = 16L, latent_nodes = 32L, weight_nodes = 16L) {

  n <- length(setup$y)
  X <- cbind(1, setup$x1)
  H <- outer(setup$cluster, 1:3, `==`) * 1
  Sigma <- diag(setup$residual_sd^2) + H %*% t(H) * setup$scale^2
  for (s in 1:3) {
    rows <- which(setup$cluster == s)
    Sigma[rows, rows] <- Sigma[rows, rows] + tcrossprod(setup$loading[rows])
  }
  patterns <- as.matrix(expand.grid(g1 = 0:1, g2 = 0:1))
  log_gate <- log(setup$gate_weights / sum(setup$gate_weights))
  n_labels <- if (setup$second_selection) 5L else 4L
  weights  <- setup$bias_weights[seq_len(n_labels)]
  log_bias <- log(weights / sum(weights))
  slopes <- list(PET = setup$se, PEESE = setup$se^2)

  # a Gaussian label and gate pattern: log mass and moments of b0 (and the slope)
  gaussian <- function(on, slope) {
    m0 <- setup$prior_mean[on]
    V  <- diag(1 / setup$prior_precision[on], nrow = length(on))
    XA <- X[, on, drop = FALSE]
    C  <- Sigma + XA %*% V %*% t(XA)
    Ci <- solve(C)
    mu0 <- drop(XA %*% m0)
    r0 <- setup$y - mu0
    log_det <- as.numeric(determinant(C)$modulus)
    k <- -0.5 * (n * log(2 * pi) + log_det + drop(t(r0) %*% Ci %*% r0))
    # E[b | y, theta] = m0 + V XA' Ci (r0 - theta s); Var = V - V XA' Ci XA V
    gain <- V %*% t(XA) %*% Ci
    post_var <- V - gain %*% XA %*% V
    i0 <- match(1L, on)
    if (is.null(slope)) {
      mean_b0 <- if (is.na(i0)) 0 else drop(m0 + gain %*% r0)[i0]
      var_b0  <- if (is.na(i0)) 0 else post_var[i0, i0]
      return(list(log_mass = k, b0 = mean_b0, b0_2 = var_b0 + mean_b0^2, t1 = 0, t2 = 0))
    }
    s <- slopes[[slope]]
    q <- drop(t(s) %*% Ci %*% s)
    r <- drop(t(s) %*% Ci %*% r0)
    # integrals of theta^j p(theta) exp(r theta - q theta^2 / 2) on [0, Inf), scaled at the mode
    top <- if (r > 0) r^2 / (2 * q) else 0
    moment <- function(j) {
      stats::integrate(function(t) t^j * .bias_slope_density(setup, slope, t) * exp(r * t - q * t^2 / 2 - top),
                       0, Inf, rel.tol = 1e-11, subdivisions = 1000L)$value
    }
    I0 <- moment(0); I1 <- moment(1); I2 <- moment(2)
    t1 <- I1 / I0; t2 <- I2 / I0
    if (is.na(i0)) {
      b0 <- 0; b0_2 <- 0
    } else {
      a0 <- drop(m0 + gain %*% r0)[i0]
      a1 <- -drop(gain %*% s)[i0]
      b0 <- a0 + a1 * t1
      b0_2 <- post_var[i0, i0] + a0^2 + 2 * a0 * a1 * t1 + a1^2 * t2
    }
    list(log_mass = k + top + log(I0), b0 = b0, b0_2 = b0_2, t1 = t1, t2 = t2)
  }

  # the selection label at one step weight: per gate pattern
  latent <- .coarse_reference_rule(latent_nodes)
  selection <- function(on, omega, anchor = NULL) {
    log_marginal <- function(B) {
      total <- numeric(nrow(B))
      for (s in 1:3) {
        rows  <- which(setup$cluster == s)
        base  <- B %*% t(X[rows, , drop = FALSE])
        means <- base[rep(seq_len(nrow(B)), each = latent_nodes), , drop = FALSE] +
          setup$scale * latent$nodes[rep(seq_len(latent_nodes), nrow(B))]
        k <- nrow(means)
        ll <- .Call("RoBMA_selnorm_cluster_step_loglik_batch", setup$y[rows], means,
                    matrix(rep(setup$residual_sd[rows], each = k), k), matrix(rep(setup$loading[rows], each = k), k),
                    setup$se[rows], matrix(rep(c(1, omega), each = k), k), c(0, -Inf), c(Inf, 0),
                    ifelse(setup$y[rows] > 0, 1L, 2L), 1L, TRUE, 1L, setup$rules$nodes, setup$rules$log_weights,
                    as.double(setup$rules$orders), as.double(model$qmc), 64L, 512L, 8L, .005, FALSE, 0L,
                    PACKAGE = "RoBMA")$log_density
        ll  <- matrix(ll, nrow = latent_nodes)
        top <- apply(ll, 2L, max)
        total <- total + top + log(colSums(latent$weights * exp(sweep(ll, 2L, top))))
      }
      total
    }
    if (length(on) == 0L) return(list(log_mass = log_marginal(matrix(0, 1, 2)), b0 = 0, b0_2 = 0))
    full <- function(beta) {
      beta <- matrix(beta, ncol = length(on))
      B <- matrix(0, nrow(beta), 2)
      B[, on] <- beta
      prior <- rep(0, nrow(beta))
      for (i in seq_along(on)) {
        prior <- prior + stats::dnorm(beta[, i], setup$prior_mean[on[i]], 1 / sqrt(setup$prior_precision[on[i]]),
                                      log = TRUE)
      }
      log_marginal(B) + prior
    }
    if (is.null(anchor)) {
      fit <- stats::optim(setup$prior_mean[on], function(beta) -full(beta), method = "BFGS", hessian = TRUE)
      # the grid's scale from the curvature at the mode (the prior's where it
      # is not positive definite)
      hessian <- fit$hessian
      if (!all(is.finite(hessian)) || any(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        hessian <- diag(setup$prior_precision[on], nrow = length(on))
      }
      anchor <- list(par = fit$par, L = t(chol(solve(hessian))))
    }
    fit <- list(par = anchor$par)
    L <- anchor$L
    rule <- .coarse_reference_rule(coefficient_nodes)
    grid <- as.matrix(expand.grid(rep(list(seq_len(coefficient_nodes)), length(on))))
    x <- matrix(rule$nodes[grid], ncol = length(on))
    w <- apply(matrix(rule$weights[grid], ncol = length(on)), 1L, prod)
    beta <- sweep(x %*% t(L), 2L, fit$par, `+`)
    log_f <- full(beta) + 0.5 * rowSums(x^2)
    top <- max(log_f)
    weight <- w * exp(log_f - top)
    mass <- sum(weight)
    i0 <- match(1L, on)
    b0 <- if (is.na(i0)) 0 else sum(weight * beta[, i0]) / mass
    b0_2 <- if (is.na(i0)) 0 else sum(weight * beta[, i0]^2) / mass
    list(log_mass = top + log(mass) + sum(log(diag(L))) + 0.5 * length(on) * log(2 * pi), b0 = b0, b0_2 = b0_2,
         anchor = anchor)
  }

  rows <- list()
  for (label in seq_len(n_labels)) {
    for (k in seq_len(nrow(patterns))) {
      on <- which(patterns[k, ] == 1)
      base <- log_bias[label] + sum(log_gate[patterns[k, ] + 1])
      cell <- if (label %in% c(2L, 5L)) {
        if (!setup$nested) {
          selection(on, if (label == 2L) setup$omega else setup$omega_second)
        } else {
          # the step weight integrated against its beta prior at the selection
          # label, on one coefficient grid adapted at the prior mean weight (a
          # refit per weight node makes the rule jump between node counts)
          rule <- .bias_legendre(weight_nodes)
          centre <- selection(on, setup$nested_shape[1] / sum(setup$nested_shape))
          cells <- lapply(rule$nodes, function(a) selection(on, a, centre$anchor))
          lm <- vapply(cells, `[[`, numeric(1), "log_mass") +
            stats::dbeta(rule$nodes, setup$nested_shape[1], setup$nested_shape[2], log = TRUE) + log(rule$weights)
          top <- max(lm)
          w <- exp(lm - top)
          list(log_mass = top + log(sum(w)), b0 = sum(w * vapply(cells, `[[`, numeric(1), "b0")) / sum(w),
               b0_2 = sum(w * vapply(cells, `[[`, numeric(1), "b0_2")) / sum(w))
        }
      } else {
        slope <- c(NA, NA, "PET", "PEESE")[label]
        gaussian(on, if (is.na(slope)) NULL else slope)
      }
      rows[[length(rows) + 1L]] <- data.frame(label = label, g1 = patterns[k, 1], g2 = patterns[k, 2],
                                              log_mass = base + cell$log_mass, b0 = cell$b0, b0_2 = cell$b0_2,
                                              t1 = if (is.null(cell$t1)) 0 else cell$t1,
                                              t2 = if (is.null(cell$t2)) 0 else cell$t2)
    }
  }
  rows <- do.call(rbind, rows)
  p <- exp(rows$log_mass - max(rows$log_mass))
  p <- p / sum(p)
  out <- c()
  for (label in seq_len(n_labels)) out[paste0("I", label)] <- sum(p[rows$label == label])
  for (label in seq_len(n_labels)) out[paste0("g1_", label)] <- sum(p[rows$label == label & rows$g1 == 1])
  for (label in seq_len(n_labels)) out[paste0("g2_", label)] <- sum(p[rows$label == label & rows$g2 == 1])
  out["pet"]    <- sum((p * rows$t1)[rows$label == 3])
  out["pet_2"]  <- sum((p * rows$t2)[rows$label == 3])
  out["peese"]  <- sum((p * rows$t1)[rows$label == 4])
  out["peese_2"] <- sum((p * rows$t2)[rows$label == 4])
  for (label in seq_len(n_labels)) out[paste0("b0_", label)] <- sum((p * rows$b0)[rows$label == label])
  for (label in seq_len(n_labels)) out[paste0("b0sq_", label)] <- sum((p * rows$b0_2)[rows$label == label])
  out
}


.bias_observables <- function(draws, n_labels = 4L) {

  label <- draws[, "bias_indicator"]
  g1 <- draws[, "g[1]"] == 2
  g2 <- draws[, "g[2]"] == 2
  b0 <- draws[, "b[1]"]
  per_label <- function(prefix, x) {
    out <- sapply(seq_len(n_labels), function(l) x * (label == l))
    colnames(out) <- paste0(prefix, seq_len(n_labels))
    out
  }
  out <- cbind(per_label("I", 1), per_label("g1_", g1), per_label("g2_", g2),
               pet = draws[, "PET"], pet_2 = draws[, "PET"]^2, peese = draws[, "PEESE"], peese_2 = draws[, "PEESE"]^2,
               per_label("b0_", b0), per_label("b0sq_", b0^2))
  out * 1
}
