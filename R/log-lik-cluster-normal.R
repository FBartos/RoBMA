# ============================================================================ #
# Normal Cluster-Unit Log-Likelihoods
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# .log_lik_cluster_norm_analytic
# ---------------------------------------------------------------------------- #
#
# Closed-form cluster-unit likelihood for normal multilevel models without
# selection, weighted or not. The normal cluster effect is integrated in the
# unit-normal latent by the stable form of .log_lik_cluster_norm_stable():
# every row has the standard deviation sqrt(sei^2 + tau_within^2) and the
# loading tau_between, and a row with likelihood power w contributes its
# Gaussian density to the power w.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_norm_analytic <- function(setup, yi, sei) {

  if (.has_native_norm_cluster_analytic()) {
    native <- .log_lik_cluster_norm_analytic_native(
      setup   = setup,
      yi      = yi,
      sei     = sei,
      row_sum = FALSE
    )
    if (all(is.finite(native))) {
      return(native)
    }
  }

  return(.log_lik_cluster_norm_analytic_r(
    setup = setup,
    yi    = yi,
    sei   = sei
  ))
}


.log_lik_cluster_norm_analytic_sum <- function(setup, yi, sei) {

  if (.has_native_norm_cluster_analytic(row_sum = TRUE)) {
    native <- .log_lik_cluster_norm_analytic_native(
      setup   = setup,
      yi      = yi,
      sei     = sei,
      row_sum = TRUE
    )
    if (all(is.finite(native))) {
      return(native)
    }
  }

  log_lik <- .log_lik_cluster_norm_analytic_r(
    setup = setup,
    yi    = yi,
    sei   = sei
  )
  return(.rowSums(log_lik, nrow(log_lik), ncol(log_lik)))
}


.log_lik_cluster_norm_analytic_native <- function(setup, yi, sei,
                                                   row_sum = FALSE) {

  cluster <- .cluster_indices_flatten(setup[["cluster"]])
  symbol  <- if (isTRUE(row_sum)) {
    "RoBMA_norm_cluster_analytic_loglik_row_sum"
  } else {
    "RoBMA_norm_cluster_analytic_loglik"
  }

  return(.Call(
    symbol,
    .native_numeric_vector(yi),
    .native_numeric_vector(sei),
    .native_numeric_matrix(setup[["mu"]]),
    .native_numeric_matrix(setup[["tau_within"]]),
    .native_numeric_matrix(setup[["tau_between"]]),
    cluster[["index"]],
    cluster[["size"]],
    .log_lik_cluster_norm_native_weights(setup),
    PACKAGE = "RoBMA"
  ))
}


.log_lik_cluster_norm_native_weights <- function(setup) {

  if (is.null(setup[["weights"]])) {
    return(NULL)
  }
  return(.native_numeric_vector(setup[["weights"]]))
}


.log_lik_cluster_norm_analytic_rho_grid_sum <- function(setup, yi, sei,
                                                         rho) {

  if (!is.loaded(
    "RoBMA_norm_cluster_analytic_rho_grid_loglik",
    PACKAGE = "RoBMA"
  )) {
    return(NULL)
  }
  cluster <- .cluster_indices_flatten(setup[["cluster"]])

  return(.Call(
    "RoBMA_norm_cluster_analytic_rho_grid_loglik",
    .native_numeric_vector(yi),
    .native_numeric_vector(sei),
    .native_numeric_matrix(setup[["mu"]]),
    .native_numeric_matrix(setup[["tau_total"]]),
    .native_numeric_vector(rho),
    cluster[["index"]],
    cluster[["size"]],
    .log_lik_cluster_norm_native_weights(setup),
    PACKAGE = "RoBMA"
  ))
}




.log_lik_cluster_norm_analytic_r <- function(setup, yi, sei) {

  cluster <- .cluster_indices_flatten(setup[["cluster"]])
  rows    <- cluster[["index"]]
  group   <- rep.int(seq_along(cluster[["size"]]), cluster[["size"]])
  weights <- setup[["weights"]][rows]
  S       <- setup[["S"]]
  log_lik <- matrix(NA_real_, nrow = S, ncol = length(cluster[["size"]]))

  for (s in seq_len(S)) {
    log_lik[s, ] <- .log_lik_cluster_norm_stable(
      residual = yi[rows] - setup[["mu"]][s, rows],
      sd       = .root_sum_squares(sei[rows], setup[["tau_within"]][s, rows]),
      loading  = setup[["tau_between"]][s, rows],
      weights  = weights,
      group    = group
    )
  }

  return(log_lik)
}



# ---------------------------------------------------------------------------- #
# .log_lik_cluster_norm_stable
# ---------------------------------------------------------------------------- #
#
# Log of the powered Gaussian cluster likelihood integrated over one shared
# unit-normal latent g per cluster,
#
#   log int prod_i N(residual_i; loading_i * g, sd_i)^weight_i N(g; 0, 1) dg,
#
# for the clusters named by `group` (integer ids 1..G; one value per cluster is
# returned). The vectors hold one entry per estimate, `weights` is NULL for unit
# powers, and entries with weight 0 are omitted. The loading is per estimate
# (a scale model gives each estimate its own). Over the loaded entries
# (loading != 0), log b_i = log|loading_i| - log(sd_i) is carried on the log
# scale with log(beta) = max log(b_i), omega_i = weight_i exp(2 (log b_i -
# log(beta))), W = sum(omega_i), and rho is the omega-weighted mean of the
# implied latent values residual_i / loading_i after one correction step. With
# log P = 2 log(beta) + log(W) and A = 1 + P, the Gaussian is completed in the
# unit-normal latent:
#
#   e_i = ((r_i - loading_i rho) + loading_i rho / A) / sd_i,
#   log L = -sum(weight (log(2 pi) + 2 log(sd))) / 2 - log(A) / 2
#           - (sum(weight e_i^2) + (P rho / A)^2) / 2,
#
# with e_i = r_i / sd_i for an unloaded entry. The residual after the latent
# projection is formed directly, so no explained quadratic is subtracted from
# the total, and no square of a scale is formed before it is needed. A cluster
# with a non-finite value or a non-positive sd gives NA; a result beyond the
# representable range, including an implied latent value residual / loading
# that overflows, is -Inf.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_norm_stable <- function(residual, sd, loading, weights = NULL,
                                         group = rep.int(1L, length(residual))) {

  n_group <- if (length(group)) max(group) else 0L
  out     <- numeric(n_group)
  if (is.null(weights)) {
    weights <- rep.int(1, length(residual))
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    return(rep(NA_real_, n_group))
  }

  active  <- weights > 0
  invalid <- active & (!is.finite(residual) | !is.finite(sd) | sd <= 0 |
                         !is.finite(loading))
  failed  <- tabulate(group[invalid], nbins = n_group) > 0L
  use     <- active & !failed[group]
  out[failed] <- NA_real_
  if (!any(use)) {
    return(out)
  }
  r <- residual[use]
  s <- sd[use]
  t <- loading[use]
  w <- weights[use]
  g <- group[use]

  norm  <- -0.5 * .group_sum(w * (log(2 * pi) + 2 * log(s)), g, n_group)
  value <- norm - 0.5 * .group_sum(w * (r / s)^2, g, n_group)
  value[failed] <- NA_real_

  loaded <- t != 0
  if (any(loaded)) {
    g_l      <- g[loaded]
    log_b    <- log(abs(t[loaded])) - log(s[loaded])
    log_beta <- .group_max(log_b, g_l, n_group)
    omega    <- w[loaded] * exp(2 * (log_b - log_beta[g_l]))
    total    <- .group_sum(omega, g_l, n_group)
    implied  <- r[loaded] / t[loaded]
    rho      <- .group_sum(omega * implied, g_l, n_group) / total
    rho      <- rho + .group_sum(omega * (implied - rho[g_l]), g_l, n_group) / total

    log_p     <- 2 * log_beta + log(total)
    log_a     <- ifelse(log_p > 0, log_p + log1p(exp(-log_p)),
                        log1p(exp(log_p)))
    inverse_a <- exp(-log_a)
    latent_rho <- exp(log_p - log_a) * rho
    e <- r / s
    e[loaded] <- ((r[loaded] - t[loaded] * rho[g_l]) +
                    t[loaded] * rho[g_l] * inverse_a[g_l]) / s[loaded]
    latent <- norm - 0.5 * log_a -
      0.5 * (.group_sum(w * e^2, g, n_group) + latent_rho^2)

    informative <- !failed & is.finite(log_beta)
    value[informative] <- latent[informative]
    # A quotient r / loading beyond the double range (finite inputs, a tiny
    # loading) leaves the implied latent value, and so the result,
    # unrepresentable: -Inf, never NA.
    value[informative & !is.finite(rho)] <- -Inf
  }

  return(value)
}


# Sums of x within the groups 1..n_group.
.group_sum <- function(x, group, n_group) {

  out <- numeric(n_group)
  if (length(x)) {
    sums <- rowsum(x, group, reorder = TRUE)
    out[as.integer(rownames(sums))] <- sums[, 1L]
  }
  return(out)
}


# Largest x within the groups 1..n_group; an empty group gives -Inf.
.group_max <- function(x, group, n_group) {

  out <- rep(-Inf, n_group)
  if (length(x)) {
    ordered <- order(group, -x)
    first   <- !duplicated(group[ordered])
    out[group[ordered][first]] <- x[ordered][first]
  }
  return(out)
}
