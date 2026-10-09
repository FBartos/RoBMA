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
    cells <- !is.finite(native)
    if (!any(cells)) {
      return(native)
    }
    # Only the cells the native kernel could not give are evaluated again.
    native[cells] <- .log_lik_cluster_norm_analytic_r(
      setup = setup,
      yi    = yi,
      sei   = sei,
      cells = cells
    )[cells]
    return(native)
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
    draws <- !is.finite(native)
    if (!any(draws)) {
      return(native)
    }
    # Only the draws whose sum the native kernel could not give are evaluated
    # again, over all their clusters.
    cells <- matrix(draws, nrow = length(native),
                    ncol = length(setup[["cluster"]]))
    log_lik <- .log_lik_cluster_norm_analytic_r(
      setup = setup,
      yi    = yi,
      sei   = sei,
      cells = cells
    )
    native[draws] <- .rowSums(log_lik[draws, , drop = FALSE],
                              sum(draws), ncol(log_lik))
    return(native)
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




# `cells` is NULL for every draw and cluster, or a logical S x G matrix naming
# the cells to evaluate; the others stay NA.
.log_lik_cluster_norm_analytic_r <- function(setup, yi, sei, cells = NULL) {

  cluster <- .cluster_indices_flatten(setup[["cluster"]])
  rows    <- cluster[["index"]]
  group   <- rep.int(seq_along(cluster[["size"]]), cluster[["size"]])
  weights <- setup[["weights"]][rows]
  S       <- setup[["S"]]
  log_lik <- matrix(NA_real_, nrow = S, ncol = length(cluster[["size"]]))
  draws   <- if (is.null(cells)) {
    seq_len(S)
  } else {
    which(.rowSums(cells, S, ncol(cells)) > 0L)
  }

  for (s in draws) {
    keep <- if (is.null(cells)) rep.int(TRUE, length(rows)) else cells[s, group]
    ids  <- unique(group[keep])
    log_lik[s, ids] <- .log_lik_cluster_norm_stable(
      residual = yi[rows[keep]] - setup[["mu"]][s, rows[keep]],
      sd       = .root_sum_squares(sei[rows[keep]],
                                   setup[["tau_within"]][s, rows[keep]]),
      loading  = setup[["tau_between"]][s, rows[keep]],
      weights  = weights[keep],
      group    = match(group[keep], ids)
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
# (a scale model gives each estimate its own). With z_i = residual_i / sd_i,
# b_i = loading_i / sd_i, A = 1 + sum(weight_i b_i^2) and the weighted latent
# location m = sum(weight_i b_i z_i) / A, the integral is
#
#   log L = -sum(weight (log(2 pi) + 2 log(sd))) / 2 - log(A) / 2
#           - (sum(weight_i (z_i - b_i m)^2) + m^2) / 2.
#
# The residual after the latent projection is formed directly, so no explained
# quadratic is subtracted from a total, and the form is evaluated in the regime
# where it neither overflows nor loses an exactly aligned residual. With
# t = max|loading| and s_min = min(sd) over the active entries of a cluster:
#
#  * t <= s_min: every |b_i| <= 1 and the form above is evaluated as written.
#  * t > s_min: the same form in raw units scaled by T, the power of two not
#    above t (so that a_i = loading_i / T is exact), with
#    p_i = weight_i (s_min / s_i)^2, D = (s_min / T)^2 + sum(p_i a_i^2),
#    c = sum(p_i a_i residual_i) / D (= T m), d_i = residual_i - a_i c,
#    log(A) = 2 (log(T) - log(s_min)) + log(D) and
#    quadratic = sum(weight_i (d_i / s_i)^2) + (c / T)^2. A cluster whose D or
#    whose numerator terms fall below 1e-290 carries a common log-scale factor
#    in D and c.
#  * t = 0: the independent Gaussian densities.
#
# No quotient by a loading is formed and a scale is squared only after it is
# divided by the one it is compared with. A cluster with a non-finite value or a
# non-positive sd gives NA; a result beyond the representable range is -Inf.
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

  z     <- r / s
  norm  <- .group_sum(-0.5 * w * (log(2 * pi) + 2 * log(s)), g, n_group)
  value <- norm - 0.5 * .group_sum(w * z * z, g, n_group)
  value[failed] <- NA_real_

  top    <- .group_max(abs(t), g, n_group)
  bottom <- -.group_max(-s, g, n_group)
  loaded <- top > 0
  scaled <- loaded & top > bottom
  plain  <- loaded & !scaled
  if (any(plain)) {
    overflow <- tabulate(g[!is.finite(z)], nbins = n_group) > 0L
    rows     <- plain[g]
    standard <- .log_lik_cluster_norm_standardized(
      z = z[rows], b = t[rows] / s[rows], w = w[rows], g = g[rows],
      n_group = n_group
    )
    latent <- norm - 0.5 * standard[["log_a"]] - 0.5 * standard[["quadratic"]]
    latent[overflow] <- -Inf
    value[plain] <- latent[plain]
  }
  if (any(scaled)) {
    rows   <- scaled[g]
    raw    <- .log_lik_cluster_norm_scaled(
      r = r[rows], s = s[rows], t = t[rows], w = w[rows], g = g[rows],
      top = ifelse(scaled, top, 1), bottom = ifelse(scaled, bottom, 1),
      groups = scaled, n_group = n_group
    )
    latent <- norm - 0.5 * raw[["log_a"]] - 0.5 * raw[["quadratic"]]
    value[scaled] <- latent[scaled]
  }

  return(value)
}


# The cluster log determinant term and quadratic of the form above for
# clusters whose loadings do not exceed their sds (|b| <= 1); the entries are
# the rows of those clusters only.
.log_lik_cluster_norm_standardized <- function(z, b, w, g, n_group) {

  information <- .group_sum(w * b * b, g, n_group)
  latent      <- .group_sum(w * b * z, g, n_group) / (1 + information)
  e           <- z - b * latent[g]

  return(list(
    log_a     = log1p(information),
    quadratic = .group_sum(w * e * e, g, n_group) + latent * latent
  ))
}


# The same for clusters with a loading above their smallest sd (`bottom`), in
# raw units scaled by the power of two not above the largest loading `top` (both
# 1 outside `groups`); `groups` marks the clusters of the rows. The returned
# log determinant term and quadratic are those of the form above.
.log_lik_cluster_norm_scaled <- function(r, s, t, w, g, top, bottom, groups,
                                         n_group) {

  scale <- 2^floor(log2(top))
  p     <- w * (bottom[g] / s)^2
  a     <- t / scale[g]
  denom <- (bottom / scale)^2 + .group_sum(p * a * a, g, n_group)
  numer <- .group_sum(p * a * r, g, n_group)
  log_k <- numeric(n_group)

  # A denominator or a numerator term that underflows is carried with a common
  # log-scale factor.
  term  <- abs(p * a * r)
  small <- tabulate(g[t != 0 & r != 0 & !(term >= 1e-290)], nbins = n_group) > 0L
  slow  <- groups & (small | !(is.finite(denom) & denom >= 1e-290))
  if (any(slow)) {
    rows  <- slow[g]
    g_s   <- g[rows]
    t_s   <- t[rows]
    r_s   <- r[rows]
    log_w <- log(w[rows]) + 2 * (log(bottom[g_s]) - log(s[rows]))
    log_a <- ifelse(t_s != 0, log(abs(t_s)) - log(scale[g_s]), -Inf)
    log_0 <- 2 * (log(bottom) - log(scale))
    log_t <- log_w + 2 * log_a
    log_c <- pmax(log_0, .group_max(log_t, g_s, n_group))
    log_x <- log_w + log_a - log_c[g_s]
    slow_denom <- exp(log_0 - log_c) +
      .group_sum(exp(log_t - log_c[g_s]), g_s, n_group)
    slow_term  <- ifelse(
      r_s == 0, 0, ifelse(
        log_x > 700, sign(t_s) * sign(r_s) * exp(log_x + log(abs(r_s))),
        sign(t_s) * exp(log_x) * r_s
      )
    )
    slow_numer <- .group_sum(slow_term, g_s, n_group)
    denom[slow] <- slow_denom[slow]
    numer[slow] <- slow_numer[slow]
    log_k[slow] <- log_c[slow]
  }

  center <- numer / denom
  latent <- center / scale
  # loading * latent equals a * center, and does not underflow with a
  shift  <- ifelse(t == 0, 0, ifelse(abs(a) >= 1e-290, a * center[g],
                                     t * latent[g]))
  d      <- (r - shift) / s

  return(list(
    log_a     = 2 * (log(scale) - log(bottom)) + log_k + log(denom),
    quadratic = .group_sum(w * d * d, g, n_group) + latent * latent
  ))
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
