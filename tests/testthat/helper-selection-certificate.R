# Fixtures, independent references and a plain-R shift-amplitude certificate for
# the product selected-normal block normalizer (factor and cluster kernels).
#
# A block is a list of groups of coincident rows. Row i of a group has the
# latent mean `mean`, residual SD `sd`, selection SE `sei` and the loading
# vector `loading` over the block's latent axes. The selection weight is the
# step function `omega` with lower cuts `z_lower`, and the normalizer is
#   Z = E[ prod_i P_i(m_i + a_i' eta) ],  eta ~ N(0, I),
#   P_i(t) = omega_B + sum_b (omega_b - omega_{b+1}) pnorm((t - z_b sei) / sd).

.test_cert_block <- function(groups, omega, z_lower = c(0, -Inf)) {

  rows <- unlist(lapply(groups, function(group) rep(1L, group[["count"]])))
  expand <- function(field) unlist(lapply(groups, function(group) {
    rep(group[[field]], group[["count"]])
  }))
  rank <- length(groups[[1L]][["loading"]])
  loading <- do.call(rbind, lapply(groups, function(group) {
    matrix(group[["loading"]], group[["count"]], rank, byrow = TRUE)
  }))

  list(
    groups  = groups,
    mean    = expand("mean"),
    sd      = expand("sd"),
    sei     = expand("sei"),
    loading = loading,
    omega   = omega,
    z_lower = z_lower,
    z_upper = c(Inf, utils::head(z_lower, -1L)),
    rows    = length(rows),
    rank    = rank
  )
}

# One group of K coincident rows on two parallel axes.
.test_cert_parallel_group <- function(count, loading, sd, mean, sei = 1) {

  list(count = count, loading = loading, sd = sd, mean = mean, sei = sei)
}

.test_cert_band_weights <- function(rho) c(1, 1 / rho, 1)
.test_cert_band_cuts    <- function() c(1, -1, -Inf)

# Native kernels. The normalizer does not depend on the observed values, so the
# observed responses are zero and every row is observed in the top bin.
.test_cert_native_factor <- function(block, tolerance = .005, initial = 256L,
                                     points = 2048L, scrambles = 8L, seed = 1L,
                                     rules = .selection_joint_factor_quadrature_rules()) {

  K   <- block[["rows"]]
  qmc <- BayesTools::selection_qmc_design(
    dimensions = 2L * block[["rank"]], points = points, scrambles = scrambles, seed = seed)

  .Call("RoBMA_selnorm_factor_step_loglik_batch", rep(0, K), matrix(block[["mean"]], 1L),
    matrix(block[["sd"]], 1L), matrix(as.double(block[["loading"]]), 1L), block[["sei"]],
    matrix(block[["omega"]], 1L), block[["z_lower"]], block[["z_upper"]], rep(1L, K), 1L,
    TRUE, 1L, rules[["nodes"]], rules[["log_weights"]], as.double(rules[["orders"]]),
    as.double(qmc), as.integer(initial), as.integer(points), as.integer(scrambles),
    tolerance, TRUE, 0L, PACKAGE = "RoBMA")
}

.test_cert_native_cluster <- function(block, tolerance = .005, initial = 256L,
                                      points = 2048L, scrambles = 8L, seed = 1L) {

  K     <- block[["rows"]]
  rules <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  qmc   <- BayesTools::selection_qmc_design(
    dimensions = 2L, points = points, scrambles = scrambles, seed = seed)

  .Call("RoBMA_selnorm_cluster_step_loglik_batch", rep(0, K), matrix(block[["mean"]], 1L),
    matrix(block[["sd"]], 1L), matrix(as.double(block[["loading"]]), 1L), block[["sei"]],
    matrix(block[["omega"]], 1L), block[["z_lower"]], block[["z_upper"]], rep(1L, K), 1L,
    TRUE, 1L, rules[["nodes"]], rules[["log_weights"]], as.double(rules[["orders"]]),
    as.double(qmc), as.integer(initial), as.integer(points), as.integer(scrambles),
    tolerance, TRUE, 0L, PACKAGE = "RoBMA")
}

# Whether the native kernel accepted a deterministic rule (no randomized
# fallback): the relative MCSE of a quadrature result is exactly zero.
.test_cert_deterministic <- function(result) identical(as.numeric(result[["relative_mcse"]]), 0)

.test_cert_relative_error <- function(log_normalizer, reference) {

  abs(expm1(log_normalizer - reference))
}

# ---------------------------------------------------------------------------- #
# Independent references (stats::integrate and composite Gauss-Legendre)
# ---------------------------------------------------------------------------- #

.test_cert_log_probability <- function(location, sd, sei, omega, z_lower) {

  bins  <- length(omega)
  value <- rep(omega[bins], length(location))
  for (bin in seq_len(bins - 1L)) {
    step <- omega[bin] - omega[bin + 1L]
    if (step == 0) next
    value <- value + step * stats::pnorm((location - z_lower[bin] * sei) / sd)
  }
  log(value)
}

.test_cert_logsumexp <- function(x) {

  top <- max(x)
  if (!is.finite(top)) return(top)
  top + log(sum(exp(x - top)))
}

.test_cert_gauss_legendre <- function(m) {

  i <- seq_len(m - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, m, m)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e[["values"]])
  list(x = e[["values"]][o], w = 2 * e[["vectors"]][1L, o]^2)
}

# Nodes and log weights of a composite Gauss-Legendre rule on [-limit, limit].
.test_cert_composite_rule <- function(limit, panel, m = 12L) {

  rule  <- .test_cert_gauss_legendre(m)
  edges <- seq(-limit, limit, by = panel)
  half  <- diff(edges) / 2
  mid   <- utils::head(edges, -1L) + half
  list(
    x     = as.vector(outer(rule[["x"]], half) + rep(mid, each = m)),
    log_w = log(as.vector(outer(rule[["w"]], half)))
  )
}

# Whether every connected set of latent axes of the block carries parallel rows
# (each row a multiple of one direction), so a one-dimensional reduction holds.
.test_cert_parallel_rows <- function(block) {

  loading <- block[["loading"]]
  active  <- loading != 0
  label   <- seq_len(ncol(loading))
  repeat {
    changed <- FALSE
    for (row in seq_len(nrow(active))) {
      columns <- which(active[row, ])
      if (!length(columns)) next
      target <- min(label[columns])
      if (any(label[columns] != target)) {
        label[label %in% label[columns]] <- target
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  for (set in unique(label[colSums(active) > 0])) {
    columns <- which(label == set)
    rows    <- which(rowSums(active[, columns, drop = FALSE]) > 0)
    a       <- loading[rows, columns, drop = FALSE]
    u       <- a[1L, ] / sqrt(sum(a[1L, ]^2))
    if (max(abs(a - outer(as.vector(a %*% u), u))) > 1e-12 * max(abs(a))) return(FALSE)
  }
  TRUE
}

# The independent reference of a block: the one-dimensional reduction of
# parallel rows, else the nested integral of two groups of coincident rows.
.test_cert_reference <- function(block) {

  if (.test_cert_parallel_rows(block)) .test_cert_reference_parallel(block)
  else .test_cert_reference_two_groups(block)
}

# Blocks whose rows are parallel within connected sets of latent axes: the
# integral factors over the sets and each factor is one-dimensional along the
# common direction u with scaled loadings c_i = a_i' u. The panel width stays
# below the smallest transition scale s_i / (8 |c_i|).
.test_cert_reference_parallel <- function(block, limit = 12, m = 12L) {

  loading <- block[["loading"]]
  active  <- loading != 0
  label   <- seq_len(ncol(loading))
  # latent axes joined by a row that loads on both belong to one set
  repeat {
    changed <- FALSE
    for (row in seq_len(nrow(active))) {
      columns <- which(active[row, ])
      if (!length(columns)) next
      target <- min(label[columns])
      if (any(label[columns] != target)) {
        label[label %in% label[columns]] <- target
        changed <- TRUE
      }
    }
    if (!changed) break
  }

  total <- 0
  for (set in unique(label[colSums(active) > 0])) {
    columns <- which(label == set)
    rows    <- which(rowSums(active[, columns, drop = FALSE]) > 0)
    a       <- loading[rows, columns, drop = FALSE]
    u       <- a[1L, ] / sqrt(sum(a[1L, ]^2))
    scale   <- as.vector(a %*% u)
    # parallel rows: every row is its scale times the common direction
    stopifnot(max(abs(a - outer(scale, u))) <= 1e-12 * max(abs(a)))
    panel <- min(.02, min(block[["sd"]][rows] / (8 * abs(scale))))
    rule  <- .test_cert_composite_rule(limit, panel, m)
    value <- stats::dnorm(rule[["x"]], log = TRUE) + rule[["log_w"]]
    # coincident rows contribute their log probability once per copy
    keys   <- paste(block[["mean"]][rows], block[["sd"]][rows], block[["sei"]][rows], scale)
    unique <- unique(keys)
    copies <- tabulate(match(keys, unique), length(unique))
    for (i in seq_along(unique)) {
      position <- match(unique[i], keys)
      row      <- rows[position]
      value    <- value + copies[i] * .test_cert_log_probability(
        block[["mean"]][row] + scale[position] * rule[["x"]],
        block[["sd"]][row], block[["sei"]][row], block[["omega"]], block[["z_lower"]])
    }
    total <- total + .test_cert_logsumexp(value)
  }
  total
}

# Two groups of coincident rows on two latent axes (any angle): nested
# composite Gauss-Legendre in coordinates aligned with the first group.
.test_cert_reference_two_groups <- function(block, limit = 10, panel = .1, m = 12L) {

  first  <- block[["groups"]][[1L]]
  second <- block[["groups"]][[2L]]
  stopifnot(length(block[["groups"]]) == 2L, block[["rank"]] == 2L)
  u <- first[["loading"]] / sqrt(sum(first[["loading"]]^2))
  v <- c(-u[2L], u[1L])
  scale_first <- sqrt(sum(first[["loading"]]^2))
  along  <- sum(second[["loading"]] * u)
  across <- sum(second[["loading"]] * v)
  rule   <- .test_cert_composite_rule(limit, panel, m)
  log_wy <- stats::dnorm(rule[["x"]], log = TRUE) + rule[["log_w"]]

  outer_value <- numeric(length(rule[["x"]]))
  chunk <- 200L
  for (start in seq(1L, length(rule[["x"]]), by = chunk)) {
    index    <- start:min(start + chunk - 1L, length(rule[["x"]]))
    location <- outer(rule[["x"]][index] * along, rule[["x"]] * across, `+`) + second[["mean"]]
    inner    <- second[["count"]] * matrix(.test_cert_log_probability(as.vector(location),
      second[["sd"]], second[["sei"]], block[["omega"]], block[["z_lower"]]), nrow(location))
    inner    <- inner + rep(log_wy, each = nrow(location))
    top      <- apply(inner, 1L, max)
    outer_value[index] <- top + log(rowSums(exp(inner - top)))
  }
  first_value <- first[["count"]] * .test_cert_log_probability(
    first[["mean"]] + scale_first * rule[["x"]], first[["sd"]], first[["sei"]],
    block[["omega"]], block[["z_lower"]])
  .test_cert_logsumexp(log_wy + first_value + outer_value)
}

# A shared first axis (root) and independent leaf axes, each leaf owning one
# group: given the root, the leaves are independent one-dimensional integrals.
.test_cert_reference_forest <- function(block) {

  groups <- block[["groups"]]
  leaf <- function(group, position, root) {
    vapply(root, function(x) {
      stats::integrate(function(z) {
        stats::dnorm(z) * exp(group[["count"]] * .test_cert_log_probability(
          group[["mean"]] + group[["loading"]][1L] * x + group[["loading"]][position] * z,
          group[["sd"]], group[["sei"]], block[["omega"]], block[["z_lower"]]))
      }, -12, 12, rel.tol = 1e-11, abs.tol = 1e-13, subdivisions = 1000L)[["value"]]
    }, numeric(1L))
  }
  value <- stats::integrate(function(x) {
    stats::dnorm(x) * Reduce(`*`, lapply(seq_along(groups), function(g) leaf(groups[[g]], g + 1L, x)))
  }, -10, 10, rel.tol = 1e-10, abs.tol = 1e-12, subdivisions = 1000L)
  log(value[["value"]])
}

# ---------------------------------------------------------------------------- #
# Plain-R shift-amplitude certificate (tensor Gauss-Hermite rung)
# ---------------------------------------------------------------------------- #

# Log of the tensor Gauss-Hermite rule of this order, evaluated at the nodes
# moved by `shift` along the axes with the exactly changed weights
# w_j exp(-shift x_j - shift^2 / 2). With a zero shift this is the plain rule.
.test_cert_tensor_log_integral <- function(block, order, shift = rep(0, block[["rank"]]),
                                           moments = FALSE) {

  rule <- .gauss_hermite_nodes(order)
  rank <- block[["rank"]]
  grid <- as.matrix(expand.grid(rep(list(seq_len(order)), rank)))
  x    <- matrix(rule[["nodes"]][grid], ncol = rank)
  lw   <- rowSums(matrix(rule[["log_weights"]][grid], ncol = rank))
  lw   <- lw - as.vector(x %*% shift) - sum(shift^2) / 2
  eta  <- sweep(x, 2L, shift, `+`)
  location <- eta %*% t(block[["loading"]]) + rep(block[["mean"]], each = nrow(eta))
  value <- lw
  for (i in seq_len(block[["rows"]])) {
    value <- value + .test_cert_log_probability(location[, i], block[["sd"]][i],
      block[["sei"]][i], block[["omega"]], block[["z_lower"]])
  }
  log_integral <- .test_cert_logsumexp(value)
  if (!moments) return(log_integral)
  list(log_integral = log_integral, mean = colSums(x * exp(value - log_integral)))
}

# The gap of the order-n nodes around a mean: R's findInterval clamped to the
# interior intervals; a mean on a node uses the gap to its right.
.test_cert_node_gap <- function(order, mean) {

  nodes <- .gauss_hermite_nodes(order)[["nodes"]]
  k     <- min(max(findInterval(mean, nodes), 1L), order - 1L)
  nodes[k + 1L] - nodes[k]
}

# Per-axis amplitudes from the half and quarter shifts at the rung's own tilted
# means; every axis is its own group.
.test_cert_amplitudes <- function(block, order) {

  base <- .test_cert_tensor_log_integral(block, order, moments = TRUE)
  q0   <- base[["log_integral"]]
  gap  <- vapply(base[["mean"]], function(m) .test_cert_node_gap(order, m), numeric(1L))
  shifted <- function(fraction) vapply(seq_len(block[["rank"]]), function(axis) {
    shift <- rep(0, block[["rank"]])
    shift[axis] <- fraction * gap[axis]
    .test_cert_tensor_log_integral(block, order, shift)
  }, numeric(1L))
  half    <- shifted(.5)
  quarter <- shifted(.25)
  e2 <- expm1(half - q0)
  e4 <- expm1(quarter - q0)
  list(
    log_integral = q0, mean = base[["mean"]], gap = gap, half = half, quarter = quarter,
    amplitude = sqrt((e2 / 2)^2 + (e2 / 2 - e4)^2)
  )
}

# ---------------------------------------------------------------------------- #
# Native certificate probe
# ---------------------------------------------------------------------------- #

# One rung of the product rule with its shift-amplitude certificate, on the
# routines the kernels use. `order` is the rung's order in the factor ladder, or
# in the rank-one cluster ladder when `cluster = TRUE`.
.test_cert_probe <- function(block, order, cluster = FALSE, tolerance = .005,
                             shift_axis = integer(), shift_delta = numeric(),
                             telescope = TRUE, precomputed = FALSE, refine = TRUE) {

  rules <- if (cluster) {
    .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  } else {
    .selection_joint_factor_quadrature_rules()
  }
  .Call("RoBMA_selnorm_certificate_probe", cluster, block[["mean"]], block[["sd"]],
    as.double(block[["loading"]]), as.integer(block[["rank"]]), block[["sei"]],
    block[["omega"]], block[["z_lower"]], block[["z_upper"]], 1L, telescope,
    rules[["nodes"]], rules[["log_weights"]], as.double(rules[["orders"]]),
    match(order, rules[["orders"]]), tolerance, as.integer(shift_axis),
    as.double(shift_delta), precomputed, refine, PACKAGE = "RoBMA")
}

# The path codes of the probe.
.test_cert_paths <- c(declined = 0L, nested = 1L, tensor_direct = 2L, tensor_log = 3L, cluster = 4L)

# ---------------------------------------------------------------------------- #
# Plain-R emulation of the acceptance ladder (independent of the native code)
# ---------------------------------------------------------------------------- #

# The exterior bound of the tail gate on the relative scale: the mass of the
# product-rule integrand outside the node box, at most omega_max per row.
.test_cert_tail_bound <- function(block, order, log_integral) {

  nodes <- .gauss_hermite_nodes(order)[["nodes"]]
  log_tail <- log(stats::pnorm(min(nodes)) + stats::pnorm(max(nodes), lower.tail = FALSE)) +
    log(block[["rank"]]) + log(max(block[["omega"]])) * block[["rows"]]
  exp(log_tail - log_integral)
}

# The deterministic ladder the kernels run for a product block: a rung needs a
# finite value, the successive changes of the rung and (factor kernel) of the
# previous rung within the tolerance, the tail gate, and the certificate with
# its combined allowance. Returns the accepted rung, or NULL when the block
# reaches the randomized fallback.
.test_cert_emulate_ladder <- function(block, tolerance = .005,
                                      orders = SELNORM_FACTOR_QUADRATURE_ORDERS,
                                      changes = 2L) {

  transitions <- any(is.finite(block[["z_lower"]]) &
    c(block[["omega"]][-length(block[["omega"]])] != block[["omega"]][-1L], FALSE))
  previous <- -Inf
  previous_change <- Inf
  for (k in seq_along(orders)) {
    value  <- .test_cert_tensor_log_integral(block, orders[k])
    difference <- previous - value
    change <- if (k == 1L || !is.finite(difference) || abs(difference) > 50) Inf else abs(expm1(difference))
    cheap <- is.finite(value) && change <= tolerance &&
      (changes == 1L || previous_change <= tolerance)
    if (cheap) {
      tail <- .test_cert_tail_bound(block, orders[k], value)
      if (tail <= tolerance) {
        amplitude <- if (transitions) sum(.test_cert_amplitudes(block, orders[k])[["amplitude"]]) else 0
        if (amplitude <= tolerance / 2 && tail + 1.3 * amplitude <= tolerance / (1 + tolerance)) {
          return(list(order = orders[k], value = value, change = change, tail = tail, amplitude = amplitude))
        }
      }
    }
    previous <- value
    previous_change <- change
  }
  NULL
}

# ---------------------------------------------------------------------------- #
# Native refinement tier
# ---------------------------------------------------------------------------- #

# The decision of the tier recomputed from the probe: the M = 4 average of the
# finest and previous affordable rung (orders 95 and 63 at rank two and for
# forests of depth one), the shift check, the tail gate and the order check.
.test_cert_tier <- function(block, tolerance = .005, finest = 95L, previous = 63L) {

  fine <- .test_cert_probe(block, finest, tolerance = tolerance)
  prior <- .test_cert_probe(block, previous, tolerance = tolerance)
  shift_change <- abs(expm1(fine[["refined_4"]] - fine[["refined_2"]]))
  order_change <- abs(expm1(fine[["refined_4"]] - prior[["refined_4"]]))
  tail <- .test_cert_tail_bound(block, finest, fine[["refined_4"]])
  list(
    value = fine[["refined_4"]], shift_change = shift_change, order_change = order_change,
    tail = tail,
    accepted = shift_change <= tolerance / 2 && tail <= tolerance && order_change <= tolerance / 2
  )
}

# The accepted order of the certificate ladder recomputed from the probe's rung
# values: finite, the two successive changes and the tail gate within the
# tolerance, and the certificate; NA when every affordable rung is rejected.
.test_cert_probe_ladder <- function(block, tolerance = .005) {

  previous <- -Inf
  previous_change <- Inf
  for (order in SELNORM_FACTOR_QUADRATURE_ORDERS) {
    probe <- .test_cert_probe(block, order, tolerance = tolerance, refine = FALSE)
    if (probe[["path"]] == .test_cert_paths[["declined"]]) break
    value <- probe[["log_normalizer"]]
    difference <- previous - value
    change <- if (!is.finite(difference) || abs(difference) > 50) Inf else abs(expm1(difference))
    if (is.finite(value) && change <= tolerance && previous_change <= tolerance &&
        probe[["tail_bound"]] <= tolerance && probe[["accepted"]]) {
      return(order)
    }
    previous <- value
    previous_change <- change
  }
  NA_integer_
}
