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
# `design` replaces the default QMC design (scrambles x points x dimensions).
.test_cert_native_factor <- function(block, tolerance = .005, initial = 256L,
                                     points = 2048L, scrambles = 8L, seed = 1L,
                                     rules = .selection_joint_factor_quadrature_rules(),
                                     design = NULL) {

  K   <- block[["rows"]]
  qmc <- if (is.null(design)) BayesTools::selection_qmc_design(
    dimensions = 2L * block[["rank"]], points = points, scrambles = scrambles, seed = seed) else design

  .Call("RoBMA_selnorm_factor_step_loglik_batch", rep(0, K), matrix(block[["mean"]], 1L),
    matrix(block[["sd"]], 1L), matrix(as.double(block[["loading"]]), 1L), block[["sei"]],
    matrix(block[["omega"]], 1L), block[["z_lower"]], block[["z_upper"]], rep(1L, K), 1L,
    TRUE, 1L, rules[["nodes"]], rules[["log_weights"]], as.double(rules[["orders"]]),
    as.double(qmc), as.integer(initial), as.integer(points), as.integer(scrambles),
    tolerance, TRUE, 0L, PACKAGE = "RoBMA")
}

.test_cert_native_cluster <- function(block, tolerance = .005, initial = 256L,
                                      points = 2048L, scrambles = 8L, seed = 1L,
                                      rules = .selection_joint_cluster_quadrature_rules(
                                        SELNORM_CLUSTER_QUADRATURE_ORDERS)) {

  K     <- block[["rows"]]
  qmc   <- BayesTools::selection_qmc_design(
    dimensions = 2L, points = points, scrambles = scrambles, seed = seed)

  .Call("RoBMA_selnorm_cluster_step_loglik_batch", rep(0, K), matrix(block[["mean"]], 1L),
    matrix(block[["sd"]], 1L), matrix(as.double(block[["loading"]]), 1L), block[["sei"]],
    matrix(block[["omega"]], 1L), block[["z_lower"]], block[["z_upper"]], rep(1L, K), 1L,
    TRUE, 1L, rules[["nodes"]], rules[["log_weights"]], as.double(rules[["orders"]]),
    as.double(qmc), as.integer(initial), as.integer(points), as.integer(scrambles),
    tolerance, TRUE, 0L, PACKAGE = "RoBMA")
}

# One block of the stored best-rule reference through the kernel of its rank
# (the cluster kernel for rank one), with the controls of that reference: the
# fallback budgets double from 64 to 512 points per scramble (`initial`, `points`).
.test_cert_native_best <- function(block,
                                   factor_rules = .selection_joint_factor_quadrature_rules(),
                                   cluster_rules = .selection_joint_cluster_quadrature_rules(
                                     SELNORM_CLUSTER_QUADRATURE_ORDERS),
                                   initial = 64L, points = 512L) {

  if (block[["rank"]] == 1L) {
    qmc <- BayesTools::selection_qmc_design(2L, points, 8L, seed = 173L)
    .Call("RoBMA_selnorm_cluster_step_loglik_batch", block[["yi"]], matrix(block[["mean"]], 1L),
      matrix(block[["sd"]], 1L), matrix(block[["loading"]][, 1L], 1L), block[["sei"]],
      matrix(block[["omega"]], 1L), block[["z_lower"]], block[["z_upper"]], block[["obs_bin"]],
      1L, TRUE, 1L, cluster_rules[["nodes"]], cluster_rules[["log_weights"]],
      as.double(cluster_rules[["orders"]]), as.double(qmc), as.integer(initial), as.integer(points),
      8L, .005, TRUE, as.integer(block[["rule"]]), PACKAGE = "RoBMA")
  } else {
    qmc <- BayesTools::selection_qmc_design(2L * block[["rank"]], points, 8L, seed = 173L)
    .Call("RoBMA_selnorm_factor_step_loglik_batch", block[["yi"]], matrix(block[["mean"]], 1L),
      matrix(block[["sd"]], 1L), matrix(as.double(block[["loading"]]), 1L), block[["sei"]],
      matrix(block[["omega"]], 1L), block[["z_lower"]], block[["z_upper"]], block[["obs_bin"]],
      1L, TRUE, 1L, factor_rules[["nodes"]], factor_rules[["log_weights"]],
      as.double(factor_rules[["orders"]]), as.double(qmc), as.integer(initial), as.integer(points),
      8L, .005, TRUE, as.integer(block[["rule"]]), PACKAGE = "RoBMA")
  }
}

# A stored output of the build before the stopping multiplier (z = 1, i.e. one standard
# error) is reproduced by a run of this build unless the run went on to a later stage of
# its nested design. Then the stored output is the stage at which the one-standard-error
# rule stops and the run is the stage at which the scaled rule stops, each stage from a run
# with that single budget (`stage_run(points)`; the design of a stage is the leading part
# of the next one). A change of proposals or of the design shows as an unreproduced stage.
.test_cert_expect_stored_fallback <- function(actual, expected, stage_run, stages, tolerance = .005,
                                              fields = c("log_normalizer", "relative_mcse", "relative_change"),
                                              info = "") {

  matches <- function(x, y) all(vapply(fields, function(field) {
    isTRUE(all.equal(as.numeric(x[[field]]), as.numeric(y[[field]]), tolerance = 1e-12))
  }, logical(1L)))
  if (matches(actual, expected)) return(invisible("unchanged"))

  results <- lapply(stages, stage_run)
  stop_stage <- function(multiplier) {
    criterion <- vapply(results, function(x) max(multiplier * x[["relative_mcse"]], x[["relative_change"]]),
                        numeric(1L))
    within <- which(criterion <= tolerance)
    if (length(within)) within[[1L]] else length(results)
  }
  old <- stop_stage(1)
  new <- stop_stage(.selection_qmc_z())
  testthat::expect_lt(old, new, label = paste(info, "stage of the one-standard-error rule"))
  for (field in fields) {
    testthat::expect_equal(as.numeric(results[[old]][[field]]), as.numeric(expected[[field]]),
      tolerance = 1e-12, info = paste(info, "stored output, stage", stages[[old]], field))
    testthat::expect_equal(as.numeric(results[[new]][[field]]), as.numeric(actual[[field]]),
      tolerance = 1e-12, info = paste(info, "run, stage", stages[[new]], field))
  }
  invisible("later stage")
}

# Whether the native kernel accepted a deterministic rule (no randomized
# fallback): the kernel's own flag, the attribute "deterministic" of its result.
# A randomized result can have a relative MCSE of exactly zero (scrambles that
# coincide), so the MCSE does not identify the route; results of builds that
# predate the flag carry none and are compared on their relative MCSE directly.
.test_cert_deterministic <- function(result) {

  flag <- attr(result, "deterministic")
  if (is.null(flag)) stop("The kernel result carries no 'deterministic' flag.", call. = FALSE)
  identical(as.logical(flag), TRUE)
}

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

# Half-width of the latent integration range of a reference. The integrand is the
# standard normal weight times a product of selection probabilities that is
# constant beyond the transition of each row, at the distance
# |cut * sei - mean| / ||loading|| of its hyperplane from the origin, so the mass
# lies within the farthest transition plus a margin of 8 (blocks of many
# coincident rows are carried by the plateau beyond a far transition). A
# reference never uses less than `floor`, the half-width it used for every block
# before the range followed the block.
.test_cert_reference_limit <- function(block, floor, margin = 8) {

  size   <- sqrt(rowSums(block[["loading"]]^2))
  cuts   <- block[["z_lower"]][is.finite(block[["z_lower"]])]
  active <- which(size > 0)
  if (!length(cuts) || !length(active)) return(floor)
  distance <- abs(outer(block[["sei"]][active], cuts) - block[["mean"]][active]) / size[active]
  max(floor, max(distance) + margin)
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
# The integrand at the ends of the range is checked to be negligible against its
# maximum, so a range that cuts off the mass fails instead of returning a value.
.test_cert_reference_parallel <- function(block, limit = .test_cert_reference_limit(block, 12),
                                          m = 12L) {

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
    if (max(value[c(1L, length(value))]) > max(value) - 30) {
      stop("The reference range [-", limit, ", ", limit, "] cuts off the mass of the block.", call. = FALSE)
    }
    total <- total + .test_cert_logsumexp(value)
  }
  total
}

# Two groups of coincident rows on two latent axes (any angle): nested
# composite Gauss-Legendre in coordinates aligned with the first group.
.test_cert_reference_two_groups <- function(block, limit = .test_cert_reference_limit(block, 10),
                                             panel = .1, m = 12L) {

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
  leaf_limit <- .test_cert_reference_limit(block, 12)
  root_limit <- .test_cert_reference_limit(block, 10)
  leaf <- function(group, position, root) {
    vapply(root, function(x) {
      stats::integrate(function(z) {
        stats::dnorm(z) * exp(group[["count"]] * .test_cert_log_probability(
          group[["mean"]] + group[["loading"]][1L] * x + group[["loading"]][position] * z,
          group[["sd"]], group[["sei"]], block[["omega"]], block[["z_lower"]]))
      }, -leaf_limit, leaf_limit, rel.tol = 1e-11, abs.tol = 1e-13, subdivisions = 1000L)[["value"]]
    }, numeric(1L))
  }
  value <- stats::integrate(function(x) {
    stats::dnorm(x) * Reduce(`*`, lapply(seq_along(groups), function(g) leaf(groups[[g]], g + 1L, x)))
  }, -root_limit, root_limit, rel.tol = 1e-10, abs.tol = 1e-12, subdivisions = 1000L)
  log(value[["value"]])
}

# Any block of groups of coincident rows on two latent axes: nested composite
# Gauss-Legendre, outer over the first axis and inner over the second. The
# panel width stays below the smallest transition width of the rows.
.test_cert_reference_rank_two <- function(block, limit = .test_cert_reference_limit(block, 8),
                                          panel = .1, m = 10L) {

  stopifnot(block[["rank"]] == 2L)
  rule <- .test_cert_composite_rule(limit, panel, m)
  log_w <- stats::dnorm(rule[["x"]], log = TRUE) + rule[["log_w"]]
  outer_value <- vapply(rule[["x"]], function(x) {
    value <- log_w
    for (group in block[["groups"]]) {
      location <- group[["mean"]] + group[["loading"]][1L] * x + group[["loading"]][2L] * rule[["x"]]
      value <- value + group[["count"]] * .test_cert_log_probability(location, group[["sd"]],
        group[["sei"]], block[["omega"]], block[["z_lower"]])
    }
    .test_cert_logsumexp(value)
  }, numeric(1L))
  .test_cert_logsumexp(log_w + outer_value)
}

# Nested groups on the prefix sums of the latent axes: group g loads equally
# (c_g) on axes 1..g, so its rows sit at mean_g + c_g S_g with the Gaussian random
# walk S_g = z_1 + ... + z_g. The normalizer is a chain of one-dimensional
# convolutions of the partial integrals on a composite Gauss-Legendre grid.
.test_cert_reference_chain <- function(block, panel = .1, m = 8L) {

  groups <- block[["groups"]]
  count  <- length(groups)
  stopifnot(block[["rank"]] == count)
  scale <- vapply(groups, function(group) group[["loading"]][1L], numeric(1L))
  for (g in seq_len(count)) {
    stopifnot(identical(groups[[g]][["loading"]], c(rep(scale[g], g), rep(0, count - g))))
  }
  rule   <- .test_cert_composite_rule(8 * sqrt(count), panel, m)
  weight <- exp(rule[["log_w"]])
  kernel <- stats::dnorm(outer(rule[["x"]], rule[["x"]], `-`))
  row_value <- function(g) {
    group <- groups[[g]]
    exp(group[["count"]] * .test_cert_log_probability(group[["mean"]] + scale[g] * rule[["x"]],
      group[["sd"]], group[["sei"]], block[["omega"]], block[["z_lower"]]))
  }
  partial <- stats::dnorm(rule[["x"]]) * row_value(1L)
  log_scale <- 0
  for (g in seq_len(count)[-1L]) {
    partial <- as.vector(kernel %*% (weight * partial)) * row_value(g)
    top <- max(partial)
    partial <- partial / top
    log_scale <- log_scale + log(top)
  }
  log_scale + log(sum(weight * partial))
}

# Groups on adjacent axes of a cycle: group g loads a_g on axis g and b_g on
# axis g + 1 (the last on axis 1), so Z = trace(prod_g diag(phi w) W_g) with
# W_g[x, y] = P_g(mean_g + a_g x + b_g y)^count_g on a composite grid.
.test_cert_reference_cycle <- function(block, limit = 8, panel = .2, m = 8L) {

  groups <- block[["groups"]]
  count  <- length(groups)
  stopifnot(block[["rank"]] == count)
  rule   <- .test_cert_composite_rule(limit, panel, m)
  weight <- exp(stats::dnorm(rule[["x"]], log = TRUE) + rule[["log_w"]])
  product <- diag(length(weight))
  for (g in seq_len(count)) {
    group <- groups[[g]]
    next_axis <- g %% count + 1L
    expected <- numeric(count)
    expected[g] <- group[["loading"]][g]
    expected[next_axis] <- group[["loading"]][next_axis]
    stopifnot(identical(group[["loading"]], expected))
    location <- group[["mean"]] + group[["loading"]][g] * rule[["x"]] +
      rep(group[["loading"]][next_axis] * rule[["x"]], each = length(weight))
    transition <- matrix(exp(group[["count"]] * .test_cert_log_probability(location,
      group[["sd"]], group[["sei"]], block[["omega"]], block[["z_lower"]])),
      length(weight))
    product <- product %*% (weight * transition)
  }
  log(sum(diag(product)))
}

# Two groups on two independent loadings whose weights vanish outside a band
# (zero outer bins). With m_g = a_g' z the normalizer is the integral over the two
# row locations of phi_2(A^-1 m) / |det A| times the two band probabilities, taken
# on the box where either band probability is not negligible (eight residual SDs
# beyond the outer cuts) with a composite Gauss-Legendre rule in each location.
.test_cert_reference_band_pair <- function(block, m = 8L) {

  groups <- block[["groups"]]
  stopifnot(length(groups) == 2L, block[["rank"]] == 2L,
    block[["omega"]][1L] == 0, block[["omega"]][length(block[["omega"]])] == 0)
  loading <- rbind(groups[[1L]][["loading"]], groups[[2L]][["loading"]])
  inverse <- solve(loading)
  cuts    <- block[["z_lower"]][is.finite(block[["z_lower"]])]
  axis <- lapply(groups, function(group) {
    low  <- min(cuts) * group[["sei"]] - 8 * group[["sd"]] - group[["mean"]]
    high <- max(cuts) * group[["sei"]] + 8 * group[["sd"]] - group[["mean"]]
    rule <- .test_cert_composite_rule((high - low) / 2, group[["sd"]] / 2, m)
    list(x = rule[["x"]] + (high + low) / 2, log_w = rule[["log_w"]], group = group)
  })
  log_band <- function(item) item[["group"]][["count"]] * .test_cert_log_probability(
    item[["group"]][["mean"]] + item[["x"]], item[["group"]][["sd"]],
    item[["group"]][["sei"]], block[["omega"]], block[["z_lower"]])
  first  <- axis[[1L]]
  second <- axis[[2L]]
  log_first  <- first[["log_w"]] + log_band(first)
  log_second <- second[["log_w"]] + log_band(second)
  grid <- as.matrix(expand.grid(first = first[["x"]], second = second[["x"]]))
  latent <- grid %*% t(inverse)
  log_density <- stats::dnorm(latent[, 1L], log = TRUE) + stats::dnorm(latent[, 2L], log = TRUE) -
    log(abs(det(loading)))
  .test_cert_logsumexp(log_density + rep(log_first, times = length(log_second)) +
    rep(log_second, each = length(log_first)))
}

# ---------------------------------------------------------------------------- #
# Plain-R shift-amplitude certificate (tensor Gauss-Hermite rung)
# ---------------------------------------------------------------------------- #

# Nodes and log terms of the tensor Gauss-Hermite rule of this order, evaluated
# at the nodes moved by `shift` along the axes with the exactly changed weights
# w_j exp(-shift x_j - shift^2 / 2): the grid rows `x` and the log of weight times
# integrand at each of them.
.test_cert_tensor_terms <- function(block, order, shift = rep(0, block[["rank"]])) {

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
  list(x = x, value = value)
}

# Log of the tensor Gauss-Hermite rule of this order (see above). With a zero
# shift this is the plain rule; `moments = TRUE` adds the tilted mean per axis.
.test_cert_tensor_log_integral <- function(block, order, shift = rep(0, block[["rank"]]),
                                           moments = FALSE) {

  terms <- .test_cert_tensor_terms(block, order, shift)
  log_integral <- .test_cert_logsumexp(terms[["value"]])
  if (!moments) return(log_integral)
  list(log_integral = log_integral,
       mean = colSums(terms[["x"]] * exp(terms[["value"]] - log_integral)))
}

# Mean and covariance of the latent coordinates under the plain tensor rule (the
# order of the loading columns).
.test_cert_tensor_moments <- function(block, order) {

  terms <- .test_cert_tensor_terms(block, order)
  p     <- exp(terms[["value"]] - .test_cert_logsumexp(terms[["value"]]))
  mean  <- colSums(terms[["x"]] * p)
  centered <- sweep(terms[["x"]], 2L, mean)
  list(mean = mean, covariance = crossprod(centered * sqrt(p)))
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

# The factor of the summed amplitude in the combined allowance, and the scales
# of the sharpness gate (the step widths the transition rule compares with a
# node gap are multiplied by them): the ladder rungs need an edge width of at
# least 0.15 of the node gap, the refinement tier 0.125 for a single row per
# axis (see .test_cert_tier_width_scale for the narrowing by coincident rows).
# Mirror the constants of src/selnorm/selnorm-mv.cc.
.test_cert_allowance_factor   <- 2
.test_cert_width_scale_ladder <- 1 / .15
.test_cert_width_scale_tier   <- 8

# Per-axis tier scales, with every row counted on each axis it loads on.
# The scalar diagnostic is their minimum, retained for the native probe.
.test_cert_tier_width_scales <- function(block) {

  omega <- block[["omega"]]
  transitions <- any(is.finite(block[["z_lower"]]) & c(omega[-length(omega)] != omega[-1L], FALSE))
  if (!transitions) return(rep(.test_cert_width_scale_tier, block[["rank"]]))
  rows <- pmax(1L, colSums(block[["loading"]] != 0))
  .test_cert_width_scale_tier / sqrt(pmax(1, 2 * log(rows)))
}

.test_cert_tier_width_scale <- function(block) min(.test_cert_tier_width_scales(block))

# Tier mirror: the geometry is unchanged, including the largest central node
# gap for axes coupled through any row. Only the width scale varies by axis.
.test_cert_tier_resolves <- function(block, order) {

  omega <- block[["omega"]]
  cuts <- which(is.finite(block[["z_lower"]]) & c(utils::head(omega, -1L) != omega[-1L], FALSE))
  if (!length(cuts)) return(TRUE)
  # as the native rule: fewer than two nodes resolve nothing
  if (order < 2L) return(FALSE)
  nodes <- .gauss_hermite_nodes(order)[["nodes"]]
  loaded <- block[["loading"]] != 0
  scales <- .test_cert_tier_width_scales(block)
  central_gaps <- diff(nodes)[utils::head(nodes, -1L) <= 3 & utils::tail(nodes, -1L) >= -3]
  # no central gap (all nodes beyond +-3) is a gap of zero, which the native
  # rule rejects for any coupled row
  central_gap <- if (length(central_gaps)) max(central_gaps) else 0
  for (axis in seq_len(block[["rank"]])) {
    rows <- which(loaded[, axis])
    coupled <- any(rowSums(loaded[rows, , drop = FALSE]) >= 2L)
    for (row in rows) {
      coefficient <- block[["loading"]][row, axis]
      if (coupled) {
        if (!(central_gap > 0) ||
            central_gap * abs(coefficient) > scales[axis] * block[["sd"]][row]) return(FALSE)
        next
      }
      for (cut in cuts) {
        center <- (block[["z_lower"]][cut] * block[["sei"]][row] - block[["mean"]][row]) / coefficient
        if (center <= nodes[1L] || center >= nodes[order]) next
        upper <- which(nodes >= center)[1L]
        gap <- nodes[upper] - nodes[upper - 1L]
        if (nodes[upper] == center && upper < order) gap <- max(gap, nodes[upper + 1L] - nodes[upper])
        if (gap * abs(coefficient) > scales[axis] * block[["sd"]][row]) return(FALSE)
      }
    }
  }
  TRUE
}

# The sharpness gate at a rung: the transition rule (its R mirror) with every
# step width multiplied by `width_scale`, on the nodes of the rung's order.
.test_cert_resolves <- function(block, order, width_scale = .test_cert_width_scale_ladder) {

  .selection_transition_rule_resolves(
    .gauss_hermite_nodes(order)[["nodes"]], block[["mean"]], block[["sd"]] * width_scale,
    block[["loading"]], block[["sei"]],
    list(omega = matrix(block[["omega"]], 1L), z_lower = block[["z_lower"]], sign = 1))
}

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
# previous rung within the tolerance, the tail gate, and, when the selection
# weight has transitions, the sharpness gate and the certificate with its
# combined allowance. Returns the accepted rung, or NULL when the block reaches
# the randomized fallback.
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
      if (tail <= tolerance && (!transitions || .test_cert_resolves(block, orders[k]))) {
        amplitude <- if (transitions) sum(.test_cert_amplitudes(block, orders[k])[["amplitude"]]) else 0
        if (amplitude <= tolerance / 2 &&
            tail + .test_cert_allowance_factor * amplitude <= tolerance / (1 + tolerance)) {
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
# Small blocks of every support structure
# ---------------------------------------------------------------------------- #

# Structures of the selection supports at test size: one axis, two independent
# roots, a root with two leaves, two trees of a root with one leaf, and partial
# overlaps (a tensor). Sharp enough that a coarse rung is far from its limit.
.test_cert_small_group <- function(count, loading, sd, mean, sei = 1) {

  list(count = count, loading = loading, sd = sd, mean = mean, sei = sei)
}

.test_cert_small_blocks <- function(omega = c(1, .2), sd_scale = .3) {

  list(
    # one latent axis
    rank_one = .test_cert_block(list(
      .test_cert_small_group(4L, .8, sd_scale * 1.1, .3)), omega),
    # two independent roots, each with rows on its own axis (a forest of depth zero)
    two_roots = .test_cert_block(list(
      .test_cert_small_group(3L, c(.8, 0), sd_scale * 1.1, .2),
      .test_cert_small_group(2L, c(0, .7), sd_scale * 1.3, -.1)), omega),
    # a shared axis with two leaves, rows on the shared axis only and a free row
    two_leaves = .test_cert_block(list(
      .test_cert_small_group(3L, c(.6, .5, 0), sd_scale * 1.1, .2),
      .test_cert_small_group(2L, c(.6, 0, .4), sd_scale * 1.3, -.1),
      .test_cert_small_group(2L, c(.7, 0, 0), sd_scale * 1.2, .1),
      .test_cert_small_group(1L, c(0, 0, 0), sd_scale * 1, 0)), omega),
    # two independent trees, each a root with one leaf
    two_trees = .test_cert_block(list(
      .test_cert_small_group(3L, c(.6, .5, 0, 0), sd_scale * 1.1, .2),
      .test_cert_small_group(2L, c(0, 0, .5, .4), sd_scale * 1.3, -.1)), omega),
    # a group on both axes with rows on either axis alone: partial overlaps, no
    # forest, the rank-two tensor rule
    crossing = .test_cert_block(list(
      .test_cert_small_group(3L, c(.7, .3), sd_scale * 1.1, .2),
      .test_cert_small_group(2L, c(.5, 0), sd_scale * 1.2, -.1),
      .test_cert_small_group(2L, c(0, .6), sd_scale * 1.3, .1)), omega)
  )
}

# ---------------------------------------------------------------------------- #
# Native refinement tier
# ---------------------------------------------------------------------------- #

# The decision of the tier recomputed from the probe: the sharpness gate at the
# finest rung, the M = 4 average of the finest and previous affordable rung
# (orders 95 and 63 at rank two and for forests of depth one), the shift check,
# the tail gate and the order check. `checks` is the decision without the gate.
.test_cert_tier <- function(block, tolerance = .005, finest = 95L, previous = 63L) {

  fine <- .test_cert_probe(block, finest, tolerance = tolerance)
  prior <- .test_cert_probe(block, previous, tolerance = tolerance)
  shift_change <- abs(expm1(fine[["refined_4"]] - fine[["refined_2"]]))
  order_change <- abs(expm1(fine[["refined_4"]] - prior[["refined_4"]]))
  tail <- .test_cert_tail_bound(block, finest, fine[["refined_4"]])
  checks <- shift_change <= tolerance / 2 && tail <= tolerance && order_change <= tolerance / 2
  list(
    value = fine[["refined_4"]], shift_change = shift_change, order_change = order_change,
    tail = tail, resolves = fine[["resolves_tier"]], checks = checks,
    accepted = fine[["resolves_tier"]] && checks
  )
}

# The accepted order of the certificate ladder recomputed from the probe's rung
# values: finite, the two successive changes and the tail gate within the
# tolerance, the sharpness gate and the certificate; NA when every affordable
# rung is rejected.
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
        probe[["tail_bound"]] <= tolerance && probe[["resolves_ladder"]] && probe[["accepted"]]) {
      return(order)
    }
    previous <- value
    previous_change <- change
  }
  NA_integer_
}

# ---------------------------------------------------------------------------- #
# Deferred geometries that reach the randomized fallback
# ---------------------------------------------------------------------------- #

# Native fallback runs. Every block here is certified by neither the ladder nor the
# refinement tier, so it
# reaches the randomized fallback; the fixed design seeds make the runs repeatable.
# The parameters come from the research battery (seeds 41 and 21) and random searches
# over forests (a shared axis and two leaves) and rank-two partial overlaps (tensor
# rule). The optimized-mode proposals of the previous build accept a minority of
# them and err by 10-1000 MCSE on several. Blocks whose groups all load on both axes
# (parallel_80, band_125, crossing_20) are forests of depth one: a root that owns no
# rows and a leaf.
.test_cert_proposal_cases <- function() {

  g <- function(count, loading, sd, mean, sei = 1) list(count = count, loading = loading, sd = sd, mean = mean, sei = sei)
  list(
    parallel_80 = list(kind = "factor", reference = .test_cert_reference_parallel, block = .test_cert_block(list(
      g(12L, c(0.014692097199723289, -0.8616553723389665), 0.13653164601012455, -5.6000270018551461, 0.34873594478420489)),
      c(1, 0.050000000000000003), c(0, -Inf))),
    band_125 = list(kind = "factor", reference = .test_cert_reference_parallel, block = .test_cert_block(list(
      g(44L, c(1.4323285021004943, 0.083058174726865039), 0.060454506795780072, 1.8231652047112927, 0.55342308823237685)),
      c(1, 1.1010918456055847, 1), c(1, -1, -Inf))),
    crossing_20 = list(kind = "factor", reference = .test_cert_reference_two_groups, block = .test_cert_block(list(
      g(77L, c(-0.91463423771209651, 0.78986739681623297), 0.28436590998117955, -2.3718989917066349),
      g(18L, c(0.015058325926115159, -0.30020325671357212), 0.019704307952072554, 1.8604332933661574)),
      c(1, 0.050000000000000003), c(1.96, -Inf))),
    # partial overlaps (a group on both axes, one group on each axis alone): the tensor rule
    tensor_45 = list(kind = "factor", reference = .test_cert_reference_rank_two, pinned = -13.835385492311222,
      block = .test_cert_block(list(
        g(27L, c(0.82706647334157291, 0.65077649884120636), 0.054795022412045077, -0.043517232008220219),
        g(11L, c(-1.3674417541470836, 0), 0.064991112511571233, 1.1165970716872602),
        g(22L, c(0, 1.3796914536846054), 0.057839652420356456, 0.27615990856989825)),
        c(1, 0.20000000000000001), c(1.96, -Inf))),
    tensor_77 = list(kind = "factor", reference = .test_cert_reference_rank_two, pinned = -2.5592954958880094,
      block = .test_cert_block(list(
        g(12L, c(-0.34458291062427782, -0.83820528250345905), 0.41598330888330703, 1.3861580619958969),
        g(25L, c(-0.69611133285596238, 0), 0.081507582329821651, 2.515650777258446),
        g(17L, c(0, -1.8009762478430924), 0.13199939157250279, 0.82995553693087865)),
        c(1, 0.01), c(1.96, -Inf))),
    forest_74 = list(kind = "factor", reference = .test_cert_reference_forest, block = .test_cert_block(list(
      g(65L, c(0.060605698492417107, 0.9655385554788849, 0), 0.57052084244334667, -2.3719209529207341),
      g(56L, c(0.060605698492417107, 0, 1.6252115421914184), 0.040203888654504902, -4.6092623941950031)),
      c(1, 0.20000000000000001), c(0, -Inf))),
    forest_84 = list(kind = "factor", reference = .test_cert_reference_forest, block = .test_cert_block(list(
      g(127L, c(0.12399440036971428, 1.9003074081263616, 0), 0.0989997504638504, -2.5188992921523741),
      g(30L, c(0.12399440036971428, 0, 1.1540407331499025), 0.030234154411339583, 1.6193386859839336)),
      c(1, 0.050000000000000003), c(1.645, -Inf))),
    # The original cluster_733 and cluster_189 now pass the cluster refinement
    # tier and are retained in its oracle tests. Halving their residual SDs
    # leaves genuinely sharp blocks for these randomized-proposal checks.
    cluster_733_sharp = list(kind = "cluster", reference = .test_cert_reference_parallel, block = .test_cert_block(list(
      g(4L, 1.8707955473336113, 0.042233265807645956 / 2, -6.8303662251056725, 0.22565989842880826)),
      c(1, 0.20000000000000001), c(0, -Inf))),
    cluster_189_sharp = list(kind = "cluster", reference = .test_cert_reference_parallel, block = .test_cert_block(list(
      g(100L, 0.65477525608481191, 0.026923363058843731 / 2, -2.5196777242881367, 0.20064200965816639)),
      c(1, 0.050000000000000003), c(1.645, -Inf)))
  )
}
