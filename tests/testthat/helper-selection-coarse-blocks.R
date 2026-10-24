# Fixtures and independent references for the coarse-corrected sampler on the
# cluster (rank one) and factor (rank 2..8) selection blocks. A block is the
# parameter set of one dselnorm_cluster_step / dselnorm_factor_step node;
# `.coarse_probe()` evaluates its surrogate and exact density through the
# distribution classes the sampler uses (RoBMA_selnorm_surrogate_probe).

# Quadrature schedule in the layout of the model data (nodes, log weights, orders).
.coarse_rules <- function(orders) {

  RoBMA:::.selection_joint_cluster_quadrature_rules(as.integer(orders))
}


# Parameter list of a block in the order of the distribution's arguments. The
# loading is a vector (cluster) or a rows x rank matrix (factor).
.coarse_block <- function(kind, mean, residual_sd, loading, sei, omega,
                          z_lower = c(0, -1e300), z_upper = c(1e300, 0),
                          obs_bin = rep(1L, length(mean)), x = rep(0, length(mean)),
                          rules = .coarse_rules(c(7, 15, 31, 63, 127)),
                          sign = 1, telescope = 1, mode = 1, vector_rule = 0,
                          initial_points = 4, max_points = 8, scrambles = 2,
                          tolerance = .005, seed = 1L) {

  loading <- as.matrix(loading)
  rank    <- if (kind == "cluster") 1L else ncol(loading)
  stopifnot(kind == "factor" || ncol(loading) == 1L, nrow(loading) == length(mean))
  qmc <- BayesTools::selection_qmc_design(2L * rank, max_points, scrambles, seed)
  parameters <- list(
    as.double(mean), as.double(residual_sd), as.double(loading), as.double(sei),
    as.double(omega), as.double(z_lower), as.double(z_upper), as.double(obs_bin),
    as.double(sign), as.double(telescope), as.double(mode),
    as.double(rules[["nodes"]]), as.double(rules[["log_weights"]]),
    as.double(rules[["orders"]]), as.double(qmc), as.double(initial_points),
    as.double(max_points), as.double(scrambles), as.double(tolerance),
    as.double(vector_rule))
  list(kind = kind, x = as.double(x), parameters = parameters, rank = rank,
       rules = rules)
}


# Replace one distribution argument by name (a 1-based index into the list).
.coarse_block_set <- function(block, index, value) {

  block[["parameters"]][[index]] <- as.double(value)
  block
}


# `unchecked` also evaluates blocks whose only refused controls are non-finite
# quadrature nodes or log weights (test use only); any other invalid control is an error.
.coarse_probe <- function(block, max_rules = 3L, unchecked = FALSE) {

  .Call("RoBMA_selnorm_surrogate_probe", block[["kind"]], block[["x"]],
        block[["parameters"]], c(.05, .05, .01, as.double(max_rules)), unchecked,
        PACKAGE = "RoBMA")
}


# Probabilists' Gauss-Hermite rule by Golub-Welsch, independent of the package's
# tables: nodes are the eigenvalues of the Jacobi matrix with off-diagonals
# sqrt(1:(n - 1)) and the weights the squared first eigenvector components.
.coarse_reference_rule <- function(n) {

  if (n == 1L) return(list(nodes = 0, weights = 1))
  jacobi <- matrix(0, n, n)
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- sqrt(seq_len(n - 1L))
  jacobi[cbind(2:n, seq_len(n - 1L))] <- sqrt(seq_len(n - 1L))
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(nodes = decomposition$values, weights = decomposition$vectors[1L, ]^2)
}


# Probability that a row with latent mean `m` (vector over quadrature points),
# residual SD `sd` and selection SE `sei` is observed under the step weights.
.coarse_reference_probability <- function(m, sd, sei, omega, z_lower, z_upper, sign) {

  lower <- ifelse(z_lower <= -1e299, -Inf, z_lower)
  upper <- ifelse(z_upper >= 1e299, Inf, z_upper)
  total <- 0
  for (bin in seq_along(omega)) {
    probability <- if (sign == 1) {
      stats::pnorm(upper[bin] * sei, m, sd) - stats::pnorm(lower[bin] * sei, m, sd)
    } else {
      stats::pnorm(-lower[bin] * sei, m, sd) - stats::pnorm(-upper[bin] * sei, m, sd)
    }
    total <- total + omega[bin] * probability
  }
  total
}


# log A_j: the log of the tensor Gauss-Hermite rule of order `order` per latent
# axis, summed over every point of the product grid (no nesting, no
# compaction), of the product of the row probabilities.
.coarse_reference_log_normalizer <- function(block, order) {

  p      <- block[["parameters"]]
  mean   <- p[[1L]]
  sd     <- p[[2L]]
  k      <- length(mean)
  rank   <- block[["rank"]]
  loading <- matrix(p[[3L]], k, rank)
  rule   <- .coarse_reference_rule(order)
  grid   <- as.matrix(expand.grid(rep(list(seq_len(order)), rank)))
  weight <- Reduce(`*`, lapply(seq_len(rank), function(axis) rule$weights[grid[, axis]]))
  latent <- matrix(rule$nodes[grid], nrow(grid), rank)
  product <- rep(1, nrow(grid))
  for (i in seq_len(k)) {
    product <- product * .coarse_reference_probability(
      mean[i] + as.numeric(latent %*% loading[i, ]), sd[i], p[[4L]][i],
      p[[5L]], p[[6L]], p[[7L]], p[[9L]])
  }
  log(sum(weight * product))
}


# Gaussian log density of x with covariance diag(sd^2) + L L' plus the observed
# log weights: the numerator of the exact density.
.coarse_reference_numerator <- function(block) {

  p       <- block[["parameters"]]
  loading <- matrix(p[[3L]], length(p[[1L]]), block[["rank"]])
  .coarse_gaussian_log_density(block[["x"]], p[[1L]], p[[2L]], loading) + sum(log(p[[5L]][p[[8L]]]))
}


# The order of rule j of the block's schedule.
.coarse_rule_order <- function(block, max_rules = 3L) {

  orders <- block[["rules"]][["orders"]]
  orders[min(max_rules, length(orders))]
}
