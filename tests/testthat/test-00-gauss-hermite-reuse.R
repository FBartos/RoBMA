context("Gauss-Hermite rule reuse")


# The previous scalar recurrence is an independent arithmetic-order reference
# for the vectorized implementation. Evaluate it only once for each rule.
.gh_reuse_reference_log_weights <- function(nodes, n) {

  log_weights <- numeric(length(nodes))
  for (node_index in seq_along(nodes)) {
    node                <- nodes[node_index]
    polynomial_previous <- 1
    polynomial          <- node
    log_scale           <- 0
    log_denominator     <- log1p(polynomial^2)
    if (n > 2L) {
      for (degree in seq_len(n - 2L)) {
        polynomial_next <- (
          node * polynomial - sqrt(degree) * polynomial_previous
        ) / sqrt(degree + 1)
        recurrence_scale <- max(abs(c(
          polynomial_previous, polynomial, polynomial_next
        )))
        if (recurrence_scale > 1e100) {
          polynomial_previous <- polynomial_previous / recurrence_scale
          polynomial          <- polynomial / recurrence_scale
          polynomial_next     <- polynomial_next / recurrence_scale
          log_scale           <- log_scale + log(recurrence_scale)
        }
        if (polynomial_next != 0) {
          log_term <- 2 * (log(abs(polynomial_next)) + log_scale)
          if (log_term > log_denominator) {
            log_denominator <- log_term +
              log1p(exp(log_denominator - log_term))
          } else {
            log_denominator <- log_denominator +
              log1p(exp(log_term - log_denominator))
          }
        }
        polynomial_previous <- polynomial
        polynomial          <- polynomial_next
      }
    }
    log_weights[node_index] <- -log_denominator
  }
  max_log_weight <- max(log_weights)
  log_normalizer <- max_log_weight +
    log(sum(exp(log_weights - max_log_weight)))
  log_weights - log_normalizer
}


.gh_reuse_reference_rule <- function(n) {

  i <- seq_len(n - 1)
  H <- diag(0, n)
  H[cbind(i, i + 1)] <- sqrt(i / 2)
  H[cbind(i + 1, i)] <- sqrt(i / 2)
  nodes <- sqrt(2) * eigen(H, symmetric = TRUE, only.values = TRUE)$values
  polynomial_previous <- rep(1, n)
  weight_denominator <- polynomial_previous^2
  if (n > 1L) {
    polynomial <- nodes
    weight_denominator <- weight_denominator + polynomial^2
    for (degree in seq_len(n - 2L)) {
      polynomial_next <- (
        nodes * polynomial - sqrt(degree) * polynomial_previous
      ) / sqrt(degree + 1)
      weight_denominator <- weight_denominator + polynomial_next^2
      polynomial_previous <- polynomial
      polynomial <- polynomial_next
    }
  }
  weights <- 1 / weight_denominator
  log_weights <- log(weights)
  if (any(!is.finite(log_weights))) {
    log_weights <- .gh_reuse_reference_log_weights(nodes, n)
    weights <- exp(log_weights)
  }
  ord <- order(nodes)
  list(nodes = nodes[ord], weights = weights[ord], log_weights = log_weights[ord])
}


test_that("full-first and nodes-first Hermite rules retain exact values and ordering", {

  cache <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(.gauss_hermite_nodes_cache = cache, .package = "RoBMA")
  orders <- c(1L, 2L, 3L, 7L, 9L, 15L, 31L, 61L, 63L, 127L, 255L, 301L, 511L, 1023L)

  for (n in orders) {
    reference <- .gh_reuse_reference_rule(n)
    full_first <- .gauss_hermite_nodes(n)
    expect_identical(full_first, reference, info = paste("full first, order", n))
    expect_identical(.gauss_hermite_nodes_only(n), reference[["nodes"]])
    expect_identical(.gauss_hermite_nodes(n), reference)

    rm(list = as.character(n), envir = cache)
    nodes_first <- .gauss_hermite_nodes_only(n)
    expect_identical(nodes_first, reference[["nodes"]], info = paste("nodes first, order", n))
    expect_null(get(as.character(n), envir = cache)[["weights"]])
    expect_identical(.gauss_hermite_nodes(n), reference, info = paste("lazy weights, order", n))
    expect_identical(.gauss_hermite_nodes_only(n), nodes_first)
    expect_identical(names(get(as.character(n), envir = cache)),
                     c("nodes", "weights", "log_weights"))
  }
})


test_that("nodes-only requests defer finest-rule weights and reuse the existing cache", {

  cache <- new.env(parent = emptyenv())
  scaled_weights <- .gauss_hermite_log_weights_scaled
  weight_calls <- 0L
  testthat::local_mocked_bindings(
    .gauss_hermite_nodes_cache = cache,
    .gauss_hermite_log_weights_scaled = function(nodes, n) {
      weight_calls <<- weight_calls + 1L
      scaled_weights(nodes, n)
    },
    .package = "RoBMA"
  )

  nodes <- .gauss_hermite_nodes_only(1023L)
  expect_identical(weight_calls, 0L)
  expect_identical(.gauss_hermite_nodes_only(1023L), nodes)
  expect_identical(weight_calls, 0L)
  rule <- .gauss_hermite_nodes(1023L)
  expect_identical(weight_calls, 1L)
  expect_identical(rule[["nodes"]], nodes)
  expect_identical(.gauss_hermite_nodes(1023L), rule)
  expect_identical(.gauss_hermite_nodes_only(1023L), nodes)
  expect_identical(weight_calls, 1L)

  rm(list = ls(cache, all.names = TRUE), envir = cache)
  expect_identical(.gauss_hermite_nodes_only(1023L), nodes)
  expect_identical(weight_calls, 1L)
  expect_identical(.gauss_hermite_nodes(1023L), rule)
  expect_identical(weight_calls, 2L)
})
