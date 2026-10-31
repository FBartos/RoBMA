context("Native cluster input validation")

test_that("cluster selection batches reject invalid observed values and bins", {

  rule <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  args <- list(c(.1, -.1), matrix(c(0, 0), 1L), matrix(c(1, 1), 1L),
    matrix(c(.2, .2), 1L), c(1, 1), matrix(c(1, .2), 1L),
    c(0, -Inf), c(Inf, 0), c(1L, 2L), 1L, TRUE, 1L,
    rule$nodes, rule$log_weights, as.double(rule$orders),
    rep(.5, 2L * 8L * 2L), 4L, 8L, 2L, .005, FALSE, 0L)
  evaluate <- function(values) {

    do.call(.Call, c(list("RoBMA_selnorm_cluster_step_loglik_batch"), values,
                     list(PACKAGE = "RoBMA")))
  }
  for (bin in c(NA_integer_, 0L, 3L)) {
    changed <- args
    changed[[9L]][[1L]] <- bin
    expect_error(evaluate(changed), "Observed selection inputs are invalid")
  }
  for (entry in c(1L, 5L)) {
    changed <- args
    changed[[entry]][[1L]] <- NaN
    expect_error(evaluate(changed), "Observed selection inputs are invalid")
  }
})


test_that("cluster ingress validates orders and Gaussian rows before workers", {

  rule <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  args <- list(c(.1, -.1), matrix(c(0, 0), 1L), matrix(c(1, 1), 1L),
    matrix(c(.2, .2), 1L), c(1, 1), matrix(c(1, .2), 1L),
    c(0, -Inf), c(Inf, 0), c(1L, 2L), 1L, TRUE, 0L,
    rule$nodes, rule$log_weights, as.double(rule$orders),
    rep(.5, 32L), 4L, 8L, 2L, .005, FALSE, 0L)
  evaluate <- function(values) {
    do.call(.Call, c(list("RoBMA_selnorm_cluster_step_loglik_batch"), values,
      list(PACKAGE = "RoBMA")))
  }
  reference <- evaluate(args)
  expect_true(all(is.finite(reference$log_density)))
  for (order in c(Inf, NaN, 2.5, .Machine$integer.max + 1, length(rule$nodes) + 1)) {
    changed <- args
    changed[[15L]][1L] <- order
    expect_error(evaluate(changed),
      "'quadrature_orders' must contain increasing positive integers matching the quadrature rules.",
      fixed = TRUE)
  }
  for (entry in 2:4) {
    changed <- args
    changed[[entry]][1L] <- Inf
    expect_error(evaluate(changed), if (entry == 4L)
      "Cluster selection loadings must be finite." else "Cluster selection location inputs are invalid.",
      fixed = TRUE)
  }
  changed <- args
  changed[[3L]][1L] <- 0
  expect_error(evaluate(changed), "Cluster selection location inputs are invalid.", fixed = TRUE)
  expect_identical(evaluate(args), reference)
})
