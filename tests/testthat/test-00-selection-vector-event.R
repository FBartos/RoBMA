context("Conditional Gaussian vector selection events")

.vector_event_fixture <- function(two_sided = FALSE, weights = c(2, 0, .7)) {

  steps <- c(.05, .3)
  cutoffs <- stats::qnorm(c(0, steps / if (two_sided) 2 else 1, 1),
                         lower.tail = FALSE)
  if (two_sided) {
    cutoffs <- c(cutoffs[1:3], 0, -rev(cutoffs[1:3]))
    weights <- c(weights, rev(weights))
  }
  list(lower = tail(cutoffs, -1L), upper = head(cutoffs, -1L),
       omega = weights, steps = steps, two_sided = two_sided)
}

.vector_event_mass <- function(mean, covariance, se, event, rule,
                               sign = 1L, lower = rep(-Inf, length(mean)),
                               upper = rep(Inf, length(mean)), mode = 1L,
                               quadrature = FALSE) {

  K <- length(mean)
  design <- if (K == 1L) numeric() else BayesTools::selection_qmc_design(
    dimensions = 2L * K, points = 4096L, scrambles = 8L, seed = 719L
  )
  .Call(
    "RoBMA_selnorm_gaussian_event_mass_batch",
    matrix(as.double(mean), 1L),
    matrix(covariance[lower.tri(covariance, diag = TRUE)], 1L),
    as.double(se), matrix(event$omega, 1L), event$lower, event$upper,
    as.integer(sign), as.integer(mode), as.integer(rule),
    matrix(as.double(lower), 1L), matrix(as.double(upper), 1L),
    as.double(design), 4096L, 8L, NULL,
    if (quadrature) .selection_joint_cluster_quadrature_rules(
      SELNORM_CLUSTER_QUADRATURE_ORDERS) else NULL,
    .005, PACKAGE = "RoBMA"
  )
}

# Independent reference: enumerate the disjoint one/two-dimensional outcome
# rectangles, and integrate the bivariate Gaussian by base-R conditioning.
# Compute the selected bin from the actual p-value of an interior point.
.vector_event_reference <- function(mean, covariance, se, event, rule,
                                    sign = 1L, lower = rep(-Inf, length(mean)),
                                    upper = rep(Inf, length(mean))) {

  K <- length(mean)
  stopifnot(K %in% 1:2)
  bins <- expand.grid(rep(list(seq_along(event$omega)), K))
  total <- 0
  for (row in seq_len(nrow(bins))) {
    index <- as.integer(bins[row, ])
    lo <- pmax(lower, if (sign == 1L) event$lower[index] * se else
      -event$upper[index] * se)
    hi <- pmin(upper, if (sign == 1L) event$upper[index] * se else
      -event$lower[index] * se)
    if (any(lo >= hi)) next
    if (rule == 0L) {
      weight <- prod(event$omega[index])
    } else {
      midpoint <- (event$lower[index] + event$upper[index]) / 2
      midpoint[is.infinite(event$upper[index])] <- event$lower[index][
        is.infinite(event$upper[index])] + 1
      midpoint[is.infinite(event$lower[index])] <- event$upper[index][
        is.infinite(event$lower[index])] - 1
      p <- if (event$two_sided) 2 * stats::pnorm(-abs(midpoint)) else
        stats::pnorm(midpoint, lower.tail = FALSE)
      weight <- event$omega[[findInterval(min(p), c(0, event$steps, 1),
                                         rightmost.closed = TRUE)]]
    }
    if (weight == 0) next
    probability <- if (K == 1L) {
      stats::pnorm(hi, mean, sqrt(covariance[[1L]])) -
        stats::pnorm(lo, mean, sqrt(covariance[[1L]]))
    } else {
      variance <- covariance[2L, 2L] - covariance[2L, 1L]^2 / covariance[1L, 1L]
      stats::integrate(function(x) {
        conditional <- mean[[2L]] + covariance[2L, 1L] / covariance[1L, 1L] *
          (x - mean[[1L]])
        stats::dnorm(x, mean[[1L]], sqrt(covariance[1L, 1L])) *
          (stats::pnorm(hi[[2L]], conditional, sqrt(variance)) -
             stats::pnorm(lo[[2L]], conditional, sqrt(variance)))
      }, lo[[1L]], hi[[1L]], rel.tol = 1e-10, abs.tol = 1e-12)$value
    }
    total <- total + weight * probability
  }
  total
}

test_that("best weights use the smallest p bin with one vector normalizer", {

  event <- list(lower = c(0, -Inf), upper = c(Inf, 0), omega = c(1, .2))
  best <- .vector_event_mass(c(0, 0), diag(2L), c(1, 1), event, 1L)
  product <- .vector_event_mass(c(0, 0), diag(2L), c(1, 1), event, 0L)
  expect_equal(exp(best$log_mass), .8, tolerance = 1e-14)
  expect_equal(exp(product$log_mass), .36, tolerance = 1e-14)
  expect_identical(best$relative_mcse, 0)
  event$omega <- c(.2, 1)
  expect_equal(exp(.vector_event_mass(c(0, 0), diag(2L), c(1, 1), event,
                                     1L)$log_mass), .4, tolerance = 1e-14)
  event$omega <- c(1, 0)
  rare <- .vector_event_mass(c(-15, -15), diag(2L), c(1, 1), event, 1L)
  expect_equal(rare$log_mass, log(2) + stats::pnorm(15, lower.tail = FALSE,
                                                  log.p = TRUE), tolerance = 1e-13)
})

test_that("Gaussian event kernels agree with independent p-bin integrals", {

  for (two_sided in c(FALSE, TRUE)) {
    event <- .vector_event_fixture(two_sided)
    for (K in 1:2) {
      mean <- c(-.23, .48)[seq_len(K)]
      se <- c(.7, 1.2)[seq_len(K)]
      covariance <- matrix(c(.8, -.19, -.19, 1.3), 2L)[seq_len(K), seq_len(K),
                                                        drop = FALSE]
      for (rule in c(0L, if (two_sided) 2L else 1L)) {
        for (sign in c(-1L, 1L)) {
          expected <- .vector_event_reference(mean, covariance, se, event, rule, sign)
          actual <- .vector_event_mass(mean, covariance, se, event, rule, sign)
          expect_equal(exp(actual$log_mass), expected,
                       tolerance = if (K == 1L) 1e-12 else 3e-4)
          expect_lt(actual$relative_mcse, 3e-4)
          box_lower <- c(-.9, -.3)[seq_len(K)]
          box_upper <- c(.6, 1.8)[seq_len(K)]
          expected_box <- .vector_event_reference(mean, covariance, se, event,
                                                  rule, sign, box_lower, box_upper)
          actual_box <- .vector_event_mass(mean, covariance, se, event, rule,
                                           sign, box_lower, box_upper)
          expect_equal(exp(actual_box$log_mass), expected_box,
                       tolerance = if (K == 1L) 1e-12 else 3e-4)
        }
      }
    }
  }
})

test_that("event endpoints and constant weights preserve Gaussian identities", {

  event <- .vector_event_fixture(TRUE, c(1.7, 1.7, 1.7))
  mean <- c(.2, -.1)
  covariance <- matrix(c(1, .3, .3, .8), 2L)
  for (rule in 0:2) {
    actual <- .vector_event_mass(mean, covariance, c(.8, .6), event, rule)
    expect_equal(actual$log_mass, log(1.7) * if (rule == 0L) 2 else 1,
                 tolerance = 1e-14)
    empty <- .vector_event_mass(mean, covariance, c(.8, .6), event, rule,
                                lower = c(0, -Inf), upper = c(0, Inf))
    expect_identical(empty$log_mass, -Inf)
    normal <- .vector_event_mass(mean, covariance, c(.8, .6), event, rule, mode = 0L)
    expect_identical(normal$log_mass, 0)
  }
})

test_that("declared rank-one events retain zero-residual and point-mass boundaries", {

  event <- list(lower = c(0, -Inf), upper = c(Inf, 0), omega = c(1, .2))
  evaluate <- function(loading, rule, lower = c(-Inf, -Inf), upper = c(Inf, Inf),
                       mean = c(0, 0)) {
    .Call(
      "RoBMA_selnorm_gaussian_event_mass_batch", matrix(mean, 1L), NULL,
      c(1, 1), matrix(event$omega, 1L), event$lower, event$upper,
      1L, 1L, as.integer(rule), matrix(lower, 1L), matrix(upper, 1L),
      numeric(), 8L, 2L, matrix(loading, 1L), NULL, .005, PACKAGE = "RoBMA"
    )
  }
  for (rule in 0:1) {
    same <- evaluate(c(1, 1), rule)
    expect_equal(exp(same$log_mass), if (rule == 0L) .52 else .6,
                 tolerance = 1e-14)
    opposite <- evaluate(c(1, -1), rule)
    expect_equal(exp(opposite$log_mass), if (rule == 0L) .2 else 1,
                 tolerance = 1e-14)
    expect_identical(opposite$relative_mcse, 0)
    restricted <- evaluate(c(1, -1), rule, upper = c(.5, .2))
    expect_equal(exp(restricted$log_mass),
                 (stats::pnorm(.5) - stats::pnorm(-.2)) * if (rule == 0L) .2 else 1,
                 tolerance = 1e-14)
    point <- evaluate(c(0, 0), rule, lower = c(0, 0), upper = c(0, 0))
    expect_identical(point$log_mass, 0)
    excluded <- evaluate(c(0, 0), rule, lower = c(.1, 0))
    expect_identical(excluded$log_mass, -Inf)
  }
})

test_that("factor and cluster normalizers share the conditional Gaussian event law", {

  mean <- c(-.1, .24)
  y <- c(.2, -.3)
  se <- c(.6, .8)
  residual <- c(.7, .9)
  loading <- matrix(c(.3, .25, -.1, .2), 2L)
  cluster_rules <- .selection_joint_cluster_quadrature_rules(
    SELNORM_CLUSTER_QUADRATURE_ORDERS
  )
  factor_rules <- .selection_joint_factor_quadrature_rules(2L)[["2"]]
  qmc <- as.double(BayesTools::selection_qmc_design(
    dimensions = 4L, points = 2048L, scrambles = 8L, seed = 417L
  ))
  cluster_qmc <- as.double(BayesTools::selection_qmc_design(
    dimensions = 2L, points = 2048L, scrambles = 8L, seed = 417L
  ))
  for (two_sided in c(FALSE, TRUE)) {
    event <- .vector_event_fixture(two_sided)
    rule <- if (two_sided) 2L else 1L
    bins <- vapply(y / se, function(value) which(value >= event$lower)[[1L]], integer(1L))
    for (method in c("cluster", "factor", "dense")) {
      covariance <- diag(residual^2) + if (method == "cluster")
        tcrossprod(loading[, 1L, drop = FALSE]) else tcrossprod(loading)
      common <- list(se, matrix(event$omega, 1L), event$lower, event$upper,
                     bins, 1L, TRUE, 1L)
      arguments <- switch(method,
        cluster = c(list(y, matrix(mean, 1L), matrix(residual, 1L),
                         matrix(loading[, 1L], 1L)), common,
                    list(cluster_rules$nodes, cluster_rules$log_weights,
                         as.double(cluster_rules$orders), cluster_qmc,
                         256L, 2048L, 8L, .005, TRUE, rule)),
        factor = c(list(y, matrix(mean, 1L), matrix(residual, 1L),
                        matrix(as.double(loading), 1L)), common,
                   list(factor_rules$nodes, factor_rules$log_weights,
                        as.double(factor_rules$orders), as.double(factor_rules$rule_counts),
                        qmc, 256L, 2048L, 8L, .005, TRUE, rule)),
        dense = c(list(y, matrix(mean, 1L),
                       matrix(covariance[lower.tri(covariance, diag = TRUE)], 1L)),
                  common, list(qmc, 2048L, 8L, .005, TRUE, rule, NULL))
      )
      symbol <- paste0("RoBMA_selnorm_", if (method == "dense") "mnorm" else method,
                       "_step_loglik_batch")
      actual <- do.call(.Call, c(list(symbol), arguments, list(PACKAGE = "RoBMA")))
      expected <- .vector_event_reference(mean, covariance, se, event, rule)
      expect_equal(exp(actual$log_normalizer), expected, tolerance = 3e-4)
      p <- if (two_sided) 2 * stats::pnorm(-abs(y / se)) else
        stats::pnorm(y / se, lower.tail = FALSE)
      weight <- event$omega[[findInterval(min(p), c(0, event$steps, 1))]]
      expected_density <- mvtnorm::dmvnorm(y, mean, covariance, log = TRUE) +
        log(weight) - log(expected)
      expect_equal(actual$log_density, expected_density, tolerance = 3e-4)
    }
  }
})

test_that("best response draws retain the full independent publication vector", {

  n <- 6000L
  covariance <- array(0, c(n, 2L, 2L))
  covariance[, 1L, 1L] <- 1
  covariance[, 2L, 2L] <- 1
  set.seed(8173)
  draws <- .Call(
    "RoBMA_selnorm_mnorm_step_rng_batch", matrix(0, n, 2L), covariance,
    c(1, 1), matrix(rep(c(1, .2), each = n), n),
    c(0, -Inf), c(Inf, 0), 1L, 1L, list(1:2), 1000L, 1L,
    PACKAGE = "RoBMA"
  )
  counts <- rowSums(draws$draws > 0)
  observed <- tabulate(counts + 1L, nbins = 3L) / n
  # Independent probability calculation: Gaussian masses (.25,.5,.25)
  # times best weights (.2,1,1), normalized by .8.
  expect_equal(observed, c(.0625, .625, .3125), tolerance = .025)
  expect_true(all(is.finite(draws$draws)))
  for (sign in c(-1L, 1L)) {
    excluded <- .Call(
      "RoBMA_selnorm_mnorm_step_rng_batch", matrix(0, 500L, 2L),
      covariance[seq_len(500L), , , drop = FALSE], c(1, 1),
      matrix(rep(c(0, 1), each = 500L), 500L), c(0, -Inf), c(Inf, 0),
      sign, 1L, list(1:2), 1000L, 1L, PACKAGE = "RoBMA"
    )
    expect_identical(excluded$failure_code, 0L)
    expect_true(all(sign * excluded$draws <= 0))
  }
})

test_that("the Gaussian event facade validates its explicit numerical contract", {

  event <- .vector_event_fixture()
  context <- list(omega = matrix(event$omega, 1L), kernel_mode = 1L,
                  vector_rule = 1L, z_lower = event$lower,
                  z_upper = event$upper, sign = 1L)
  plan <- list(points_per_scramble = 8L, scrambles = 2L, seed = 519L)
  evaluate <- function(...) .selection_gaussian_event_mass(
    mean = matrix(c(0, 0), 1L), selection_se = c(1, 1),
    selection_context = context, execution_plan = plan, ...
  )
  actual <- evaluate(covariance_lower = matrix(c(1, 0, 1), 1L))
  reference <- .vector_event_reference(c(0, 0), diag(2L), c(1, 1), event, 1L)
  expect_equal(exp(actual$log_mass), reference, tolerance = 1e-14)
  expect_error(evaluate(covariance_lower = matrix(c(1, 0, 1), 1L), lower = NA_real_),
               "'lower' must contain numeric Gaussian event limits.", fixed = TRUE)
  expect_error(evaluate(covariance_lower = matrix(c(1, 0, 1), 1L), upper = 1),
               "'upper' must have one value per observation or match 'mean'.", fixed = TRUE)
  expect_error(evaluate(covariance_lower = matrix(c(1, 0, 1), 1L),
                        rank_one_loading = matrix(c(1, 1), 1L)),
               "'rank_one_loading' and 'covariance_lower' are alternative Gaussian covariance representations.",
               fixed = TRUE)
  plan$points_per_scramble <- 1.5
  expect_error(evaluate(covariance_lower = matrix(c(1, 0, 1), 1L)),
               "Gaussian selection event controls are invalid.", fixed = TRUE)
})

test_that("JAGS vector distributions accept the same explicit best event rule", {

  skip_if_not_installed("rjags")
  quadrature <- .selection_joint_factor_quadrature_rules(2L)[["2"]]
  cluster_quadrature <- .selection_joint_cluster_quadrature_rules(
    SELNORM_CLUSTER_QUADRATURE_ORDERS
  )
  qmc <- BayesTools::selection_qmc_design(
    dimensions = 4L, points = 2048L, scrambles = 8L, seed = 519L
  )
  cluster_qmc <- BayesTools::selection_qmc_design(
    dimensions = 2L, points = 2048L, scrambles = 8L, seed = 519L
  )
  for (method in c("cluster", "factor", "mnorm")) {
    for (two_sided in c(FALSE, TRUE)) {
      event <- .vector_event_fixture(two_sided)
      data <- list(
        y = c(.2, -.3), offset = c(-.1, .24), sei = c(.6, .8),
        omega = event$omega,
        z_lower = ifelse(is.infinite(event$lower), -1e300, event$lower),
        z_upper = ifelse(is.infinite(event$upper), 1e300, event$upper),
        bins = vapply(c(.2, -.3) / c(.6, .8), function(value) {
          which(value >= event$lower)[[1L]]
        }, integer(1L)),
        rule = if (two_sided) 2L else 1L
      )
      covariance_arguments <- switch(method,
        cluster = "mu[],residual_sd[],loading[],",
        factor = "mu[],residual_sd[],loading[,],",
        mnorm = "mu[],covariance[],"
      )
      integration_arguments <- switch(method,
        cluster = "nodes[],log_weights[],orders[],qmc[,,],256,2048,8,.005,rule",
        factor = "nodes[],log_weights[],orders[],rule_counts[],qmc[,,],256,2048,8,.005,rule",
        mnorm = "qmc[,,],2048,8,.005,rule,nodes[],log_weights[],orders[]"
      )
      if (method != "mnorm") {
        selected_quadrature <- if (method == "cluster") cluster_quadrature else quadrature
        data <- c(data, selected_quadrature[c("nodes", "log_weights", "orders")])
        data$residual_sd <- c(.7, .9)
        data$loading <- if (method == "cluster") c(.3, .25) else
          matrix(c(.3, .25, -.1, .2), 2L)
      } else {
        data$covariance <- c(.59, .055, .9125)
        data <- c(data, cluster_quadrature[c("nodes", "log_weights", "orders")])
      }
      data$qmc <- if (method == "cluster") cluster_qmc else qmc
      if (method == "factor") data$rule_counts <- quadrature$rule_counts
      syntax <- paste0("model { beta ~ dnorm(0,1)\n",
        "for(i in 1:2) { mu[i] <- beta + offset[i] }\n",
        "y[1:2] ~ dselnorm_", method, "_step(", covariance_arguments,
        "sei[],omega[],z_lower[],z_upper[],bins[],1,1,1,",
        integration_arguments, ") }")
      connection <- textConnection(syntax)
      model <- tryCatch(rjags::jags.model(
        connection, data = data, n.chains = 1L, n.adapt = 0L, quiet = TRUE,
        inits = list(beta = 0, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 519L)
      ), finally = close(connection))
      draws <- rjags::coda.samples(model, "beta", n.iter = 4L, progress.bar = "none")
      expect_s3_class(model, "jags")
      expect_true(all(is.finite(as.matrix(draws))))
    }
  }
})


test_that("full-space event quadrature preserves trivariate Gaussian sign identities", {

  sei <- c(.2, .35, .6)
  covariance <- outer(sei, sei) *
    matrix(c(1, .7, .5, .7, 1, .5, .5, .5, 1), 3L)
  diag(covariance) <- diag(covariance) + .1^2
  threshold <- stats::qnorm(.975)
  mean <- threshold * sei
  correlation <- stats::cov2cor(covariance)
  pair_sum <- sum(asin(correlation[lower.tri(correlation)]))
  for (omega in list(c(1, 0), c(1, 2))) {
    event <- list(lower = c(threshold, -Inf), upper = c(Inf, threshold),
      omega = omega)
    # Centering each response on its original selection threshold reduces the
    # product to c + (J/2) sign(Y_i). Odd sign moments vanish; pair moments are
    # 2*asin(rho)/pi. Candidate SDs differ from the original threshold SEs.
    center <- mean(omega)
    jump <- omega[[1L]] - omega[[2L]]
    reference <- log(center^3 + center * jump^2 * pair_sum / (2 * pi))
    for (sign in c(-1L, 1L)) {
      actual <- .vector_event_mass(sign * mean, covariance, sei, event,
        0L, sign = sign, quadrature = TRUE)
      expect_lt(abs(actual$log_mass - reference), 5e-4)
      expect_identical(actual$relative_mcse, 0)
      expect_gt(actual$relative_quadrature_error, 0)
      expect_lte(actual$relative_quadrature_error, .005)
    }
  }
})

test_that("event quadrature preserves fallbacks and reports its own error", {

  event <- list(lower = c(0, -Inf), upper = c(Inf, 0), omega = c(1, .2))
  for (singular in c(FALSE, TRUE)) {
    covariance <- if (singular) matrix(1, 2L, 2L) else
      matrix(c(1, .4, .4, 1), 2L)
    lower <- if (singular) c(-Inf, -Inf) else c(-.3, -Inf)
    actual <- .vector_event_mass(c(0, 0), covariance, c(1, 1), event,
      0L, lower = lower, quadrature = TRUE)
    original <- .vector_event_mass(c(0, 0), covariance, c(1, 1), event,
      0L, lower = lower)
    expect_identical(actual, original)
    expect_identical(actual$relative_quadrature_error, 0)
  }
  plan <- list(points_per_scramble = 4L, max_points_per_scramble = 4L,
    relative_tolerance = .005)
  expect_error(.selection_joint_checked_event(function(plan, rows) {
    list(log_mass = -.2, relative_mcse = 0, relative_quadrature_error = .02)
  }, plan, "Test event", "Increase the integration budget."),
    "Test event was rejected by diagnostics: relative quadrature error was 0.02. Inspect the integration diagnostics.",
    fixed = TRUE)
})
