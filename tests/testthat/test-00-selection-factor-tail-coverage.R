test_that("factor quadrature agreement cannot hide an uncovered selection tail", {

  # Y_i = -10 + Z + E_i, with independent standard normals. Expanding
  # [w + (1-w) I(Y_1 > 0)] [w + (1-w) I(Y_2 > 0)] gives an independent
  # reference in terms of a univariate marginal and a bivariate tail.
  # Piecewise adaptive integration resolves the latter's displaced mode.
  weight <- 1e-9
  breaks <- c(-Inf, seq(-2, 14, by = 2), Inf)
  joint_tail <- sum(vapply(seq_len(length(breaks) - 1L), function(index) {
    stats::integrate(function(z) {
      stats::dnorm(z) * stats::pnorm(z - 10)^2
    }, lower = breaks[[index]], upper = breaks[[index + 1L]],
    rel.tol = 1e-10, abs.tol = 0, subdivisions = 200L)[["value"]]
  }, numeric(1L)))
  reference <- log(weight^2 + 2 * weight * (1 - weight) *
    stats::pnorm(-10 / sqrt(2)) + (1 - weight)^2 * joint_tail)

  # All deliberately short rules miss the distant mode and agree near w^2.
  # The exterior-mass bound must route these cases to the ordinary QMC
  # fallback. The public error criterion and integration budget are unchanged.
  quadrature <- .selection_joint_cluster_quadrature_rules(c(1L, 2L, 3L))
  for (rank in 1:2) {
    loading <- if (rank == 1L) matrix(1, 2L, 1L) else
      matrix(rep(c(.8, .6), each = 2L), 2L, 2L)
    qmc <- as.double(BayesTools::selection_qmc_design(
      dimensions = 2L * rank, points = 8192L, scrambles = 8L, seed = 173L
    ))
    for (sign in c(-1L, 1L)) {
      common <- list(
        c(0, 0), matrix(rep(-10 * sign, 2L), 1L), matrix(c(1, 1), 1L),
        matrix(as.double(loading * sign), 1L), c(1, 1), matrix(c(1, weight), 1L),
        c(0, -Inf), c(Inf, 0), c(1L, 1L), sign, TRUE, 1L,
        quadrature[["nodes"]], quadrature[["log_weights"]],
        as.double(quadrature[["orders"]])
      )
      if (rank == 2L) common <- c(common, list(3))
      arguments <- c(common, list(qmc, 512L, 8192L, 8L, .005, TRUE, 0L))
      symbol <- paste0("RoBMA_selnorm_", if (rank == 1L) "cluster" else "factor",
                       "_step_loglik_batch")
      observed <- do.call(.Call, c(list(symbol), arguments, list(PACKAGE = "RoBMA")))
      expect_gt(observed[["relative_mcse"]], 0)
      expect_lte(observed[["relative_mcse"]], .005)
      expect_lte(observed[["relative_change"]], .005)
      expect_equal(observed[["log_normalizer"]], reference, tolerance = .005)
      expect_gt(abs(observed[["log_normalizer"]] - 2 * log(weight)), 1)
    }
  }
})

.selection_envelope_test_call <- function(
    covariance, mean, sei, omega, sign = 1L, vector_rule = 0L,
    z_lower = c(stats::qnorm(.975), -Inf), z_upper = c(Inf, stats::qnorm(.975)),
    orders = SELNORM_CLUSTER_QUADRATURE_ORDERS) {

  size <- length(sei)
  yi <- rep(0, size)
  bins <- vapply(sign * yi / sei, function(value) {
    which(value >= z_lower)[[1L]]
  }, integer(1L))
  design <- BayesTools::selection_qmc_design(2L * size, 512L, 8L, seed = 173L)
  .Call("RoBMA_selnorm_mnorm_step_loglik_batch",
    yi, matrix(rep_len(mean, size), 1L),
    matrix(covariance[lower.tri(covariance, diag = TRUE)], 1L), sei,
    matrix(omega, 1L), z_lower, z_upper, bins, sign, TRUE, SELKERNEL_STEP,
    as.double(design), 512L, 8L, .005, TRUE, vector_rule,
    .selection_joint_cluster_quadrature_rules(orders),
    PACKAGE = "RoBMA")
}

test_that("ordinary covariance envelopes preserve full-V density and bounded error", {

  data <- data.frame(vi = c(.04, .09, .16), study = 1L,
                     type = c("a", "a", "b"), estimate = 1:3)
  covariance <- as.matrix(metafor::vcalc(vi, cluster = study, type = type,
    obs = estimate, rho = c(.7, .5), data = data))
  # A genuine nonzero departure from compound structure remains authoritative.
  covariance[1L, 3L] <- covariance[3L, 1L] <- covariance[1L, 3L] + 1e-6
  sei <- sqrt(data[["vi"]])
  mean <- c(.1, .3, -.2)
  weight <- .2
  threshold <- stats::qnorm(.975) * sei
  # Expand the product into eight orthant probabilities. Each Gaussian
  # probability has at most three dimensions and uses the independent TVPACK
  # implementation, rather than the factor quadrature under test.
  reference <- sum(vapply(0:7, function(mask) {
    selected <- which(as.logical(intToBits(mask)[1:3]))
    probability <- if (!length(selected)) 1 else if (length(selected) == 1L) {
      stats::pnorm(threshold[selected], mean[selected], sei[selected], lower.tail = FALSE)
    } else {
      as.numeric(mvtnorm::pmvnorm(lower = threshold[selected],
        upper = rep(Inf, length(selected)), mean = mean[selected],
        sigma = covariance[selected, selected, drop = FALSE],
        algorithm = mvtnorm::TVPACK(abseps = 1e-10)))
    }
    weight^(3L - length(selected)) * (1 - weight)^length(selected) * probability
  }, numeric(1L)))
  for (sign in c(-1L, 1L)) {
    observed <- .selection_envelope_test_call(covariance, mean * sign, sei,
                                              c(1, weight), sign)
    diagnostic <- observed[["integration_diagnostics"]][1L, ]
    expect_equal(diagnostic[["used_covariance_envelope"]], 1)
    expect_lte(diagnostic[["covariance_width"]] +
      2 * diagnostic[["quadrature_change"]] + diagnostic[["tail_bound"]], .005)
    expect_gt(diagnostic[["covariance_width"]], 0)
    expect_lt(abs(observed[["log_normalizer"]] - log(reference)), 5e-4)
    expected_numerator <- mvtnorm::dmvnorm(rep(0, 3L), mean * sign, covariance,
                                          log = TRUE) + 3 * log(weight)
    expect_equal(observed[["log_density"]] + observed[["log_normalizer"]],
                 expected_numerator, tolerance = 1e-12)
  }
  best <- .selection_envelope_test_call(covariance, mean, sei,
                                        c(1, weight), vector_rule = 1L)
  expect_equal(unname(best[["integration_diagnostics"]][1L, 1L]), 0)
})

test_that("ordinary Assink V retains a distant normalizer missed by coarse rules", {

  skip_if_not_installed("metadat")
  data("dat.assink2016", package = "metadat", envir = environment())
  study <- dat.assink2016[dat.assink2016[["study"]] == 11L, , drop = FALSE]
  covariance <- as.matrix(metafor::vcalc(vi, cluster = study, type = deltype,
    obs = esid, rho = c(.7, .5), data = study))
  observed <- .selection_envelope_test_call(covariance, -2.2,
    sqrt(study[["vi"]]), c(1, .001))
  # Independent statmod nested-GH reference: orders 128 and 256 agreed within
  # 1.14e-13 in log A. Orders 16 and 24 falsely agreed while missing 87.8% of A.
  reference <- -149.865411056598
  expect_equal(unname(observed[["integration_diagnostics"]][1L, 1L]), 1)
  expect_lt(abs(observed[["log_normalizer"]] - reference), 5e-4)
  expect_gt(abs(observed[["log_normalizer"]] - -151.970449476753), 2)
})

test_that("monotone covariance bounds cover trivariate arcsine normalizers", {

  sei <- c(.2, .35, .6)
  for (departure in c(-2e-4, 2e-4)) {
    correlation <- matrix(c(1, .7, .5, .7, 1, .5, .5, .5, 1), 3L)
    correlation[1L, 3L] <- correlation[3L, 1L] <- .5 + departure
    covariance         <- outer(sei, sei) * correlation
    actual_correlation <- stats::cov2cor(covariance)
    pair_sum           <- sum(asin(actual_correlation[lower.tri(actual_correlation)]))

    # At zero thresholds, each row weight is c + (J/2) sign(Y_i). Odd sign
    # moments vanish and E[sign(Y_i) sign(Y_j)] = 2*asin(rho_ij)/pi. This
    # reference is independent of covariance envelopes and quadrature.
    # Weights above one check the dimension-dependent factor in the bound.
    for (omega in list(c(3, 1), c(1, 3))) {
      center    <- mean(omega)
      jump      <- omega[[1L]] - omega[[2L]]
      reference <- center^3 + center * jump^2 * pair_sum / (2 * pi)
      observed <- .selection_envelope_test_call(
        covariance, 0, sei, omega, z_lower = c(0, -Inf), z_upper = c(Inf, 0)
      )
      diagnostic <- observed[["integration_diagnostics"]][1L, ]
      combined   <- diagnostic[["covariance_width"]] +
        2 * diagnostic[["quadrature_change"]] + diagnostic[["tail_bound"]]

      expect_gt(diagnostic[["covariance_width"]], 0)
      expect_lte(combined, .005)
      # Use the estimate's denominator, as the diagnostics do. The allowance
      # covers floating-point evaluation of the independent closed form.
      relative_error <- abs(expm1(log(reference) - observed[["log_normalizer"]]))
      expect_lte(relative_error, combined + 1e-12)

      # Y=0 is in the first bin. The numerator retains the full supplied V.
      numerator <- mvtnorm::dmvnorm(rep(0, 3L), sigma = covariance, log = TRUE) +
        3 * log(omega[[1L]])
      expect_equal(observed[["log_density"]] + observed[["log_normalizer"]],
                   numerator, tolerance = 1e-12)
    }
  }
})


test_that("absolute covariance bounds cover nonmonotone and two-sided products", {

  sei <- c(.2, .35, .6)
  mean <- c(.06, -.15, .3)
  corners <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 3L)))
  combinations <- as.matrix(expand.grid(rep(list(1:3), 3L)))
  cases <- list(
    list(omega = c(3, 0, 2), z_lower = c(.7, -.5, -Inf),
      z_upper = c(Inf, .7, -.5)),
    list(omega = c(1, .2, 1), z_lower = c(1, -1, -Inf),
      z_upper = c(Inf, 1, -1))
  )
  for (departure in c(-2e-4, 2e-4)) {
    correlation <- matrix(c(1, .7, .5, .7, 1, .5, .5, .5, 1), 3L)
    correlation[1L, 3L] <- correlation[3L, 1L] <- .5 + departure
    covariance <- outer(sei, sei) * correlation

    # TVPACK independently evaluates semi-infinite Gaussian rectangles in
    # up to three dimensions. Inclusion-exclusion gives each finite rectangle;
    # weighting all 27 rectangles then gives the product normalizer directly.
    # This reference does not use factorization or a signed weight expansion.
    orthant <- function(upper) {
      if (any(upper == -Inf)) return(0)
      finite <- which(is.finite(upper))
      if (!length(finite)) return(1)
      if (length(finite) == 1L) {
        return(stats::pnorm(upper[finite], mean[finite], sei[finite]))
      }
      as.numeric(mvtnorm::pmvnorm(lower = rep(-Inf, length(finite)),
        upper = upper[finite], mean = mean[finite],
        sigma = covariance[finite, finite, drop = FALSE],
        algorithm = mvtnorm::TVPACK(abseps = 1e-10)))
    }
    for (case in cases) {
      omega <- case[["omega"]]
      reference <- sum(vapply(seq_len(nrow(combinations)), function(index) {
        bins <- combinations[index, ]
        weight <- prod(omega[bins])
        if (weight == 0) return(0)
        lower <- case[["z_lower"]][bins] * sei
        upper <- case[["z_upper"]][bins] * sei
        probability <- sum(vapply(seq_len(nrow(corners)), function(corner) {
          use_upper <- corners[corner, ]
          (-1)^sum(!use_upper) * orthant(ifelse(use_upper, upper, lower))
        }, numeric(1L)))
        weight * probability
      }, numeric(1L)))
      for (sign in c(-1L, 1L)) {
        observed <- .selection_envelope_test_call(covariance, mean * sign,
          sei, omega, sign, z_lower = case[["z_lower"]],
          z_upper = case[["z_upper"]])
        diagnostic <- observed[["integration_diagnostics"]][1L, ]
        combined <- diagnostic[["covariance_width"]] +
          2 * diagnostic[["quadrature_change"]] + diagnostic[["tail_bound"]]
        relative_error <- abs(expm1(log(reference) - observed[["log_normalizer"]]))
        expect_equal(diagnostic[["used_covariance_envelope"]], 1)
        expect_gt(diagnostic[["covariance_width"]], 0)
        expect_lte(combined, .005)
        # TVPACK is requested at 1e-10 absolute accuracy. The allowance covers
        # the weighted rectangle sums and their division by the normalizer.
        expect_lte(relative_error, combined + 1e-7)
        numerator <- mvtnorm::dmvnorm(rep(0, 3L), mean * sign,
          covariance, log = TRUE) + 3 * log(omega[[2L]])
        expect_equal(observed[["log_density"]] + observed[["log_normalizer"]],
          numerator, tolerance = 1e-12)
      }
    }
  }
})


test_that("a loose nonmonotone bound cannot use a covariance sandwich", {

  sei <- c(.2, .35, .6)
  correlation <- matrix(c(1, .7, .50002, .7, 1, .5, .50002, .5, 1), 3L)
  covariance <- outer(sei, sei) * correlation
  # Almost every candidate lies in the middle bin, so the lower and upper
  # covariance integrals agree numerically. The absolute Price bound ignores
  # threshold locations and exceeds .005 relative to this tiny normalizer.
  # Their apparent agreement cannot certify a nonmonotone covariance sandwich.
  # TV(weights)^2 = 1.8^2 gives a relative bound of about .012; replacing
  # total variation by the range .9 would wrongly reduce it below .005.
  observed <- .selection_envelope_test_call(covariance, 0, sei,
    c(1, .1, 1), z_lower = c(8, -8, -Inf), z_upper = c(Inf, 8, -8))
  expect_equal(unname(observed[["integration_diagnostics"]][1L,
    "used_covariance_envelope"]), 0)
  # A union bound gives at most 6*pnorm(-8) mass outside the middle bins;
  # its contribution relative to .1^3 is less than 4e-12.
  expect_lt(abs(observed[["log_normalizer"]] - 3 * log(.1)), 1e-10)
})


test_that("Assink two-sided defaults agree with independent adaptive integration", {

  skip_if_not_installed("metadat")
  data("dat.assink2016", package = "metadat", envir = environment())
  study <- dat.assink2016[dat.assink2016[["study"]] == 11L, , drop = FALSE]
  sampling <- as.matrix(metafor::vcalc(vi, cluster = study, type = deltype,
    obs = esid, rho = c(.7, .5), data = study))
  sei <- sqrt(study[["vi"]])
  types <- as.integer(factor(study[["deltype"]]))
  expect_equal(tabulate(types), c(16L, 6L))
  # The two type effects have covariance [[.7,.5],[.5,.7]], with independent
  # residual variance .3*vi. Check every supplied covariance entry before using
  # this separate two-factor representation in the numerical reference.
  type_covariance <- matrix(c(.7, .5, .5, .7), 2L)
  represented <- outer(sei, sei) * type_covariance[types, types] +
    diag(.3 * sei^2)
  expect_lt(max(abs(sampling - represented)), 1e-14)
  tau <- .1
  covariance <- sampling + diag(tau^2, length(sei))
  conditional_sd <- sqrt(.3 + tau^2 / sei^2)
  cases <- list(
    list(weights = c(1, .001), steps = .05),
    list(weights = c(1, .1, .001), steps = c(.05, .1))
  )
  for (case in cases) {
    thresholds <- stats::qnorm(1 - case[["steps"]] / 2)
    product_weight <- function(h, type) {
      indices <- which(types == type)
      mean_matrix <- matrix(h, length(indices), length(h), byrow = TRUE)
      q <- matrix(tail(case[["weights"]], 1L), length(indices), length(h))
      for (index in seq_along(thresholds)) {
        tail_probability <- stats::pnorm((-thresholds[[index]] - mean_matrix) /
            conditional_sd[indices]) + stats::pnorm((mean_matrix - thresholds[[index]]) /
            conditional_sd[indices])
        q <- q + (case[["weights"]][[index]] - case[["weights"]][[index + 1L]]) *
          tail_probability
      }
      exp(colSums(log(q)))
    }
    largest_inner_error <- 0
    reference <- stats::integrate(function(x) {
      vapply(x, function(value) {
        h_general <- sqrt(.7) * value
        inner <- stats::integrate(function(z) {
          h_overt <- (.5 / .7) * h_general + sqrt(.7 - .5^2 / .7) * z
          stats::dnorm(z) * product_weight(h_overt, 2L)
        }, -10, 10, rel.tol = 1e-12, abs.tol = 1e-13)
        largest_inner_error <<- max(largest_inner_error, inner[["abs.error"]])
        stats::dnorm(value) * product_weight(h_general, 1L) * inner[["value"]]
      }, numeric(1L))
    }, -10, 10, rel.tol = 1e-12, abs.tol = 1e-13)
    # Weights are at most one: omitted mass outside either independent
    # standard-normal coordinate is bounded by 4*pnorm(-10). Combine that
    # bound with the outer and largest reported inner adaptive error estimates.
    reference_error <- (reference[["abs.error"]] + largest_inner_error +
      4 * stats::pnorm(-10)) / reference[["value"]]
    expect_lt(reference_error, 1e-7)
    z_lower <- c(thresholds, -rev(thresholds), -Inf)
    z_upper <- c(Inf, head(z_lower, -1L))
    omega <- c(case[["weights"]], rev(head(case[["weights"]], -1L)))
    observed <- .selection_envelope_test_call(covariance, 0, sei, omega,
      z_lower = z_lower, z_upper = z_upper)
    diagnostic <- observed[["integration_diagnostics"]][1L, ]
    combined <- diagnostic[["covariance_width"]] +
      2 * diagnostic[["quadrature_change"]] + diagnostic[["tail_bound"]]
    relative_error <- abs(expm1(log(reference[["value"]]) -
      observed[["log_normalizer"]]))
    expect_equal(diagnostic[["used_covariance_envelope"]], 1)
    expect_lte(combined, .005)
    expect_lte(relative_error, combined + reference_error)
    numerator <- mvtnorm::dmvnorm(rep(0, length(sei)), sigma = covariance,
      log = TRUE) + length(sei) * log(tail(case[["weights"]], 1L))
    expect_equal(observed[["log_density"]] + observed[["log_normalizer"]],
      numerator, tolerance = 1e-12)
  }
})


test_that("repeated row geometry and singleton children retain their Gaussian law", {

  # Rows 1 and 3 are identical but remain distinct factors in the product.
  # The second fixture also has a child factor affecting only row 2. Original
  # selection SEs differ from candidate SDs after adding estimate variation.
  for (types in list(c("a", "a", "a"), c("a", "b", "a"))) {
    data <- data.frame(vi = c(.04, .09, .04), study = 1L,
      type = types, estimate = 1:3)
    covariance <- as.matrix(metafor::vcalc(vi, cluster = study, type = type,
      obs = estimate, rho = c(.7, .5), data = data))
    diag(covariance) <- diag(covariance) + .1^2
    sei <- sqrt(data[["vi"]])
    mean <- stats::qnorm(.975) * sei
    correlation <- stats::cov2cor(covariance)
    pair_sum <- sum(asin(correlation[lower.tri(correlation)]))
    for (omega in list(c(1, 0), c(1, 2))) {
      # Centering on the original selection thresholds gives the independent
      # trivariate Gaussian sign-moment identity, with no numerical integral.
      center <- mean(omega)
      jump <- omega[[1L]] - omega[[2L]]
      reference <- log(center^3 + center * jump^2 * pair_sum / (2 * pi))
      for (sign in c(-1L, 1L)) {
        observed <- .selection_envelope_test_call(covariance, sign * mean,
          sei, omega, sign)
        diagnostic <- observed[["integration_diagnostics"]][1L, ]
        expect_equal(diagnostic[["used_covariance_envelope"]], 1)
        expect_lt(abs(observed[["log_normalizer"]] - reference), 5e-4)
        expect_lte(diagnostic[["covariance_width"]] +
          2 * diagnostic[["quadrature_change"]] + diagnostic[["tail_bound"]], .005)
      }
    }
  }
})


test_that("the early dense rule retains distant tails and hard-zero transitions", {

  skip_if_not_installed("metadat")
  data("dat.assink2016", package = "metadat", envir = environment())
  study <- dat.assink2016[dat.assink2016[["study"]] == 11L, , drop = FALSE]
  sampling <- as.matrix(metafor::vcalc(vi, cluster = study, type = deltype,
    obs = esid, rho = c(.7, .5), data = study))
  sei <- sqrt(study[["vi"]])
  cases <- list(
    list(mean = -2.2, tau = 0, omega = .001, reference = -149.865411056598),
    list(mean = stats::qnorm(.975) * sei +
      rep(c(-.35, .35), length.out = length(sei)) * sqrt(.3 * sei^2 + 1e-20),
      tau = 1e-10, omega = 0, reference = -2.59746606554819)
  )
  # Independent statmod/base-R references: the distant tail agreed at 128/256
  # within 1.14e-13; the rowwise transition agreed at 96/128/192 within 7.11e-15.
  # The latter is the largest-error case from the additional GH7 stress grid.
  for (case in cases) {
    covariance <- sampling + diag(case[["tau"]]^2, length(sei))
    observed <- .selection_envelope_test_call(covariance, case[["mean"]], sei,
      c(1, case[["omega"]]), orders = c(7L, SELNORM_CLUSTER_QUADRATURE_ORDERS))
    diagnostic <- observed[["integration_diagnostics"]][1L, ]
    expect_equal(diagnostic[["used_covariance_envelope"]], 1)
    expect_lt(abs(observed[["log_normalizer"]] - case[["reference"]]), 5e-4)
    expect_lte(diagnostic[["covariance_width"]] +
      2 * diagnostic[["quadrature_change"]] + diagnostic[["tail_bound"]], .005)
    if (case[["omega"]] > 0) {
      numerator <- mvtnorm::dmvnorm(rep(0, length(sei)),
        rep_len(case[["mean"]], length(sei)), covariance, log = TRUE) +
        length(sei) * log(case[["omega"]])
      expect_equal(observed[["log_density"]] + observed[["log_normalizer"]],
        numerator, tolerance = 1e-12)
    }
  }
})
