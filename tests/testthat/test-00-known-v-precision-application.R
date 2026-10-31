# The precision application B^-1 (r - U c) of a tight known covariance and the
# GLS residual variance M_ii - x_i' G^-1 x_i subtract nearly equal numbers. The
# references below are 80-digit mpmath solutions of the dense systems.

.precision_application <- function(y, sampling_covariance, plans, states,
                                   block_indices = list(seq_along(y)),
                                   extra_variances = rep(0, length(y))) {

  .marglik_covariance_plan_precision_residual_batch(
    cache                    = new.env(parent = emptyenv()),
    y                        = y,
    means                    = matrix(0, 1L, length(y)),
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = plans,
    random_covariance_states = list(states),
    block_indices            = block_indices,
    extra_variances          = matrix(extra_variances, 1L)
  )
}

.group_plan <- function(model_matrix, group_map = rep(1L, nrow(model_matrix))) {

  list(type = "group", model_matrix = model_matrix, group_map = group_map,
       coefficient_structure = "diagonal")
}

.identity_state <- function(rank) {

  list(coefficient_factor = diag(1, nrow = rank, ncol = rank))
}


test_that("the rank-one precision application keeps the constant direction", {

  # M = d I + 1 1': the precision applied to (1, 1) is 1 / (d + 2) in every
  # entry; subtracting U c from r left 2 at d = 1e-10 (0.5 expected).
  references <- list(
    "1e-6" = list(
      d = 1e-6,
      "1,1" = c(0.49999975000012499994, 0.49999975000012499994),
      "1,0" = c(500000.24999987500006, -499999.75000012499994),
      "gen" = c(999999.65000017499991, -1000000.3499998250001)),
    "1e-10" = list(
      d = 1e-10,
      "1,1" = c(0.499999999975, 0.499999999975),
      "1,0" = c(5000000000.25, -4999999999.75),
      "gen" = c(9999999999.65, -10000000000.35)),
    "1e-16" = list(
      d = 1e-16,
      "1,1" = c(0.499999999999999975, 0.499999999999999975),
      "1,0" = c(5000000000000000.25, -4999999999999999.75),
      "gen" = c(9999999999999999.65, -10000000000000000.35)),
    "2^-60" = list(
      d = 2^-60,
      "1,1" = c(0.49999999999999999978, 0.49999999999999999978),
      "1,0" = c(576460752303423488.25, -576460752303423487.75),
      "gen" = c(1152921504606846975.6, -1152921504606846976.3))
  )
  vectors <- list("1,1" = c(1, 1), "1,0" = c(1, 0), "gen" = c(.3, -1.7))
  for (reference in references) {
    for (name in names(vectors)) {
      observed <- .precision_application(
        vectors[[name]], diag(reference$d, 2L),
        list(.group_plan(matrix(1, 2L, 1L))), list(.identity_state(1L))
      )
      expect_equal(as.numeric(observed), reference[[name]], tolerance = 1e-14,
        info = paste("d =", reference$d, "y = (", name, ")"))
    }
  }
})


test_that("block-base precision applications keep tight sampling blocks", {

  # Two sampling blocks V = d [[1, .5], [.5, 1]] and one group effect shared by
  # all four rows, M = V + 1 1'.
  references <- list(
    list(d = 1e-6, ones = rep(0.24999990625003515624, 4L),
      general = c(1599999.975000009375, -2400000.024999990625,
                  1199999.975000009375, -400000.024999990625)),
    list(d = 1e-10, ones = rep(0.249999999990625, 4L),
      general = c(15999999999.975, -24000000000.025, 11999999999.975,
                  -4000000000.025)),
    list(d = 1e-16, ones = rep(0.24999999999999999062, 4L),
      general = c(15999999999999999.975, -24000000000000000.025,
                  11999999999999999.975, -4000000000000000.025)),
    list(d = 2^-60, ones = rep(0.24999999999999999992, 4L),
      general = c(1844674407370955161.6, -2767011611056432742.4,
                  1383505805528216371.2, -461168601842738790.42))
  )
  for (reference in references) {
    sampling <- matrix(0, 4L, 4L)
    sampling[1:2, 1:2] <- sampling[3:4, 3:4] <- reference$d *
      matrix(c(1, .5, .5, 1), 2L)
    for (case in c("ones", "general")) {
      y <- if (case == "ones") rep(1, 4L) else c(.3, -1.7, .9, .1)
      observed <- .precision_application(
        y, sampling, list(.group_plan(matrix(1, 4L, 1L))),
        list(.identity_state(1L)), block_indices = list(1:4)
      )
      expect_equal(as.numeric(observed), reference[[case]], tolerance = 1e-13,
        info = paste("d =", reference$d, case))
    }
  }
})


test_that("sparse latent systems of grouped effects keep tight sampling variances", {

  # Three groups of two rows share a group effect, M = d I + blockdiag(1 1'); the
  # latent system K = I + U' U / d is diagonal and stored sparsely.
  references <- list(
    list(d = 1e-6, value = c(999999.65000017499991, -1000000.3499998250001,
      400000.24999987500006, -399999.75000012499994, 800000.19999990000005,
      -799999.80000009999995)),
    list(d = 1e-10, value = c(9999999999.65, -10000000000.35, 4000000000.25,
      -3999999999.75, 8000000000.2, -7999999999.8)),
    list(d = 1e-16, value = c(9999999999999999.65, -10000000000000000.35,
      4000000000000000.25, -3999999999999999.75, 8000000000000000.2,
      -7999999999999999.8)),
    list(d = 2^-60, value = c(1152921504606846975.6, -1152921504606846976.3,
      461168601842738790.65, -461168601842738790.15, 922337203685477581.0,
      -922337203685477580.6))
  )
  y <- c(.3, -1.7, .9, .1, 1.2, -.4)
  plan <- .group_plan(matrix(1, 6L, 1L), c(1L, 1L, 2L, 2L, 3L, 3L))
  for (reference in references) {
    observed <- .precision_application(
      y, diag(reference$d, 6L), list(plan), list(.identity_state(1L))
    )
    expect_equal(as.numeric(observed), reference$value, tolerance = 1e-14,
      info = paste("d =", reference$d))
  }
})


test_that("the condition gate does not depend on the scale of the covariance", {

  # The sums of absolute entries that estimate the 1-norm overflow for entries
  # near the largest double although the matrix is well conditioned (a
  # reciprocal condition number of 0.12): the gate refused these as
  # ill-conditioned. 900-digit mpmath solutions of (V + diag(e)) x = y.
  sampling <- 6.3e307 * matrix(c(1, .8, .7, .8, 1, .75, .7, .75, 1), 3L)
  extra <- c(1e307, 2e307, 3e307)
  references <- list(
    list(y = c(1, 2, 3), value = c(-1.6708740750427293e-308,
      1.5994128945279828e-308, 3.2055192198165288e-308)),
    list(y = c(1e300, 2e300, 3e300), value = c(-1.6708740750427294e-8,
      1.5994128945279829e-8, 3.205519219816529e-8))
  )
  for (reference in references) {
    observed <- .precision_application(reference$y, sampling, list(), list(),
      extra_variances = extra)
    expect_equal(as.numeric(observed), reference$value, tolerance = 1e-12)
  }

  # The latent system K = I + U'U of two columns U = 2^511 T, T = [1 1; 1 1; 1 0],
  # has a condition number of 10 and column sums above the largest double (its
  # entries are 1 + 3 2^1022 and 2^1023); M = I + U U' is not conditioned
  # for the dense route. 900-digit mpmath solution of M x = (1, 2, 3).
  model_matrix <- 2^511 * rbind(c(1, 1), c(1, 1), c(1, 0))
  observed <- .precision_application(c(1, 2, 3), diag(3L),
    list(.group_plan(model_matrix)), list(.identity_state(2L)))
  reference <- c(-.5, .5, 1.0012832363282406e-307)
  expect_equal(as.numeric(observed)[1:2], reference[1:2], tolerance = 1e-12)
  # The third component results from a cancellation of unit-size terms and is
  # zero to double precision next to the others. A relative tolerance on the
  # whole vector tolerates an error of 1e-12 there, so it is compared with its
  # own absolute tolerance: the application returns 1e-30, whereas a loss of
  # precision in the scaled latent system leaves errors of 1e-16 and above.
  expect_lt(abs(observed[[3L]] - reference[[3L]]), 1e-20)
})


test_that("the condition gate of 2 x 2 blocks uses the exact reciprocal condition number", {

  # Blocks of order one and two have the 1-norm reciprocal condition number in
  # closed form, det(M) / ||M||_1^2 (order one is the latent system of a single
  # coefficient, tested above, whose value is 1). The gate refuses a block at or
  # below eps / 1e-6, and reports the reciprocal condition number of the
  # refused block, which is the exact value up to the rounding of its Cholesky
  # factor (what the LAPACK estimate it replaced returned). The covariance is
  # V + diag(extra) with unequal extra variances, which keeps the block on the
  # dense route (equal ones are handled by the eigensystem of V).
  threshold <- .Machine$double.eps / 1e-6
  reciprocal_condition <- function(sampling, extra) {

    outcome <- tryCatch(
      .precision_application(c(1, 1), sampling, list(), list(),
        extra_variances = extra),
      RoBMA_known_v_precision_unavailable = function(condition) condition
    )
    if (inherits(outcome, "condition")) outcome$reciprocal_condition else NA_real_
  }

  # M = s diag(1, 1 / r): the reciprocal condition number is 1 / r exactly, for
  # eigenvalue ratios r from 1 to 1e18 (10^9.64 and 10^9.67 straddle the gate,
  # 1 / threshold = 10^9.65) and scales from 1e-12 to 1e12.
  for (scale in c(1e-12, 1, 1e12)) {
    for (exponent in c(0, 3, 6, 9, 9.64, 9.67, 10, 12, 14, 16, 18)) {
      reciprocal <- 10^-exponent
      observed <- reciprocal_condition(
        scale * diag(c(1, reciprocal / 2), 2L), scale * c(0, reciprocal / 2)
      )
      info <- paste("scale", scale, "ratio 1e", exponent)
      if (reciprocal > threshold) {
        expect_true(is.na(observed), info = info)
      } else {
        expect_equal(observed, reciprocal, tolerance = 1e-12, info = info)
      }
    }
  }

  # M = s [[b^2 + delta, b], [b, 1]] and its permutation with b = 2^12: the
  # entries are integers times a power of two, so det(M) = delta s^2 and
  # ||M||_1 = (b^2 + delta + b) s are exact and the reciprocal condition number
  # is delta / (b^2 + delta + b)^2. The Cholesky factor cancels
  # l22^2 = 1 - b^2 / (b^2 + delta), a relative error of at most
  # eps b^2 / delta in the reciprocal condition number. V is M less 2^-30 on
  # the second diagonal entry (positive definite: delta > b^2 2^-30).
  b <- 2^12
  shift <- 2^-30
  for (delta in c(1, 7, 100, 5000, 60000, 65000)) {
    exact <- delta / (b^2 + delta + b)^2
    for (scale in c(2^-40, 1, 2^40)) {
      for (swap in c(FALSE, TRUE)) {
        sampling <- matrix(c(b^2 + delta, b, b, 1 - shift), 2L)
        extra <- c(0, shift)
        if (swap) {
          sampling <- sampling[2:1, 2:1]
          extra <- extra[2:1]
        }
        observed <- reciprocal_condition(scale * sampling, scale * extra)
        info <- paste("delta", delta, "scale", scale, "swapped", swap)
        if (exact > threshold) {
          expect_true(is.na(observed), info = info)
        } else {
          expect_equal(observed, exact,
            tolerance = 16 * .Machine$double.eps * b^2 / delta + 1e-12,
            info = info)
        }
      }
    }
  }
})


test_that("block-base precision applications refuse ill-conditioned sampling blocks", {

  # Two sampling blocks [[1, rho], [rho, 1]] and one group effect shared by the
  # four rows, M = V + 1 1'. The base blocks have the reciprocal condition
  # numbers about (1 - rho) / 2; the application was returned with no check on
  # them, with norm-wise errors of 5.6e-8 (1 - rho = 1e-9), 5.6e-5 (1e-12) and
  # 5.6e-4 (1e-13). References: 900-digit mpmath solutions of M x = y.
  apply_blocks <- function(gap, y) {

    rho <- 1 - gap
    block <- matrix(c(1, rho, rho, 1), 2L)
    sampling <- matrix(0, 4L, 4L)
    sampling[1:2, 1:2] <- sampling[3:4, 3:4] <- block
    .precision_application(y, sampling, list(.group_plan(matrix(1, 4L, 1L))),
      list(.identity_state(1L)), block_indices = list(1:4))
  }
  expected <- list(
    "1e-09" = list(ones = rep(0.16666666669444444, 4L),
      split = c(0.33333333347222222, 0.33333333347222222,
                -0.16666666677777777, -0.16666666677777777)))
  for (case in c("ones", "split")) {
    y <- if (case == "ones") rep(1, 4L) else c(1, 1, 0, 0)
    observed <- apply_blocks(1e-9, y)
    reference <- expected[["1e-09"]][[case]]
    expect_lte(max(abs(observed - reference)) / max(abs(reference)), 1e-6)
  }
  for (gap in c(1e-11, 1e-12, 1e-13)) {
    for (y in list(rep(1, 4L), c(1, 1, 0, 0))) {
      expect_error(apply_blocks(gap, y),
        class = "RoBMA_known_v_precision_unavailable")
    }
  }
})


test_that("ill-conditioned precision applications are accurate or unavailable", {

  # One row at 1e-18, the others at `others`, and the uncentred loadings
  # (1, x) with x in [2000, 2020]: M = B + U U' is nearly singular. A returned
  # application is accurate to the condition gate (condition number times the
  # machine epsilon at most 1e-6; 9.5e-10 is observed at the condition number
  # 2.9e9 of the unit-scale case); otherwise it is a classed unavailability,
  # never a rounded number (the subtraction returned errors up to 73).
  x <- 2000 + 20 * (0:5) / 5
  y <- c(.5, -1.25, .75, 2, -.3, 1)
  references <- list(
    "1" = c(0.31016088706462401499, -1.7510099775524440739,
            0.24798004489511185193, 1.4969700673426677777,
            -0.80403991020977629646, 0.49495011223777962935),
    "1e-14" = c(99536405778175.357205, -179629124624827.99182,
                15731797109766.198819, 136092718844360.38946,
                -98546359421045.419897, 26814562313548.770744),
    "1e-10" = c(9954542248.3930292189, -17963636085.041374241,
                1572727730.3718290335, 13609091545.785032308,
                -9854544638.8017644168, 2681819176.611438858),
    "1e-6" = c(992352.9856412144184, -1796158.3709877856852,
               157683.25802343627654, 1361524.8870346582383,
               -984633.48395411979991, 269208.14505710216187)
  )
  for (name in names(references)) {
    outcome <- tryCatch(
      .precision_application(
        y, diag(c(1e-18, rep(as.numeric(name), 5L))),
        list(.group_plan(cbind(1, x))), list(.identity_state(2L))
      ),
      RoBMA_known_v_precision_unavailable = function(condition) condition
    )
    if (inherits(outcome, "condition")) {
      expect_s3_class(outcome, "RoBMA_numerical_unavailable")
      expect_s3_class(outcome, "error")
      expect_true(all(outcome$reciprocal_condition < 1e-9))
      expect_match(conditionMessage(outcome), "unavailable", fixed = TRUE)
      # the positions are those among the draws the covariance was applied to
      expect_match(conditionMessage(outcome),
        "at draw 1 among the posterior draws it was applied to", fixed = TRUE)
      expect_false(identical(name, "1"), info = "the unit-scale case is applied")
    } else {
      scale <- max(abs(references[[name]]))
      expect_lte(max(abs(as.numeric(outcome) - references[[name]])) / scale,
        1e-6)
      expect_identical(name, "1")
    }
  }
})


test_that("the tight GLS residual variance has no cancellation", {

  # M = d I + 1 1' with an intercept: the residual variance is
  # M_ii - 1 / (1' M^-1 1) = d / 2, which the subtraction returned as 3.8e-6
  # at d = 1e-10.
  for (d in c(1e-6, 1e-10, 1e-16, 2^-60)) {
    covariance_diagonal <- matrix(1 + d, 1L, 2L)
    plan_data <- list(
      sampling_covariance      = diag(d, 2L),
      random_covariance_plans  = list(.group_plan(matrix(1, 2L, 1L))),
      random_covariance_states = list(list(.identity_state(1L))),
      block_indices            = list(1:2),
      extra_variances          = matrix(0, 1L, 2L),
      covariance_diagonal      = covariance_diagonal
    )
    testthat::local_mocked_bindings(
      .known_v_marginal_factor_plan = function(...) plan_data,
      .package = "RoBMA"
    )
    X <- matrix(1, 2L, 1L)
    attr(X, "rank") <- 1L
    projection <- .known_v_factor_gls_projection_batch(
      object = list(), posterior_samples = matrix(0, 1L, 1L), known_V = list(),
      X = X, y = c(0, sqrt(d)), return_se = TRUE, return_resid = TRUE
    )
    expect_equal(as.numeric(projection$residual_variance), rep(d / 2, 2L),
      tolerance = 1e-12, info = paste("d =", d))
    expect_equal(as.numeric(projection$H_diag), c(.5, .5), tolerance = 1e-12)
  }
})


test_that("the tight GLS residual variance forms no K x K identity per draw", {

  # An indicator column for row 1 (a singleton level) makes row 1 a tight row
  # in every draw. The directions e_1 - H_1 are K values; a K x K identity per
  # draw (32 MB at K = 2000) must not be allocated for them.
  K <- 2000L
  S <- 3L
  sampling <- diag(K)
  plan_data <- list(
    sampling_covariance      = sampling,
    random_covariance_plans  = list(),
    random_covariance_states = rep(list(list()), S),
    block_indices            = list(seq_len(K)),
    extra_variances          = matrix(0, S, K),
    covariance_diagonal      = matrix(1, S, K)
  )
  recorded <- list()
  testthat::local_mocked_bindings(
    .known_v_marginal_factor_plan = function(...) plan_data,
    .known_v_factor_quadratic = function(plan_data, draw, vectors) {
      recorded[[draw]] <<- vectors
      rep(0, ncol(vectors))
    },
    .package = "RoBMA"
  )
  X <- cbind(1, c(1, rep(0, K - 1L)))
  attr(X, "rank") <- 2L
  before <- gc(reset = TRUE)
  projection <- .known_v_factor_gls_projection_batch(
    object = list(), posterior_samples = matrix(0, S, 1L), known_V = list(),
    X = X, y = rep(0, K), return_se = TRUE
  )
  growth <- (gc()["Vcells", "max used"] - before["Vcells", "used"]) * 8 / 1024^2
  expect_lt(growth, 8)

  # With W = I the hat column of row 1 is X (X'X)^-1 x_1 = e_1, so the
  # direction e_1 - H_1 vanishes up to rounding.
  hat_column <- as.numeric(X %*% solve(crossprod(X), X[1L, ]))
  expected <- -hat_column
  expected[1L] <- expected[1L] + 1
  expect_length(recorded, S)
  for (draw in seq_len(S)) {
    expect_identical(dim(recorded[[draw]]), c(K, 1L))
    expect_lt(max(abs(as.numeric(recorded[[draw]]) - expected)), 1e-12)
  }
  expect_equal(projection$H_diag[1L, 1L], 1, tolerance = 1e-12)
})


test_that("the factor-form residual variance matches the dense projection", {

  # Group and known-group effects, a row-scaled factor and a dense sampling
  # block with one tight row: w' M w in factor form equals the dense quadratic.
  set.seed(20261009)
  K <- 6L
  sampling <- crossprod(matrix(rnorm(K * K), K)) / K + diag(.1, K)
  basis <- matrix(rnorm(K * 2L), K)
  plans <- list(
    .group_plan(basis, c(1L, 1L, 2L, 2L, 3L, 3L)),
    list(type = "known_group", model_matrix = matrix(1, K, 1L),
         group_map = c(1L, 2L, 3L, 1L, 2L, 3L),
         group_covariance = crossprod(matrix(rnorm(9), 3L)) / 3,
         coefficient_structure = "diagonal"),
    list(type = "row_group", model_matrix = matrix(1, K, 1L),
         group_map = rep(1L, K), coefficient_structure = "diagonal")
  )
  states <- list(
    list(coefficient_factor = matrix(c(.7, .2, 0, .5), 2L)),
    list(coefficient_factor = matrix(.4, 1L, 1L)),
    list(coefficient_factor = matrix(.9, 1L, 1L), row_scale = seq(.5, 1.5, length.out = K))
  )
  plan_data <- list(
    sampling_covariance      = sampling,
    random_covariance_plans  = plans,
    random_covariance_states = list(states),
    extra_variances          = matrix(seq(.01, .06, length.out = K), 1L)
  )
  random <- matrix(0, K, K)
  for (index in seq_along(plans)) {
    plan  <- plans[[index]]
    state <- states[[index]]
    factor <- plan$model_matrix %*% state$coefficient_factor
    if (identical(plan$type, "row_group")) factor <- factor * state$row_scale
    mask <- if (identical(plan$type, "known_group")) {
      plan$group_covariance[plan$group_map, plan$group_map]
    } else {
      outer(plan$group_map, plan$group_map, "==")
    }
    random <- random + tcrossprod(factor) * mask
  }
  covariance <- sampling + diag(plan_data$extra_variances[1L, ]) + random
  vectors <- matrix(rnorm(K * 3L), K)
  expect_equal(
    .known_v_factor_quadratic(plan_data, 1L, vectors),
    colSums(vectors * (covariance %*% vectors)),
    tolerance = 1e-12
  )
})


test_that("rstandard keeps the tight GLS residual variance of a known-V factor", {

  # y = (0, sqrt(d)) with the factor covariance M = d I + 1 1' and an intercept:
  # the residuals are -/+ sqrt(d) / 2 with the variance d / 2, so the
  # standardized residuals are -/+ 1 / sqrt(2). The subtraction returned the
  # standard error 0.00195 at d = 1e-10 (7.07e-6 expected). The residuals
  # themselves come from the intercept X'M^-1 y, whose sum cancels by a factor
  # 1 / d, so z is accurate to about 1e-16 / d.
  old_max_cores <- RoBMA.get_option("max_cores")
  on.exit(RoBMA.options(max_cores = old_max_cores), add = TRUE)
  RoBMA.options(max_cores = 1L)
  for (d in c(1e-6, 1e-10)) {
    object <- brma.mv(
      yi = c(0, sqrt(d)), V = known_v_factor(rep(d, 2L), matrix(1, 2L, 1L)),
      measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE
    )
    object$fit <- coda::mcmc(matrix(c(0, 0), 2L, dimnames = list(NULL, "mu")))
    standardized <- rstandard(object)
    expect_equal(unname(standardized[, "se"]), rep(sqrt(d / 2), 2L),
      tolerance = 1e-12, info = paste("d =", d))
    expect_equal(unname(standardized[, "z"]), c(-1, 1) / sqrt(2),
      tolerance = 1e-16 / d * 100, info = paste("d =", d))
  }
})
