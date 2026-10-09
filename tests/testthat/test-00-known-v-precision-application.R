# The precision application B^-1 (r - U c) of a tight known covariance and the
# GLS residual variance M_ii - x_i' G^-1 x_i subtract nearly equal numbers. The
# references below are 80-digit mpmath solutions of the dense systems.

.precision_application <- function(y, sampling_covariance, plans, states,
                                   block_indices = list(seq_along(y))) {

  .marglik_covariance_plan_precision_residual_batch(
    cache                    = new.env(parent = emptyenv()),
    y                        = y,
    means                    = matrix(0, 1L, length(y)),
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = plans,
    random_covariance_states = list(states),
    block_indices            = block_indices,
    extra_variances          = matrix(0, 1L, length(y))
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
