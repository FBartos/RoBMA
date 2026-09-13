.postfit_selection_context <- function(side, rule, S = 2L) {

  prior <- BayesTools::prior_weightfunction(
    side, c(.05, .3), BayesTools::wf_fixed(c(1, 2, .4)),
    model = BayesTools::selection_model(
      other_random_effects = "integrate", known_sampling_variance = "integrate",
      weight_rule = rule, group = "paper"
    )
  )
  context <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0),
                             c(.8, 1.2), "positive")
  context$omega <- matrix(context$fixed_omega[1L, ], S,
                          length(context$p_cuts) - 1L, byrow = TRUE)
  context$kernel_mode <- c(SELKERNEL_STEP, rep(SELKERNEL_NORMAL, S - 1L))
  context$vector_rule <- rep(if (rule == "product") 0L else
    if (side == "two-sided") 2L else 1L, S)
  context$alpha <- numeric(S)
  context$phack_kind <- integer(S)
  context$use_normal <- context$kernel_mode == SELKERNEL_NORMAL
  context
}


test_that("retained estimate CDF keeps deletion available and names unsupported marginal mixtures", {

  object <- bselmodel(yi = c(.2, .5), sei = c(.2, .3), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE,
    selection = BayesTools::selection_model(estimate_random_effects = "condition"))
  expected <- matrix(c(.3, .7), 1L)
  testthat::local_mocked_bindings(.cdf_lik_estimate.brma = function(object) expected,
    .package = "RoBMA")
  expect_equal(unname(.cdf.brma(object, conditioning_depth = "estimate")), expected)
  error <- tryCatch(.cdf.brma(object, conditioning_depth = "marginal"), error = identity)
  expect_identical(conditionMessage(error), paste0(
    "Joint-selection CDF evaluation is unavailable at this conditioning depth. ",
    "Use the estimate-deletion CDF for LOO-PIT or 'as_zplot()' for marginal selected projections."
  ))
  expect_null(conditionCall(error))
})


test_that("fixed-effect mv selection has no estimate variance source", {

  dat <- data.frame(yi = c(.2, .5), vi = c(.04, .09))
  samples <- matrix(c(.1, -.2), 2L, dimnames = list(NULL, "mu"))
  fixed <- matrix(samples[, "mu"], 2L, 2L)
  for (mode in c("integrate", "condition")) {
    object <- bselmodel.mv(yi = yi, V = diag(dat$vi), random = NULL, data = dat,
      measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
      selection = BayesTools::selection_model(estimate_random_effects = mode))
    parts <- .predict_joint_selection_gaussian_parts(object, object$data, samples,
      fixed_mu = fixed, within = matrix(.7, 2L, 2L), between = matrix(.8, 2L, 2L),
      fitted_context = TRUE)
    expect_false(.selection_retains_estimate(object$data))
    expect_false(.selection_integrates_estimate(object$data))
    expect_identical(parts$latent_means, fixed)
    expect_identical(parts$random_covariance, array(0, c(2L, 2L, 2L)))
    for (s in 1:2) expect_equal(parts$covariance[s, , ], diag(dat$vi), tolerance = 0)
  }
})


test_that("partial selected vectors keep the full original publication event", {

  plan <- list(points_per_scramble = 256L, max_points_per_scramble = 4096L,
               scrambles = 8L, seed = 197L, relative_tolerance = .005,
               designs = list())
  means <- matrix(c(.1, -.2, .3, .4), 2L, 2L)
  covariance <- matrix(c(.64, .24, .24, 1.44), 2L)
  lower <- matrix(covariance[lower.tri(covariance, diag = TRUE)], 2L, 3L,
                   byrow = TRUE)
  sei <- c(.8, 1.2)
  observed <- .6
  for (side in c("one-sided", "two-sided")) {
    for (rule in c("product", "best")) {
      selection <- .postfit_selection_context(side, rule)
      result <- .selection_joint_event_numerator(
        means, lower, sei, selection, plan, 1L, observed
      )
      shared <- .selection_joint_event_numerator(
        means, lower[1L, , drop = FALSE], sei, selection, plan, 1L, observed
      )
      expect_equal(shared, result, tolerance = 1e-12)
      expect_equal(.selection_joint_condition_event_context(
        selection, matrix(observed, 1L), sei[1L]
      ), .selection_joint_condition_event_context(
        selection, matrix(observed, 2L), sei[1L]
      ), tolerance = 0)
      # Independent scalar Gaussian integral after conditioning Y2 on Y1.
      # Select the bin using the actual p-value, not any native bin code.
      p_value <- function(y, se) if (side == "two-sided") {
        2 * stats::pnorm(-abs(y / se))
      } else stats::pnorm(y / se, lower.tail = FALSE)
      weight <- function(p) c(1, 2, .4)[findInterval(p, c(0, .05, .3, 1),
                                                     rightmost.closed = TRUE)]
      conditional_mean <- means[1, 2] + covariance[1, 2] / covariance[1, 1] *
        (observed - means[1, 1])
      conditional_sd <- sqrt(covariance[2, 2] - covariance[1, 2]^2 / covariance[1, 1])
      cuts <- sort(unique(c(-Inf, Inf, sei[2] * stats::qnorm(
        c(.05, .3) / if (side == "two-sided") 2 else 1, lower.tail = FALSE),
        if (side == "two-sided") -sei[2] * stats::qnorm(c(.05, .3) / 2,
                                                       lower.tail = FALSE))))
      expected_mass <- sum(vapply(seq_len(length(cuts) - 1L), function(i) {
        stats::integrate(function(y) {
          p1 <- p_value(observed, sei[1])
          p2 <- p_value(y, sei[2])
          w <- if (rule == "product") weight(p1) * weight(p2) else weight(pmin(p1, p2))
          stats::dnorm(y, conditional_mean, conditional_sd) * w
        }, cuts[i], cuts[i + 1L], rel.tol = 1e-10)$value
      }, numeric(1L)))
      expected <- stats::dnorm(observed, means[1, 1], sqrt(covariance[1, 1]),
                               log = TRUE) + log(expected_mass)
      expect_equal(result$log_numerator[1], expected, tolerance = 1e-10)
      expect_equal(result$log_numerator[2], stats::dnorm(
        observed, means[2, 1], sqrt(covariance[1, 1]), log = TRUE), tolerance = 1e-12)
      expect_identical(result$relative_mcse, c(0, 0))
      # The same group can span independent Gaussian subblocks; retaining an
      # outcome still changes best selection for the other outcome.
      mapped <- .selection_joint_deleted_row_context(
        selection, c(observed, -1), sei, list(1:2), 2L
      )
      grid <- c(-3, -.2, 1.5, 3)
      got <- .selection_joint_log_weight(matrix(grid, 1L), rep(sei[2], 4L),
                                         BayesTools::selection_context_subset_rows(mapped, 1L))
      expected_weight <- if (rule == "best") weight(pmin(p_value(observed, sei[1]),
                                                         p_value(grid, sei[2]))) else
        weight(p_value(grid, sei[2]))
      expect_equal(got, sum(log(expected_weight)), tolerance = 1e-12)
    }
  }
})


test_that("declared rank-one partial events condition the scalar latent exactly", {

  plan <- list(points_per_scramble = 16L, scrambles = 2L, seed = 3L)
  mean <- matrix(c(.1, .3), 1L)
  loading <- matrix(c(.8, -1.2), 1L)
  sei <- c(.8, 1.2)
  observed <- .6
  missing <- .3 - 1.2 * (observed - .1) / .8
  for (rule in c("product", "best")) {
    context <- .postfit_selection_context("one-sided", rule, 1L)
    value <- .selection_joint_event_numerator(
      mean, selection_se = sei, selection_context = context,
      execution_plan = plan, observed = 1L, values = observed,
      rank_one_loading = loading
    )
    p <- stats::pnorm(c(observed, missing) / sei, lower.tail = FALSE)
    weight <- c(1, 2, .4)[findInterval(p, c(0, .05, .3, 1))]
    expected <- stats::dnorm(observed, .1, .8, log = TRUE) +
      if (rule == "product") sum(log(weight)) else log(weight[which.min(p)])
    expect_equal(value$log_numerator, expected, tolerance = 1e-12)
    expect_identical(value$relative_mcse, 0)
    excluded <- .selection_joint_event_numerator(
      mean, selection_se = sei, selection_context = context,
      execution_plan = plan, observed = 1L, values = observed,
      upper = missing - .1, rank_one_loading = loading
    )
    expect_identical(excluded$log_numerator, -Inf)
  }
  error <- tryCatch(.selection_joint_event_numerator(
    mean, selection_se = sei, selection_context = context,
    execution_plan = plan, observed = 1:2, values = c(observed, missing),
    rank_one_loading = loading
  ), error = identity)
  expect_identical(conditionMessage(error), paste0(
    "Selected partial-vector density is unavailable because the observed ",
    "rank-one law has no Lebesgue density."
  ))
  expect_null(conditionCall(error))
})


test_that("best row deletion retains physical cutoff equality in mixed branches", {

  se <- seq(.001, 2, length.out = 10000L)[4L]
  context <- .postfit_selection_context("one-sided", "best", 2L)
  context$vector_rule <- c(1L, 0L)
  context$kernel_mode <- rep(SELKERNEL_STEP, 2L)
  y <- c(se * stats::qnorm(.05, lower.tail = FALSE), 0)
  context$obs_bin <- .selection_obs_bin(y, c(se, se), context$p_cuts, 1)
  mapped <- .selection_joint_deleted_row_context(context, y, c(se, se), list(1:2), 1L)
  expect_identical(mapped$obs_bin_by_sample, c(1L, context$obs_bin[1L]))
  expect_identical(.selection_joint_log_weight(matrix(y, 1L), c(se, se),
    BayesTools::selection_context_subset_rows(context, 1L)), 0)
})


test_that("integrated true effects use the Gaussian posterior conditional on retained context", {

  # The two nonzero retained context contributions belong to the fixed mean;
  # only U is updated. This identity applies to all four conditioning cells.
  Q <- matrix(c(.5, .1, .1, .3), 2L)
  R <- matrix(c(.8, .2, .2, .6), 2L)
  fixed <- c(.2, -.1)
  h_random <- c(.7, -.4)
  h_sampling <- c(-.3, .2)
  y <- c(.4, .9)
  for (retain_random in c(FALSE, TRUE)) for (retain_sampling in c(FALSE, TRUE)) {
    random_context <- tcrossprod(h_random)
    sampling_context <- tcrossprod(h_sampling)
    random_covariance <- Q + if (retain_random) 0 else random_context
    sampling_covariance <- R + if (retain_sampling) 0 else sampling_context
    latent_mean <- fixed + if (retain_random) h_random else 0
    mean <- latent_mean + if (retain_sampling) h_sampling else 0
    parts <- list(
      means = matrix(mean, 1L), latent_means = matrix(latent_mean, 1L),
      random_covariance = array(random_covariance, c(1L, 2L, 2L)),
      sampling_covariance = array(sampling_covariance, c(1L, 2L, 2L)),
      covariance = array(random_covariance + sampling_covariance, c(1L, 2L, 2L))
    )
    expected <- latent_mean + as.vector(random_covariance %*%
      solve(random_covariance + sampling_covariance, y - mean))
    set.seed(147)
    before <- .Random.seed
    expect_equal(as.vector(.predict_joint_selection_source_posterior(parts, y, FALSE)),
                 expected, tolerance = 1e-12)
    expect_identical(.Random.seed, before)
  }
  singleton <- list(means = matrix(.2), latent_means = matrix(.3),
    random_covariance = array(.4, c(1L, 1L, 1L)),
    sampling_covariance = array(.6, c(1L, 1L, 1L)),
    covariance = array(1, c(1L, 1L, 1L)))
  expect_identical(dim(.predict_joint_selection_source_posterior(singleton, .5, FALSE)), c(1L, 1L))
  expect_equal(as.numeric(.predict_joint_selection_source_posterior(singleton, .5, FALSE)), .42)
  # A declared singular sampling law still has a deterministic true effect
  # when no true-effect source is integrated.
  singular <- list(means = matrix(c(.2, .4), 1L),
    latent_means = matrix(c(.1, .3), 1L),
    random_covariance = array(0, c(1L, 2L, 2L)),
    sampling_covariance = array(c(1, -2, -2, 4), c(1L, 2L, 2L)),
    covariance = array(c(1, -2, -2, 4), c(1L, 2L, 2L)))
  set.seed(126)
  before <- .Random.seed
  expect_identical(.predict_joint_selection_source_posterior(
    singular, c(.5, -.2), TRUE), singular$latent_means)
  expect_identical(.Random.seed, before)
})


test_that("declared singleton source reconstruction preserves Gaussian draws and boundaries", {

  make_parts <- function(q, d, blocks = lapply(seq_len(ncol(q)), identity)) {

    S <- nrow(q)
    K <- ncol(q)
    Q <- D <- array(0, c(S, K, K))
    for (k in seq_len(K)) {
      Q[, k, k] <- q[, k]
      D[, k, k] <- d[, k]
    }
    list(means = matrix(.2, S, K), latent_means = matrix(.3, S, K),
         random_covariance = Q, sampling_covariance = D, covariance = Q + D,
         dependency_blocks = blocks)
  }
  original_sampler <- .outcome_rng.norm_known_v_covariance
  sampler_calls <- 0L
  testthat::local_mocked_bindings(
    .outcome_rng.norm_known_v_covariance = function(...) {
      sampler_calls <<- sampler_calls + 1L
      original_sampler(...)
    }, .package = "RoBMA"
  )
  q <- rbind(c(.4, .7), c(0, 0), c(.3, .8))
  d <- rbind(c(.6, 1.3), c(0, 0), c(0, 0))
  parts <- make_parts(q, d)
  y <- rbind(c(.5, -.1), c(.2, .2), c(-.2, .8))
  set.seed(146)
  fast <- .predict_joint_selection_source_posterior(parts, y)
  fast_rng <- .Random.seed
  expect_identical(sampler_calls, 0L)
  reference <- parts
  reference$dependency_blocks <- NULL
  set.seed(146)
  slow <- .predict_joint_selection_source_posterior(reference, y)
  expect_identical(sampler_calls, 2L)
  expect_equal(fast, slow, tolerance = 1e-12)
  expect_identical(.Random.seed, fast_rng)
  expect_identical(as.vector(fast[2L, ]), c(.3, .3))

  # Independent scalar conditional law: mean .42, variance .24. These
  # coefficients retain the two original Gaussian phases, each of length one.
  scalar <- make_parts(matrix(.4), matrix(.6))
  set.seed(147)
  z <- stats::rnorm(2L)
  scalar_rng <- .Random.seed
  expected <- .42 + .6 * sqrt(.4) * z[1L] - .4 * sqrt(.6) * z[2L]
  set.seed(147)
  actual <- .predict_joint_selection_source_posterior(scalar, .5)
  expect_identical(dim(actual), c(1L, 1L))
  expect_equal(as.numeric(actual), expected, tolerance = 1e-12)
  expect_identical(.Random.seed, scalar_rng)
  expect_equal(as.numeric(.predict_joint_selection_source_posterior(
    scalar, .5, FALSE)), .42, tolerance = 1e-12)
  expect_identical(.Random.seed, scalar_rng)

  # Structural nonsingletons stay generic even when an evaluated matrix is diagonal.
  sampler_calls <- 0L
  nonsingleton <- make_parts(matrix(c(.4, .7), 1L), matrix(c(.6, 1.3), 1L), list(1:2))
  .predict_joint_selection_source_posterior(nonsingleton, c(.5, -.1))
  expect_identical(sampler_calls, 2L)
  # Mixed semidefinite sources retain the full-matrix spectral policy, including
  # the tiny positive diagonal beside an exact zero.
  mixed <- make_parts(matrix(c(1, 1e-22, 0), 1L), matrix(1, 1L, 3L))
  sampler_calls <- 0L
  set.seed(148)
  mixed_draw <- .predict_joint_selection_source_posterior(mixed, rep(.4, 3L))
  mixed_rng <- .Random.seed
  expect_identical(sampler_calls, 2L)
  mixed$dependency_blocks <- NULL
  set.seed(148)
  expect_identical(.predict_joint_selection_source_posterior(mixed, rep(.4, 3L)), mixed_draw)
  expect_identical(.Random.seed, mixed_rng)

  invalid <- scalar
  invalid$covariance[] <- 0
  error <- tryCatch(.predict_joint_selection_source_posterior(invalid, .5, FALSE),
                    error = identity)
  expect_identical(conditionMessage(error),
    "Selected latent posterior covariance must be positive definite.")
  expect_null(conditionCall(error))
  invalid <- nonsingleton
  invalid$dependency_blocks <- list(1L)
  error <- tryCatch(.predict_joint_selection_source_posterior(invalid, c(.5, -.1)),
                    error = identity)
  expect_identical(conditionMessage(error), "Known-V block metadata must partition the fitted rows.")
  expect_null(conditionCall(error))
})


test_that("diagonal source reconstruction preserves the full conditional Gaussian law", {

  S <- 5000L
  q <- c(.4, .2, .6)
  V <- matrix(c(.8, .15, .1, .15, .7, .2, .1, .2, .9), 3L)
  Q <- diag(q)
  mean <- c(.2, -.1, .3)
  latent <- c(.1, -.2, .4)
  y <- c(.5, .3, -.2)
  parts <- list(
    means = matrix(mean, S, 3L, byrow = TRUE),
    latent_means = matrix(latent, S, 3L, byrow = TRUE),
    random_covariance = array(rep(Q, each = S), c(S, 3L, 3L)),
    sampling_covariance = array(rep(V, each = S), c(S, 3L, 3L)),
    covariance = array(rep(Q + V, each = S), c(S, 3L, 3L)),
    random_diagonal = matrix(q, S, 3L, byrow = TRUE),
    sampling_covariance_matrix = V,
    dependency_blocks = list(1:3)
  )
  expected_mean <- latent + as.vector(Q %*% solve(Q + V, y - mean))
  expected_covariance <- Q - Q %*% solve(Q + V, Q)
  conditional <- .predict_joint_selection_source_posterior(parts, y, FALSE)
  expect_equal(conditional, matrix(expected_mean, S, 3L, byrow = TRUE), tolerance = 1e-12)
  withr::local_seed(362)
  actual <- .predict_joint_selection_source_posterior(parts, y)
  rng <- .Random.seed
  general <- parts
  general$random_diagonal <- general$sampling_covariance_matrix <- NULL
  set.seed(362)
  reference <- .predict_joint_selection_source_posterior(general, y)
  expect_equal(actual, reference, tolerance = 1e-12)
  expect_identical(.Random.seed, rng)
  expect_true(all(abs(colMeans(actual) - expected_mean) <
                    6 * sqrt(diag(expected_covariance) / S)))
  expect_equal(stats::cov(actual), expected_covariance, tolerance = .025)

  # A fixed sampling factor is identical to its repeated covariance law;
  # singular rank-one support retains the existing fallback representation.
  loading <- c(.3, -.4, .2)
  singular <- tcrossprod(loading)
  set.seed(363)
  fixed <- .outcome_rng.norm_known_v_covariance(matrix(0, 5L, 3L), singular)
  set.seed(363)
  repeated <- .outcome_rng.norm_known_v_covariance(matrix(0, 5L, 3L),
    array(rep(singular, each = 5L), c(5L, 3L, 3L)))
  expect_equal(fixed, repeated, tolerance = 1e-15)
  expect_equal(fixed[, 1L] / loading[1L], fixed[, 2L] / loading[2L], tolerance = 1e-14)
})


test_that("same-design marginal selection retains certified integrated covariance", {

  dat <- data.frame(yi = c(.2, -.1, .4), study = c("a", "a", "b"), esid = 1:3,
                    sei = sqrt(c(.09, .12, .1)))
  V <- matrix(c(.09, .02, 0, .02, .12, 0, 0, 0, .1), 3L)
  object <- bselmodel.mv(yi = yi, V = V, random = ~ 1 | study/esid, data = dat,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    prior_heterogeneity = BayesTools::prior_random(
      sd = BayesTools::prior("point", list(location = .4))),
    selection = BayesTools::selection_model(group = "study"))
  object$fit <- structure(list(), formula_design = object$formula_design,
    prior_list = c(object$formula_design$mu$prior_list,
                   .create_fit_priors(object$data, object$priors)))
  samples <- matrix(c(.1, -.2), 2L, dimnames = list(NULL, "mu_intercept"))
  fixed <- matrix(samples[, 1L], 2L, 3L)
  original_factors <- .selection_joint_random_factor_samples
  original_covariance <- .brma_mv_random_effects_marginal_vcov
  use_factors <- TRUE
  factor_calls <- 0L
  level_policy <- NULL
  testthat::local_mocked_bindings(
    .selection_joint_random_factor_samples = function(...) {
      factor_calls <<- factor_calls + 1L
      if (use_factors) original_factors(...) else NULL
    },
    .brma_mv_random_effects_marginal_vcov = function(...) {
      level_policy <<- list(...)$new_levels
      original_covariance(...)
    }, .package = "RoBMA")
  evaluate <- function(data = object$data, sampling = NULL) {

    .predict_joint_selection_gaussian_parts(object, data, samples, fixed,
      matrix(0, 2L, 3L), matrix(0, 2L, 3L), known_V_new = sampling,
      fitted_context = FALSE, draw_context = TRUE)
  }
  withr::local_seed(371)
  fast <- evaluate()
  fast_seed <- .Random.seed
  expect_equal(fast$random_diagonal, matrix(.4^2, 2L, 3L), tolerance = 1e-14)
  use_factors <- FALSE
  set.seed(371)
  reference <- evaluate()
  expect_identical(.Random.seed, fast_seed)
  expect_null(reference$random_diagonal)
  expect_identical(level_policy, "sample")
  for (field in c("means", "latent_means")) {
    expect_equal(fast[[field]], reference[[field]], tolerance = 1e-12)
  }
  # The dense BayesTools backend labels its covariance axes; the diagonal
  # calculation uses internal positional arrays. Compare their values and
  # dimensions while checking the backend labels and dependency order explicitly.
  for (field in c("random_covariance", "covariance")) {
    expect_identical(dim(fast[[field]]), c(2L, 3L, 3L))
    expect_identical(dim(reference[[field]]), dim(fast[[field]]))
    expect_identical(dimnames(reference[[field]]),
      list(draw = NULL, row = rownames(dat), column = rownames(dat)))
    expect_equal(as.numeric(fast[[field]]), as.numeric(reference[[field]]),
                 tolerance = 1e-12)
  }
  expect_identical(fast$dependency_blocks, reference$dependency_blocks)
  # New estimate identities can repeat, producing non-diagonal integrated
  # covariance. A fitted diagonal certificate must not be reused for them.
  new_rows <- data.frame(study = rep("new", 3L), esid = c(10L, 10L, 11L),
    sei = dat$sei)
  new_data <- .prepare_newdata(object, new_rows, type = "estimate",
    bias_adjusted = FALSE, include_scale = TRUE, include_random = TRUE)
  new_V <- .known_v_newdata_prepare(diag(dat$sei^2), k = 3L)
  use_factors <- TRUE
  before <- factor_calls
  changed <- evaluate(new_data, new_V)
  expect_identical(factor_calls, before)
  expect_identical(level_policy, "sample")
  expect_null(changed$random_diagonal)
  expected <- .4^2 * matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 1), 3L)
  expect_identical(dim(changed$random_covariance), c(2L, 3L, 3L))
  expect_identical(dimnames(changed$random_covariance),
    list(draw = NULL, row = rownames(new_rows), column = rownames(new_rows)))
  for (draw in 1:2) {
    expect_equal(as.numeric(changed$random_covariance[draw, , ]), as.numeric(expected),
                 tolerance = 1e-12)
  }
})


test_that("selected prediction chunks retain posterior rows and fitted-context means", {

  dat <- data.frame(yi = c(.8, -.1), vi = c(.04, .09), paper = c("a", "a"))
  object <- bselmodel(
    yi = yi, vi = vi, cluster = paper, data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .05, BayesTools::wf_fixed(c(1, .4)),
      model = BayesTools::selection_model(group = "paper", weight_rule = "best")
    )
  )
  samples <- cbind(mu = seq(.1, .2, length.out = 7L), tau = .4, rho = .3,
                    "gamma[1]" = seq(-.3, .4, length.out = 7L))
  rownames(samples) <- paste0("draw", 1:7)
  K <- nrow(dat)
  withr::local_options(RoBMA.known_v_covariance_max_bytes =
    4 * .known_v_covariance_peak_bytes(2L, K))
  original_parts <- .predict_joint_selection_gaussian_parts
  sizes <- integer()
  testthat::local_mocked_bindings(
    .conditional_effect_parameters = function(object) "mu",
    .predict_joint_selection_gaussian_parts = function(...) {
      result <- original_parts(...)
      sizes <<- c(sizes, nrow(result$means))
      expect_null(result$context_covariance)
      result
    }, .package = "RoBMA"
  )
  set.seed(127)
  before <- .Random.seed
  fitted <- blup(object, .posterior_samples = samples)
  expect_identical(.Random.seed, before)
  context_mean <- samples[, "mu"] + .4 * sqrt(.3) * samples[, "gamma[1]"]
  expected <- matrix(context_mean, 7L, K)
  expected <- expected + .4^2 * (1 - .3) /
    matrix(.4^2 * (1 - .3) + dat$vi, 7L, K, byrow = TRUE) *
    (matrix(dat$yi, 7L, K, byrow = TRUE) - expected)
  expect_equal(as.numeric(fitted), as.numeric(expected), tolerance = 1e-12)
  expect_identical(rownames(fitted), rownames(samples))
  expect_identical(sizes, c(2L, 2L, 2L, 1L))
  set.seed(128)
  first <- predict(object, type = "estimate", quiet = TRUE,
                   .posterior_samples = samples)
  set.seed(128)
  second <- predict(object, type = "estimate", quiet = TRUE,
                    .posterior_samples = samples)
  expect_identical(as.numeric(first), as.numeric(second))
  expect_identical(dim(first), c(7L, 2L))
  expect_true(all(is.finite(first)))
  expect_equal(attr(first, "RoBMA_target"), attr(second, "RoBMA_target"),
               tolerance = 0)
})


test_that("deleting part of a publication retains its original best-weight event", {

  dat <- data.frame(yi = c(.4, -.2, .3), sei = c(.8, 1.1, .6),
                    paper = rep("a", 3L))
  weights <- c(1, 2, .4)
  object <- bselmodel(
    yi = yi, sei = sei, data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", c(.05, .3), BayesTools::wf_fixed(weights),
      model = BayesTools::selection_model(group = "paper", weight_rule = "best")
    )
  )
  samples <- matrix(c(.2, .4), 1L, dimnames = list(NULL, c("mu", "tau")))
  setup <- .log_lik_posterior_setup(
    object$fit, samples, object$data, object$priors, "estimate", NULL,
    condition_local_effects = FALSE
  )
  actual <- .selection_joint_deletion_loglik_from_setup(setup, list(1:2))[, 1L]
  # Y1 and Y2 are independent before selection, conditional on retained Y3.
  # Sum their exact Gaussian rectangle probabilities over the p-value bins.
  sd <- sqrt(.4^2 + dat$sei^2)
  p <- stats::pnorm(dat$yi / dat$sei, lower.tail = FALSE)
  observed_bin <- findInterval(p, c(0, .05, .3, 1))
  boundary <- stats::qnorm(c(0, .05, .3, 1), lower.tail = FALSE)
  bin_mass <- lapply(1:2, function(row) {
    -diff(stats::pnorm(dat$sei[row] * boundary, .2, sd[row]))
  })
  normalizer <- sum(outer(1:3, 1:3, Vectorize(function(i, j) {
    bin_mass[[1]][i] * bin_mass[[2]][j] * weights[min(i, j, observed_bin[3])]
  })))
  expected <- sum(stats::dnorm(dat$yi[1:2], .2, sd[1:2], log = TRUE)) +
    log(weights[min(observed_bin)]) - log(normalizer)
  expect_equal(actual, expected, tolerance = 1e-12)
})


test_that("bridge setup integrates known V without requesting a retained decomposition", {

  dat <- data.frame(yi = c(.1, -.2), paper = c("a", "a"))
  V <- matrix(c(.04, .01, .01, .09), 2L)
  object <- bselmodel.mv(
    yi = yi, V = V, random = ~1 | paper, data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .05, BayesTools::wf_fixed(c(1, .4)),
      model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "integrate",
        group = "paper"
      )
    )
  )
  fit_priors <- .create_fit_priors(object$data, object$priors)
  result <- .marglik_sampling_latent_setup(object$data, object$priors, fit_priors)
  expect_true(result$marginalized)
  expect_identical(result$fit_priors, fit_priors)
  expect_identical(result$diagnostics$included, "known V")
  metadata <- .selection_postfit_target_metadata(object$data)
  expect_identical(metadata$sampling_structure$policy,
                   "full_sampling_covariance_integration")
  expect_identical(metadata$sampling_structure$source, "whole_sampling_error")
  expect_identical(metadata$sampling_structure$covariance, "V")
  expect_identical(metadata$selection_model$known_sampling_variance, "integrate")
  expect_identical(metadata$selection_model$groups$requested, "paper")
  expect_length(metadata$selection_model$groups$row_blocks, 1L)
})


test_that("product deletion batches preserve the conditional row law and log tails", {

  dat <- data.frame(yi = c(-.6, .2), sei = c(.8, 1.2), paper = c("a", "a"))
  object <- bselmodel(
    yi = yi, sei = sei, data = dat, measure = "GEN", only_priors = TRUE,
    prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", c(.05, .3), BayesTools::wf_fixed(c(1, 2, .4)),
      model = BayesTools::selection_model(other_random_effects = "integrate",
        known_sampling_variance = "integrate", group = "paper")
    )
  )
  setup <- list(data = object$data, S = 3L, K = 2L, selection_sei = dat$sei)
  y <- -dat$yi
  base_mean <- matrix(c(.1, .3, -.4, -.2, .4, .2), 3L, 2L)
  # Independent two-variable conditional identity for Sigma=(.64,.24;.24,1.44).
  mean <- cbind(base_mean[, 1L] + .24 / 1.44 * (y[2L] - base_mean[, 2L]),
                base_mean[, 2L] + .24 / .64 * (y[1L] - base_mean[, 1L]))
  variance <- matrix(c(.64 - .24^2 / 1.44, 1.44 - .24^2 / .64), 3L, 2L,
                     byrow = TRUE)
  context <- .postfit_selection_context("one-sided", "product", 3L)
  context$kernel_mode <- c(SELKERNEL_STEP, SELKERNEL_NORMAL, SELKERNEL_STEP)
  context$use_normal <- context$kernel_mode == SELKERNEL_NORMAL
  context$omega[3L, 2L] <- 0
  context$obs_bin <- .selection_obs_bin(y, dat$sei, context$p_cuts, 1)
  conditional <- list(y = y, means = mean, variance = variance, lower_tail = FALSE,
    selection_context = context, groups = list(1:2), dependency_blocks = list(1:2))
  testthat::local_mocked_bindings(
    .selection_joint_conditional_summary_from_setup = function(setup) conditional,
    .package = "RoBMA"
  )
  got <- .selection_joint_estimate_targets(setup)
  expected <- lapply(got, function(x) matrix(NA_real_, 3L, 2L))
  for (s in 1:3) for (i in 1:2) {
    heights <- if (s == 2L) c(1, 1, 1) else if (s == 3L) c(1, 0, .4) else c(1, 2, .4)
    weight <- function(x) heights[findInterval(
      stats::pnorm(x / dat$sei[i], lower.tail = FALSE), c(0, .05, .3, 1),
      rightmost.closed = TRUE
    )]
    cuts <- sort(c(-Inf, dat$sei[i] * stats::qnorm(c(.05, .3), lower.tail = FALSE), Inf))
    mass <- function(lower = -Inf, upper = Inf, power = 0L) {
      sum(vapply(seq_len(length(cuts) - 1L), function(j) {
        a <- max(lower, cuts[j])
        b <- min(upper, cuts[j + 1L])
        if (a >= b) return(0)
        stats::integrate(function(x) x^power * weight(x) *
          stats::dnorm(x, mean[s, i], sqrt(variance[s, i])), a, b,
          rel.tol = 1e-10, abs.tol = 1e-12)$value
      }, numeric(1L)))
    }
    normalizer <- mass()
    first <- mass(power = 1L) / normalizer
    expected$log_density[s, i] <- stats::dnorm(y[i], mean[s, i], sqrt(variance[s, i]),
      log = TRUE) + log(weight(y[i])) - log(normalizer)
    expected$cdf[s, i] <- mass(lower = y[i]) / normalizer
    expected$log_lower[s, i] <- log(expected$cdf[s, i])
    expected$log_upper[s, i] <- log(mass(upper = y[i]) / normalizer)
    expected$mean[s, i] <- -first
    expected$variance[s, i] <- mass(power = 2L) / normalizer - first^2
  }
  attr(expected, "dependency_blocks") <- list(1:2)
  expect_equal(got, expected, tolerance = 1e-9)
  expect_identical(got$log_density[3L, 1L], -Inf)
  for (component in names(got)) {
    expect_identical(.selection_joint_estimate_targets(setup, component)[[component]],
                     got[[component]])
  }
  conditional$selection_context$omega[1L, ] <- 0
  error <- tryCatch(.selection_joint_estimate_targets(setup, "log_density"), error = identity)
  expect_identical(conditionMessage(error), paste0(
    "Selected row deletion is unavailable because its conditional selection event ",
    "cannot be normalized."
  ))
  expect_null(conditionCall(error))
})


test_that("selection conditionals reuse declared covariance roots and reject invalid states", {

  dat <- data.frame(yi = c(.3, -.2, .1), paper = rep("a", 3L), time = 1:3)
  D <- c(.09, .16, .25)
  B <- matrix(c(.1, -.2, .15), 3L)
  object <- bselmodel.mv(
    yi = yi, V = known_v_factor(D, B), random = ~ ar1(time | paper),
    data = dat, measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    prior_heterogeneity = BayesTools::prior_random(
      sd = BayesTools::prior("point", list(location = .4)),
      covariance = BayesTools::random_covariance(
        cor = BayesTools::prior("point", list(location = .3)), cor_scale = "cor"
      )
    ),
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .05, BayesTools::wf_fixed(c(1, .4)),
      model = BayesTools::selection_model(other_random_effects = "integrate",
        known_sampling_variance = "integrate", group = "paper")
    )
  )
  fit <- structure(list(), formula_design = object$formula_design,
    prior_list = c(object$formula_design$mu$prior_list,
                   .create_fit_priors(object$data, object$priors)))
  samples <- matrix(c(.1, -.2), 2L, dimnames = list(NULL, "mu_intercept"))
  setup <- .log_lik_evaluated_setup(fit, object$data, object$priors,
    "estimate", NULL, matrix(samples[, 1L], 2L, 3L), matrix(0, 2L, 3L), NULL,
    samples, random_effects_conditioning = "included_in_mu")
  fixed_mean <- matrix(samples[, 1L], 2L, 3L)
  reference <- function(rho) {
    sigma <- diag(D) + tcrossprod(B) + .4^2 * rho^abs(outer(1:3, 1:3, "-"))
    mean <- variance <- matrix(NA_real_, 2L, 3L)
    for (i in 1:3) {
      other <- setdiff(1:3, i)
      coefficient <- solve(sigma[other, other], sigma[other, i])
      mean[, i] <- fixed_mean[, i] + as.vector(
        (matrix(dat$yi[other], 2L, 2L, byrow = TRUE) - fixed_mean[, other]) %*% coefficient)
      variance[, i] <- sigma[i, i] - sum(sigma[i, other] * coefficient)
    }
    list(mean = mean, variance = variance)
  }
  original_plan <- .brma_mv_random_effects_marginal_factor_plan
  original_summary <- .marglik_covariance_plan_conditional_summary_batch
  mode <- "interior"
  route <- NULL
  testthat::local_mocked_bindings(
    .brma_mv_random_effects_marginal_factor_plan = function(...) {
      value <- original_plan(...)
      if (mode == "boundary") {
        # The compiler's public fixed-correlation support is unchanged. Test
        # the declared internal root contract using its analytic endpoint.
        root <- matrix(0, 3L, 3L)
        root[, 1L] <- .4
        value$factor_states[[2L]][[1L]]$coefficient_factor <- root
        value$factor_states[[2L]][[1L]]$markov_transition <- c(1, 1)
        value$factor_states[[2L]][[1L]]$markov_innovation_variance <- c(0, 0)
      } else if (mode == "unsupported") {
        value$factor_plans[[1L]]$coefficient_structure <- "unsupported"
      } else if (mode == "invalid") {
        value$factor_states[[1L]][[1L]]$markov_innovation_variance[1L] <- -1
      }
      value
    },
    .marglik_covariance_plan_conditional_summary_batch = function(...) {
      route <<- list(...)$random_covariance_plans[[1L]]$coefficient_structure
      original_summary(...)
    }, .package = "RoBMA"
  )
  expected <- reference(.3)
  for (current in c("interior", "boundary", "unsupported")) {
    mode <- current
    route <- NULL
    actual <- .selection_joint_conditional_summary_from_setup(setup)
    target <- expected
    if (mode == "boundary") {
      boundary <- reference(1)
      target$mean[2L, ] <- boundary$mean[2L, ]
      target$variance[2L, ] <- boundary$variance[2L, ]
    }
    expect_equal(actual$means, target$mean, tolerance = 1e-11)
    expect_equal(actual$variance, target$variance, tolerance = 1e-11)
    expect_identical(route, switch(mode, interior = "markov", boundary = "dense", NULL))
  }
  mode <- "invalid"
  error <- tryCatch(.selection_joint_conditional_summary_from_setup(setup), error = identity)
  expect_identical(conditionMessage(error), "Bridge-marginalized random-effect Markov state is invalid.")
  expect_null(conditionCall(error))
  state <- list(coefficient_scale = rep(.4, 3L), markov_transition = c(1, 1),
                markov_innovation_variance = c(0, 0))
  expect_error(.marglik_validate_random_covariance_markov_state(state, 3L),
               "Bridge-marginalized random-effect Markov state is invalid.", fixed = TRUE)
  expect_identical(.marglik_validate_random_covariance_markov_state(state, 3L, TRUE), state)
  state$markov_transition[2L] <- NaN
  error <- tryCatch(.marglik_validate_random_covariance_markov_state(state, 3L, TRUE), error = identity)
  expect_identical(conditionMessage(error), "Bridge-marginalized random-effect Markov state is invalid.")
  expect_null(conditionCall(error))
})
