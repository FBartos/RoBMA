test_that("Gaussian z densities retain scaled mixtures and extreme tail mass", {

  z    <- c(-Inf, -10, -1, 0, 1, 10, Inf)
  sei  <- c(.1, 1, 10)
  mean <- matrix(c(-1, 0, 1, .3, .8, -.5), nrow = 2)
  sd   <- matrix(c(.05, .2, 1, 2, 8, 20), nrow = 2)
  reference <- vapply(z, function(value) {

    rowMeans(stats::dnorm(
      value,
      mean = sweep(mean, 2L, sei, "/"),
      sd   = sweep(sd, 2L, sei, "/")
    ))
  }, numeric(2))
  expect_equal(
    .zplot_normal_density_matrix(z, mean, sd, sei),
    reference,
    tolerance = 1e-14
  )

  # A tiny z-scale SD recovers a finite density even when dnorm(score) would
  # underflow or become subnormal before multiplication by the Jacobian.
  scores <- c(38, 40)
  expected <- exp(stats::dnorm(scores, log = TRUE) - log(1e-308))
  observed <- .zplot_normal_density_matrix(
    scores * 1e-308, matrix(0), matrix(1e-8), 1e300
  )
  expect_equal(as.numeric(observed), expected, tolerance = 1e-12)
  expect_equal(
    as.numeric(.zplot_normal_density_matrix(
      1, matrix(-1e308), matrix(1e308), 1e308
    )),
    stats::dnorm(2),
    tolerance = 1e-14
  )
})


test_that("zplot total SD matches replicated-SE evaluation", {

  tau_within <- matrix(
    c(0, .1, 1, 10, .25, .5, 2, 20, .75, 1.5, 3, 30),
    nrow = 4,
    ncol = 3
  )
  sei     <- c(.05, .2, 2)
  sei_mat <- matrix(sei, nrow = nrow(tau_within), ncol = ncol(tau_within),
                    byrow = TRUE)
  expected <- .root_sum_squares(tau_within, sei_mat)

  expect_identical(.zplot_total_sd(tau_within, sei), expected)
  expect_equal(
    .zplot_total_sd(tau_within, sei),
    sqrt(sweep(tau_within^2, 2, sei^2, "+")),
    tolerance = 1e-15
  )
})


test_that("joint z marginals match an independently integrated bivariate law", {

  cutoff <- stats::qnorm(.975)
  sei <- c(.7, 1.3)
  prior <- BayesTools::prior_weightfunction(
    "one-sided", .025, BayesTools::wf_fixed(c(1, .2)),
    model = BayesTools::selection_model(
      other_random_effects = "integrate", known_sampling_variance = "integrate", group = "paper"
    )
  )
  spec <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0), sei,
                          effect_direction = "positive", signed_data = FALSE)
  spec$omega       <- matrix(c(1, .2, 1, .2), 2, 2, byrow = TRUE)
  spec$alpha       <- c(0, 0)
  spec$phack_kind  <- c(0L, 0L)
  spec$kernel_mode <- c(SELKERNEL_STEP, SELKERNEL_NORMAL)
  spec$vector_rule <- rep(0L, nrow(spec$omega))
  spec$use_normal  <- c(FALSE, TRUE)
  rho <- .7
  covariance <- matrix(rep(c(sei[1]^2, prod(sei) * rho, sei[2]^2), each = 2), 2)
  weight <- function(z) ifelse(z > cutoff, 1, .2)
  other_weight <- function(z) {
    .2 + .8 * stats::pnorm(cutoff, rho * z, sqrt(1 - rho^2), lower.tail = FALSE)
  }
  numerator <- function(z) stats::dnorm(z) * weight(z) * other_weight(z)
  normalizer <- stats::integrate(numerator, -Inf, cutoff, rel.tol = 1e-10)$value +
    stats::integrate(numerator, cutoff, Inf, rel.tol = 1e-10)$value
  z <- c(-3, -.5, 0, 1, 2.5, 4)
  control <- set_selection_likelihood_control(points_per_scramble = 256,
    max_points_per_scramble = 16384, relative_tolerance = .002)
  result <- .zplot_joint_block(z, matrix(0, 2, 2),
    covariance, sei, spec, FALSE, control)
  expect_equal(result$density[1, ], numerator(z) / normalizer, tolerance = .002)
  expect_equal(result$density[2, ], stats::dnorm(z), tolerance = 1e-12)
  expect_equal(exp(-result$log_density), c(normalizer, 1), tolerance = .002)
  factor <- .zplot_joint_block(z, matrix(0, 2, 2),
    covariance, sei, spec, FALSE, control,
    factors = list(residual_sd = matrix(rep(sqrt(1 - rho) * sei, each = 2), 2),
                   loading = matrix(rep(sqrt(rho) * sei, each = 2), 2)))
  expect_equal(factor$density[1, ], numerator(z) / normalizer, tolerance = .002)
  probability <- .zplot_joint_block(cutoff, matrix(0, 2, 2),
    covariance, sei, spec, TRUE, control)
  expected_p <- (stats::integrate(numerator, -Inf, -cutoff)$value +
    stats::integrate(numerator, cutoff, Inf)$value) / normalizer
  expect_equal(probability$density[, 1], c(expected_p, .05), tolerance = .002)
  whole_line <- .zplot_joint_block(0, matrix(.3, 2, 2),
    covariance, sei, spec, TRUE, control)
  expect_equal(as.numeric(whole_line$density), c(1, 1), tolerance = 1e-12)
  # The reciprocal remains a useful normalizer check, independent of the
  # normalized pre-selection reference used by the public vector display.
  expect_gt(1 / normalizer - 1, 0)
  spec$sign <- -1L
  reflected <- .zplot_joint_block(-z, matrix(0, 2, 2),
    covariance, sei, spec, FALSE, control)
  expect_equal(reflected$density[1, ], numerator(z) / normalizer, tolerance = .002)

  object <- bselmodel.mv(
    yi                  = c(0, 0),
    data                = data.frame(paper = c("p1", "p1")),
    V                   = known_v_factor((1 - rho) * sei^2,
      matrix(sqrt(rho) * sei, ncol = 1)),
    prior_bias          = prior,
    measure             = "GEN",
    prior_unit_information_sd = 1,
    only_priors         = TRUE,
    silent              = TRUE
  )
  expect_identical(
    .data_selection_execution_plan(object$data)$block_methods, "rank_one"
  )
  testthat::local_mocked_bindings(
    selection_qmc_design = function(...) stop("Rank-one metadata must use quadrature."),
    .package = "BayesTools"
  )
  spec$sign <- 1L
  routed <- .zplot_joint_marginal(
    object, matrix(0, 2, 1),
    list(mu = matrix(0, 2, 2), mu_extrapolated = matrix(0, 2, 2), sei = sei),
    spec, z, FALSE, control
  )
  expect_equal(routed$fitted[1, ], numerator(z) / normalizer, tolerance = .002)
  expect_equal(routed$fitted[2, ], stats::dnorm(z), tolerance = 1e-12)
  expect_identical(routed$weights, c(1, 1))
  expect_equal(routed$extrapolated, matrix(rep(stats::dnorm(z), each = 2), 2), tolerance = 1e-12)
  extrapolated <- .zplot_joint_marginal(
    object, matrix(0, 2, 1),
    list(mu = matrix(0, 2, 2), mu_extrapolated = matrix(0, 2, 2), sei = sei),
    spec, cutoff, TRUE, control, extrapolate_only = TRUE
  )
  expect_null(extrapolated$fitted)
  expect_identical(extrapolated$weights, c(1, 1))
  expect_equal(extrapolated$EDR, c(.05, .05), tolerance = 1e-12)
})


test_that("exact z normalizers retain zero observed weights and normal branches", {

  cutoff <- stats::qnorm(.975)
  sei    <- c(.7, 1.3)
  rho    <- .7
  prior  <- BayesTools::prior_weightfunction(
    "one-sided", .025, BayesTools::wf_fixed(c(1, 0))
  )
  context <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0), sei,
    effect_direction = "positive", signed_data = FALSE)
  context$omega       <- matrix(c(1, 0, 1, 0), 2, 2, byrow = TRUE)
  context$alpha       <- c(0, 0)
  context$phack_kind  <- c(0L, 0L)
  context$kernel_mode <- c(SELKERNEL_STEP, SELKERNEL_NORMAL)
  context$vector_rule <- rep(0L, nrow(context$omega))
  context$use_normal  <- c(FALSE, TRUE)
  means <- matrix(rep(cutoff * sei, each = 2), 2)
  covariance <- matrix(rep(c(sei[1]^2, rho * prod(sei), sei[2]^2), each = 2), 2)
  control <- set_selection_likelihood_control(relative_tolerance = .002)
  plans <- lapply(c("rank_one", "factor", "dense"), function(method) {

    .selection_joint_execution_plan(
      row_blocks             = list(1:2),
      block_methods          = method,
      factor_ranks           = switch(method, rank_one = 1L, factor = 2L, 0L),
      selection_control      = control,
      sampling               = NULL,
      sampling_factor_blocks = NULL,
      random_covariance      = NULL
    )
  })
  # Each factorization has marginal variance one and correlation rho. The
  # selected region begins at its mean, giving the Gaussian quadrant identity.
  expected <- .25 + asin(rho) / (2 * pi)
  rank_one <- .selection_joint_cluster_loglik_block(
    c(0, 0), means, matrix(rep(sqrt(.3) * sei, each = 2), 2),
    matrix(rep(sqrt(.7) * sei, each = 2), 2), sei, context, plans[[1]],
    return_normalizer = TRUE
  )
  loading <- cbind(sqrt(.8) * sei, sqrt(.1) * sei * c(1, -1))
  factor <- .selection_joint_factor_loglik_block(
    c(0, 0), means, matrix(rep(sqrt(.1) * sei, each = 2), 2),
    matrix(rep(as.numeric(loading), each = 2), 2), sei, context, plans[[2]], 1L,
    return_normalizer = TRUE
  )
  dense <- .selection_joint_dense_loglik_block(
    c(0, 0), means, covariance, sei, context, plans[[3]], 2L,
    return_normalizer = TRUE
  )
  for (result in list(rank_one, factor, dense)) {
    expect_identical(result$log_density[[1]], -Inf)
    expect_identical(result$log_normalizer[[2]], 0)
    expect_equal(exp(result$log_normalizer), c(expected, 1), tolerance = .002)
  }
  singleton <- .selection_joint_dense_loglik_block(
    0, means[, 1, drop = FALSE], covariance[, 1, drop = FALSE], sei[1],
    BayesTools::selection_context_subset_observations(context, 1L),
    plans[[3]], 1L, return_normalizer = TRUE
  )
  expect_identical(singleton$log_density[[1]], -Inf)
  expect_equal(singleton$log_normalizer, c(log(.5), 0), tolerance = 1e-12)
})


test_that("declared factor projections match independent Gaussian rectangles", {

  cutoff <- stats::qnorm(.025, lower.tail = FALSE)
  rectangle_normalizer <- function(mean, covariance, sei) {

    K <- length(mean)
    regions <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), K)))
    sum(vapply(seq_len(nrow(regions)), function(index) {

      significant <- regions[index, ]
      sign <- ifelse(significant, -1, 1)
      .2^sum(!significant) * as.numeric(mvtnorm::pmvnorm(
        lower = rep(-Inf, K), upper = sign * cutoff * sei,
        mean = sign * mean, sigma = covariance * tcrossprod(sign),
        algorithm = if (K <= 3) mvtnorm::TVPACK(abseps = 1e-10) else
          mvtnorm::Miwa(steps = 128L)
      ))
    }, numeric(1)))
  }
  prior <- BayesTools::prior_weightfunction(
    "one-sided", .025, BayesTools::wf_fixed(c(1, .2))
  )
  z <- c(-3, -.5, 0, 1, 2.5, 4)
  # Count the actual native factor ranks without replacing the calculation.
  native <- get(".Call", envir = environment(.zplot_joint_block))
  project <- .zplot_joint_block
  environment(project) <- new.env(parent = environment(project))
  environment(project)$.zplot_joint_block <- project
  environment(project)$.Call <- function(.NAME, ..., PACKAGE = NULL) {

    arguments <- list(...)
    native_ranks <<- c(native_ranks,
      ncol(arguments[[18L]][[2L]]) %/% ncol(arguments[[2L]]))
    native(.NAME, ..., PACKAGE = PACKAGE)
  }
  # Structural supports remain authoritative for inactive loading columns.
  testthat::local_mocked_bindings(
    selection_qmc_design = function(...) stop("Declared factors must use quadrature."),
    .package = "BayesTools"
  )
  for (case in c("1", "2", "3", "4", "nonnested")) {
    rank <- if (case == "nonnested") 2L else as.integer(case)
    K    <- max(3L, rank)
    loading <- rbind(c(.3, .1, 0, 0), c(.4, .2, .1, 0),
                     c(.5, .25, .15, .1), c(.6, .35, .2, .15))
    loading <- loading[seq_len(K), seq_len(rank), drop = FALSE]
    if (case == "nonnested") loading[cbind(c(1L, 3L), c(2L, 1L))] <- 0
    mean <- c(-.2, .1, .3, .5)[seq_len(K)]
    sei  <- c(.7, .8, 1, 1.2)[seq_len(K)]
    sd   <- c(.8, .9, 1, 1.1)[seq_len(K)]
    loadings <- list(loading, matrix(0, K, rank))
    if (rank > 1L) {
      partial <- loading
      partial[, rank] <- 0
      loadings[[3L]] <- partial
    }
    S <- length(loadings)
    covariances <- lapply(loadings, function(value) diag(sd^2) + tcrossprod(value))
    normalizers <- vapply(covariances, function(covariance) {

      rectangle_normalizer(mean, covariance, sei)
    }, numeric(1))
    references <- lapply(seq_along(covariances), function(state) {

      covariance <- covariances[[state]]
      vapply(z, function(value) {

        mean(vapply(seq_len(K), function(row) {

          other <- setdiff(seq_len(K), row)
          cross <- covariance[other, row]
          conditional_mean <- mean[other] + cross / covariance[row, row] *
            (value * sei[row] - mean[row])
          conditional_covariance <- covariance[other, other, drop = FALSE] -
            tcrossprod(cross) / covariance[row, row]
          stats::dnorm(value, mean[row] / sei[row],
                       sqrt(covariance[row, row]) / sei[row]) *
            ifelse(value > cutoff, 1, .2) * rectangle_normalizer(
              conditional_mean, conditional_covariance, sei[other]
            ) / normalizers[state]
        }, numeric(1)))
      }, numeric(1))
    })
    context <- .selection_spec(list(outcome = list(bias = prior)),
      rep(0, K), sei, effect_direction = "positive", signed_data = FALSE)
    context$omega       <- matrix(rep(c(1, .2), S), S, byrow = TRUE)
    context$alpha       <- numeric(S)
    context$phack_kind  <- integer(S)
    context$kernel_mode <- rep(SELKERNEL_STEP, S)
    context$vector_rule <- rep(0L, nrow(context$omega))
    context$use_normal  <- rep(FALSE, S)
    arguments <- list(
      z = z, mean = matrix(rep(mean, each = S), S), covariance_lower = do.call(rbind,
        lapply(covariances, function(value) value[lower.tri(value, diag = TRUE)])),
      sei = sei, selection = context, probability = FALSE,
      control = set_selection_likelihood_control(relative_tolerance = 1e-6),
      factors = list(
        residual_sd = matrix(rep(sd, each = S), S),
        loading = do.call(rbind, lapply(loadings, as.numeric)),
        loading_support = loading != 0
      )
    )
    native_ranks <- integer()
    result <- do.call(project, arguments)
    expect_setequal(native_ranks, if (rank > 1L) c(rank, rank - 1L) else rank)
    expect_equal(result$density, do.call(rbind, references), tolerance = 2e-6)
    expect_equal(exp(-result$log_density), normalizers, tolerance = 2e-6)
    row_normalizer <- .2 + .8 * stats::pnorm(cutoff * sei, mean, sd,
                                            lower.tail = FALSE)
    independent <- vapply(z, function(value) {

      mean(stats::dnorm(value, mean / sei, sd / sei) *
             ifelse(value > cutoff, 1, .2) / row_normalizer)
    }, numeric(1))
    expect_equal(result$density[2, ], independent, tolerance = 1e-12)
    expect_equal(exp(-result$log_density[2]), prod(row_normalizer),
                 tolerance = 1e-12)
    # A seq with a selection cutoff has several different rounded spacings.
    # Its marginal law must agree at the independently certified points.
    arguments$z <- sort(unique(c(seq(-3, 4, by = .05), z, cutoff)))
    grid <- do.call(project, arguments)
    expect_equal(grid$density[, match(z, arguments$z)], do.call(rbind, references),
                 tolerance = 2e-6)
    arguments$z <- c(0, cutoff)
    arguments$probability <- TRUE
    probability <- do.call(project, arguments)
    expected_tail <- mean((
      .2 * stats::pnorm(-cutoff * sei, mean, sd) +
        stats::pnorm(cutoff * sei, mean, sd, lower.tail = FALSE)
    ) / row_normalizer)
    expect_equal(probability$density[, 1], rep(1, S), tolerance = 1e-12)
    expect_equal(probability$density[2, 2], expected_tail, tolerance = 1e-12)
    if (case == "3") {
      invalid <- arguments
      invalid$factors$residual_sd[3, 1] <- -1
      expect_error(do.call(project, invalid),
        "'residual_sd' must be finite and positive.", fixed = TRUE)
      invalid <- arguments
      invalid$factors$loading <- cbind(invalid$factors$loading, 0)
      expect_error(do.call(project, invalid), "'loading' has invalid dimensions.",
                   fixed = TRUE)
    }
  }
})


test_that("nested density grids retain rare selection and recovered tail mass", {

  coefficient <- .01
  variance    <- 1 + coefficient^2
  loading     <- cbind(rep(coefficient, 3), 0, 0)
  covariance  <- diag(1, 3) + tcrossprod(loading)
  z <- sort(unique(c(seq(-1, 1.5, by = .05), -1e-8, 0, 1e-8)))
  for (case in list(c(-10, 1e-30), c(-38, 1e-300), c(-10, 0))) {
    location <- case[1]
    omega    <- case[2]
    prior <- BayesTools::prior_weightfunction(
      "one-sided", .5, BayesTools::wf_fixed(c(1, omega))
    )
    context <- .selection_spec(list(outcome = list(bias = prior)),
      rep(0, 3), rep(1, 3), effect_direction = "positive", signed_data = FALSE)
    context$omega       <- matrix(c(1, omega), 1)
    context$alpha       <- 0
    context$phack_kind  <- 0L
    context$kernel_mode <- SELKERNEL_STEP
    context$vector_rule <- rep(0L, nrow(context$omega))
    context$use_normal  <- FALSE
    arguments <- list(
      z = z, mean = matrix(location, 1, 3),
      covariance_lower = matrix(covariance[lower.tri(covariance, diag = TRUE)], 1),
      sei = rep(1, 3), selection = context, probability = FALSE,
      control = set_selection_likelihood_control(relative_tolerance = 1e-7),
      factors = list(residual_sd = matrix(1, 1, 3), loading = matrix(loading, 1),
                     loading_support = matrix(TRUE, 3, 3))
    )
    # Condition the three-observation Gaussian law on one observation. Integrate
    # the other two selection weights over its conditional shared factor.
    log_scale <- max(log(omega), stats::pnorm(0, location, sqrt(variance),
                                              lower.tail = FALSE, log.p = TRUE))
    log_normalizer <- function(u) {

      if (omega == 0) {
        stats::pnorm(0, location + coefficient * u, 1, lower.tail = FALSE, log.p = TRUE)
      } else {
        log(omega + (1 - omega) *
              stats::pnorm(0, location + coefficient * u, 1, lower.tail = FALSE))
      }
    }
    area <- stats::integrate(function(u) {
      exp(stats::dnorm(u, log = TRUE) + 3 * (log_normalizer(u) - log_scale))
    }, -Inf, Inf, rel.tol = 1e-10, abs.tol = 0)$value
    reference <- vapply(z, function(value) {

      conditional_mean <- coefficient / variance * (value - location)
      conditional_sd   <- sqrt(1 / variance)
      other <- stats::integrate(function(u) {
        exp(stats::dnorm(u, log = TRUE) + 2 *
              (log_normalizer(conditional_mean + conditional_sd * u) - log_scale))
      }, -Inf, Inf, rel.tol = 1e-10, abs.tol = 0)$value
      exp(stats::dnorm(value, location, sqrt(variance), log = TRUE) +
            log(ifelse(value >= 0, 1, omega)) - log_scale + log(other / area))
    }, numeric(1))
    result <- do.call(.zplot_joint_block, arguments)
    positive <- reference > 0
    expect_equal(as.numeric(result$density)[positive] / reference[positive],
                 rep(1, sum(positive)), tolerance = 1e-7)
    expect_identical(as.numeric(result$density)[!positive], reference[!positive])
    arguments$z <- -rev(z)
    arguments$mean <- -arguments$mean
    arguments$selection$sign <- -1L
    reflected <- do.call(.zplot_joint_block, arguments)
    expect_equal(as.numeric(reflected$density), rev(as.numeric(result$density)),
                 tolerance = 1e-12)
  }
})


test_that("conditional z marginals integrate conditional selection normalizers", {

  cutoff <- stats::qnorm(.975)
  prior <- BayesTools::prior_weightfunction(
    "one-sided", .025, BayesTools::wf_fixed(c(1, .2))
  )
  spec <- .selection_spec(list(outcome = list(bias = prior)), 0, 1,
                          effect_direction = "positive", signed_data = FALSE)
  spec$omega       <- matrix(c(1, .2), 1)
  spec$alpha       <- 0
  spec$phack_kind  <- 0L
  spec$kernel_mode <- SELKERNEL_STEP
  spec$vector_rule <- rep(0L, nrow(spec$omega))
  spec$use_normal  <- FALSE
  z <- c(-3, 0, 1, 2.5, 4)
  normalizer <- function(u) .2 + .8 * stats::pnorm(cutoff, .3 + .8 * u, 1,
                                                 lower.tail = FALSE)
  reference <- vapply(z, function(value) {
    stats::integrate(function(u) {
      stats::dnorm(u) * stats::dnorm(value, .3 + .8 * u, 1) / normalizer(u)
    }, -Inf, Inf, rel.tol = 1e-10)$value
  }, numeric(1))
  result <- .zplot_latent_mixture(z, matrix(.3), matrix(1), matrix(.8), 1,
    spec, FALSE, set_selection_likelihood_control(relative_tolerance = 1e-7))
  expect_equal(as.numeric(result$extrapolated), reference, tolerance = 1e-7)
  expect_equal(as.numeric(result$fitted), reference * ifelse(z > cutoff, 1, .2),
               tolerance = 1e-7)
  threshold <- .zplot_latent_mixture(cutoff, matrix(.3), matrix(1), matrix(.8), 1,
    spec, TRUE, set_selection_likelihood_control(relative_tolerance = 1e-7))
  expected_area <- stats::integrate(function(u) stats::dnorm(u) / normalizer(u),
                                   -Inf, Inf, rel.tol = 1e-10)$value
  expect_equal(threshold$weights, expected_area, tolerance = 1e-7)
  expect_equal(result$weights, expected_area, tolerance = 1e-7)
  whole_line <- .zplot_latent_mixture(0, matrix(.3), matrix(1), matrix(.8), 1,
    spec, TRUE, set_selection_likelihood_control(relative_tolerance = 1e-7))
  expect_equal(as.numeric(whole_line$fitted), 1, tolerance = 1e-12)
  expect_equal(whole_line$EDR, 1, tolerance = 1e-12)
})


test_that("latent z densities cannot accept a missed narrow Gaussian kernel", {

  cutoff <- stats::qnorm(.025, lower.tail = FALSE)
  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .2)))
  context <- .selection_spec(list(outcome = list(bias = prior)), .100137, .1,
    effect_direction = "positive", signed_data = FALSE)
  context$omega       <- matrix(c(1, .2), 1L)
  context$alpha       <- 0
  context$phack_kind  <- 0L
  context$kernel_mode <- SELKERNEL_STEP
  context$vector_rule <- 0L
  context$use_normal  <- FALSE
  z <- sort(unique(c(seq(-6, 6, by = .05), cutoff)))
  mu <- .100137
  candidate_sd <- 1e-7
  retained_sd <- .1
  total_sd <- sqrt(candidate_sd^2 + retained_sd^2)
  # Independent adaptive integration, including the nontrivial value exactly
  # at the selection jump. The original population-node quadrature returned
  # an all-zero curve at both orders 15 and 31 and falsely accepted it.
  reference <- vapply(z, function(value) {
    y <- value * .1
    conditional_mean <- mu * candidate_sd^2 / total_sd^2 +
      y * retained_sd^2 / total_sd^2
    conditional_sd <- candidate_sd * retained_sd / total_sd
    inverse <- stats::integrate(function(u) {
      location <- conditional_mean + conditional_sd * u
      stats::dnorm(u) / (.2 + .8 * stats::pnorm(
        (location - .1 * cutoff) / candidate_sd))
    }, -Inf, Inf, rel.tol = 1e-10)$value
    .1 * stats::dnorm(y, mu, total_sd) *
      ifelse(value >= cutoff, 1, .2) * inverse
  }, numeric(1L))
  result <- .zplot_latent_mixture(z, matrix(mu), matrix(candidate_sd),
    matrix(retained_sd), .1, context, FALSE,
    set_selection_likelihood_control(), fitted_only = TRUE)
  expect_equal(as.numeric(result$fitted), reference, tolerance = 1e-6)
})


test_that("latent z densities preserve step weights, directions, and normal rows", {

  mu <- .3
  candidate_sd <- .8
  cutoffs <- stats::qnorm(c(.025, .25), lower.tail = FALSE)
  control <- set_selection_likelihood_control(relative_tolerance = 1e-7)
  for (retained_sd in c(.4, 1.1)) for (weight in list(c(1, 3, 1e-8), c(1, 0, .2))) {
    z <- if (retained_sd <= candidate_sd) {
      sort(unique(c(seq(-6, 6, by = .05), cutoffs)))
    } else c(-3, -.3, .3, 1, 3)
    prior <- BayesTools::prior_weightfunction("one-sided", c(.025, .25),
      BayesTools::wf_fixed(weight))
    context <- .selection_spec(list(outcome = list(bias = prior)), mu, 1,
      effect_direction = "positive", signed_data = FALSE)
    context$omega       <- matrix(rep(weight, 2L), 2L, byrow = TRUE)
    context$alpha       <- numeric(2L)
    context$phack_kind  <- integer(2L)
    context$kernel_mode <- c(SELKERNEL_STEP, SELKERNEL_NORMAL)
    context$vector_rule <- integer(2L)
    context$use_normal  <- c(FALSE, TRUE)
    acceptance <- function(location) {
      lower <- stats::pnorm(cutoffs[2L], location, candidate_sd)
      middle <- stats::pnorm(cutoffs[1L], location, candidate_sd) - lower
      upper <- stats::pnorm(cutoffs[1L], location, candidate_sd, lower.tail = FALSE)
      weight[1L] * upper + weight[2L] * middle + weight[3L] * lower
    }
    # Adaptive integration over the original retained source is independent
    # of the conditional-Gaussian density implementation.
    reference <- vapply(z, function(value) {
      observed <- weight[if (value >= cutoffs[1L]) 1L else
        if (value >= cutoffs[2L]) 2L else 3L]
      if (observed == 0) return(0)
      observed * stats::integrate(function(u) {
        location <- mu + retained_sd * u
        stats::dnorm(u) * stats::dnorm(value, location, candidate_sd) /
          acceptance(location)
      }, -Inf, Inf, rel.tol = 1e-10)$value
    }, numeric(1L))
    for (direction in c(1, -1)) {
      context$sign <- as.integer(direction)
      actual <- .zplot_latent_mixture(direction * z,
        matrix(direction * mu, 2L), matrix(candidate_sd, 2L),
        matrix(retained_sd, 2L), 1, context, FALSE, control, fitted_only = TRUE)
      expect_equal(as.numeric(actual$fitted[1L, ]), reference, tolerance = 1e-7)
      expect_equal(as.numeric(actual$fitted[2L, ]),
        stats::dnorm(z, mu, sqrt(candidate_sd^2 + retained_sd^2)), tolerance = 1e-12)
      expect_identical(as.numeric(actual$fitted[1L, reference == 0]),
        reference[reference == 0])
    }
  }
})


test_that("retained Gaussian projection reuses only unchanged grid geometry", {

  cutoff <- stats::qnorm(.025, lower.tail = FALSE)
  sei <- c(.1, .3, .7)
  location <- matrix(c(-.1, .2, .4, -.3, .2, .6), 2L)
  sd <- matrix(c(.4, .6, .5, .8, .9, 1.1), 2L)
  retained <- matrix(c(.25, .1, .2, .3, .5, .4), 2L)
  omega <- c(.2, 3)
  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, omega[[1L]])))
  context <- .selection_spec(list(outcome = list(bias = prior)), rep(0, 3L), sei,
    effect_direction = "positive", signed_data = FALSE)
  context$omega <- cbind(1, omega)
  context$alpha <- numeric(2L)
  context$phack_kind <- integer(2L)
  context$kernel_mode <- rep(SELKERNEL_STEP, 2L)
  context$vector_rule <- integer(2L)
  context$use_normal <- rep(FALSE, 2L)
  rule <- .gauss_hermite_nodes(15L)
  grids <- list(
    sort(unique(c(seq(-6, 6, length.out = 241L), cutoff, cutoff + 1e-6, -sqrt(2)))),
    c(1.2, -.3, cutoff, .11, -sqrt(2), 4.8)
  )
  for (z in grids) {
    # Direct finite-mixture densities use base R Gaussian functions and the
    # analytic one-cutoff acceptance, independently of the grid recurrence.
    expected <- vapply(z, function(value) vapply(seq_len(2L), function(row) {
      mean(vapply(seq_along(sei), function(column) {
        nodes <- location[row, column] + retained[row, column] * rule$nodes
        acceptance <- omega[row] + (1 - omega[row]) * stats::pnorm(
          cutoff * sei[column], nodes, sd[row, column], lower.tail = FALSE)
        sum(rule$weights * sei[column] * stats::dnorm(value * sei[column],
          nodes, sd[row, column]) / acceptance)
      }, numeric(1L))) * if (value >= cutoff) 1 else omega[row]
    }, numeric(1L)), numeric(2L))
    actual <- .zplot_selnorm_density_matrix(z, location, sd, sei, context, FALSE,
      latent_sd = retained, quadrature = rule)
    expect_equal(actual, expected, tolerance = 1e-12)
  }
})


test_that("conditional projection refinement is independent across posterior rows", {

  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .2)))
  context <- .selection_spec(list(outcome = list(bias = prior)), 0, 1,
    effect_direction = "positive", signed_data = FALSE)
  context$omega       <- matrix(rep(c(1, .2), 3), 3, byrow = TRUE)
  context$alpha       <- numeric(3)
  context$phack_kind  <- integer(3)
  context$kernel_mode <- c(SELKERNEL_STEP, SELKERNEL_STEP, SELKERNEL_NORMAL)
  context$vector_rule <- rep(0L, nrow(context$omega))
  context$use_normal  <- c(FALSE, FALSE, TRUE)
  z       <- c(-3, 0, 1, 2.5, 4)
  mean    <- matrix(.3, 3, 1)
  sd      <- matrix(1, 3, 1)
  latent  <- matrix(c(.8, 2, .8), 3, 1)
  control <- set_selection_likelihood_control(relative_tolerance = 1e-7)
  result <- .zplot_latent_mixture(z, mean, sd, latent, 1, context, FALSE, control)
  individual <- lapply(seq_len(3), function(row) {

    .zplot_latent_mixture(z, mean[row, , drop = FALSE], sd[row, , drop = FALSE],
      latent[row, , drop = FALSE], 1,
      BayesTools::selection_context_subset_rows(context, row), FALSE, control)
  })
  expect_identical(result$fitted, do.call(rbind, lapply(individual, `[[`, "fitted")))
  expect_identical(result$extrapolated,
    do.call(rbind, lapply(individual, `[[`, "extrapolated")))
  expect_identical(result$weights, unlist(lapply(individual, `[[`, "weights")))

  cutoff <- stats::qnorm(.025, lower.tail = FALSE)
  for (row in 1:2) {
    normalizer <- function(u) {

      .2 + .8 * stats::pnorm(cutoff, .3 + latent[row, 1] * u, 1, lower.tail = FALSE)
    }
    reference <- vapply(z, function(value) {

      stats::integrate(function(u) stats::dnorm(u) *
        stats::dnorm(value, .3 + latent[row, 1] * u, 1) / normalizer(u),
        -Inf, Inf, rel.tol = 1e-10)$value
    }, numeric(1))
    expect_equal(result$extrapolated[row, ], reference, tolerance = 1e-7)
    expect_equal(result$fitted[row, ], reference * ifelse(z >= cutoff, 1, .2),
      tolerance = 1e-7)
  }
  expect_equal(result$fitted[3, ], stats::dnorm(z, .3, sqrt(1 + .8^2)),
    tolerance = 1e-12)
  thresholds <- .zplot_latent_mixture(cutoff, mean, sd, latent, 1, context,
    TRUE, control)
  individual <- lapply(seq_len(3), function(row) {

    .zplot_latent_mixture(cutoff, mean[row, , drop = FALSE], sd[row, , drop = FALSE],
      latent[row, , drop = FALSE], 1,
      BayesTools::selection_context_subset_rows(context, row), TRUE, control)
  })
  expect_identical(thresholds$fitted, do.call(rbind, lapply(individual, `[[`, "fitted")))
  expect_identical(thresholds$EDR, unlist(lapply(individual, `[[`, "EDR")))
  expect_identical(thresholds$weights, unlist(lapply(individual, `[[`, "weights")))
})


test_that("factor importance integration agrees with a one-factor reference", {

  cutoff <- stats::qnorm(.975)
  prior <- BayesTools::prior_weightfunction(
    "one-sided", .025, BayesTools::wf_fixed(c(1, .05))
  )
  spec <- .selection_spec(list(outcome = list(bias = prior)), rep(0, 3), rep(1, 3),
    effect_direction = "positive", signed_data = FALSE)
  spec$omega       <- matrix(c(1, .05), 1)
  spec$alpha       <- 0
  spec$phack_kind  <- 0L
  spec$kernel_mode <- SELKERNEL_STEP
  spec$vector_rule <- rep(0L, nrow(spec$omega))
  spec$use_normal  <- FALSE
  loading <- matrix(0, 3, 3)
  loading[, 1] <- sqrt(.7)
  testthat::local_mocked_bindings(
    .gauss_hermite_nodes = function(...) stop("Missing factor supports require QMC."),
    .package = "RoBMA"
  )
  covariance <- diag(.3, 3) + tcrossprod(loading)
  normalizer <- function(u) .05 + .95 * stats::pnorm(cutoff, sqrt(.7) * u,
    sqrt(.3), lower.tail = FALSE)
  area <- stats::integrate(function(u) stats::dnorm(u) * normalizer(u)^3,
    -Inf, Inf, rel.tol = 1e-10)$value
  z <- c(-1, 0, 2.5, 4)
  reference <- vapply(z, function(value) {
    stats::integrate(function(u) stats::dnorm(u) * normalizer(u)^2 *
      stats::dnorm(value, sqrt(.7) * u, sqrt(.3)), -Inf, Inf,
      rel.tol = 1e-10)$value * ifelse(value > cutoff, 1, .05) / area
  }, numeric(1))
  result <- .zplot_joint_block(z, matrix(0, 1, 3),
    matrix(covariance[lower.tri(covariance, diag = TRUE)], 1), rep(1, 3),
    spec, FALSE, set_selection_likelihood_control(relative_tolerance = .002,
      max_points_per_scramble = 32768),
    factors = list(residual_sd = matrix(sqrt(.3), 1, 3), loading = matrix(loading, 1),
                   loading_support = NULL))
  expect_equal(as.numeric(result$density), reference, tolerance = .002)
  expect_equal(exp(-result$log_density), area, tolerance = .002)
})


test_that("zplot reports unavailable targets and failed integration explicitly", {

  control <- set_selection_likelihood_control(points_per_scramble = 128,
    max_points_per_scramble = 128, relative_tolerance = 1e-12)
  testthat::local_mocked_bindings(
    .selection_context = function(...) list(),
    .selection_require_step_evaluable = function(...) NULL,
    .package = "RoBMA"
  )
  expect_error(
    .zplot_selection_marginal(NULL, matrix(0), 0, NULL, "estimate", control),
    'Conditional zplot selection targets are unavailable. Use \'conditioning_depth = "marginal"\'.',
    fixed = TRUE
  )

  prior <- BayesTools::prior_weightfunction(
    "one-sided", .025, BayesTools::wf_fixed(c(1, .2))
  )
  spec <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0), c(1, 1),
    effect_direction = "positive", signed_data = FALSE)
  spec$omega       <- matrix(c(1, .2), 1)
  spec$alpha       <- 0
  spec$phack_kind  <- 0L
  spec$kernel_mode <- SELKERNEL_STEP
  spec$vector_rule <- rep(0L, nrow(spec$omega))
  spec$use_normal  <- FALSE
  failure <- tryCatch(
    .zplot_joint_block(0, matrix(0, 1, 2), matrix(c(1, .7, 1), 1),
      c(1, 1), spec, FALSE, control),
    error = identity
  )
  expect_s3_class(failure, "error")
  message <- conditionMessage(failure)
  observed <- sub("^.*error was ([^ ]+)\\. Increase.*$", "\\1", message)
  expect_gt(as.numeric(observed), control$relative_tolerance)
  expect_identical(message,
    paste0("Zplot marginal integration was rejected by diagnostics: ",
      "relative integration error was ", observed,
      ". Increase 'max_points_per_scramble' ",
      "in 'integration_control = set_selection_likelihood_control()'.")
  )

  # A normalizer-only request uses the shared factor fallback and must honor
  # both of its diagnostics, even when the other criterion is satisfied.
  testthat::local_mocked_bindings(
    .selection_joint_cluster_loglik_block = function(...) {

      list(log_density = 0, relative_mcse = diagnostics[1L],
           relative_change = diagnostics[2L], log_normalizer = log(.8))
    },
    .package = "RoBMA"
  )
  for (diagnostics in list(c(.02, 0), c(0, .02))) {
    failure <- tryCatch(
      .zplot_normalizer_block(
        c(0, 0), matrix(0, 1L, 2L), matrix(c(1, .7, 1), 1L),
        c(1, 1), spec, set_selection_likelihood_control(),
        list(block_methods = "rank_one"), 1L, new.env(),
        list(residual_sd = matrix(1, 1L, 2L), loading = matrix(0, 1L, 2L))
      ),
      error = identity
    )
    expect_s3_class(failure, "error")
    expect_null(conditionCall(failure))
    expect_identical(conditionMessage(failure), paste0(
      "Zplot marginal integration was rejected by diagnostics: ",
      "relative integration error was 0.02. Increase 'max_points_per_scramble' ",
      "in 'integration_control = set_selection_likelihood_control()'."
    ))
  }
})


test_that("zplot reuses invariant predictive components", {

  calls <- new.env(parent = emptyenv())
  calls$predictive <- 0L
  calls$density    <- 0L
  calls$paired     <- 0L
  calls$selection  <- FALSE

  predictive <- list(
    mu         = matrix(c(.1, .2, .3, .4), nrow = 2),
    tau_within = matrix(c(.05, .1, .15, .2), nrow = 2),
    sei        = c(.1, .2)
  )
  object <- list(fit = structure(list(), class = "BayesTools_fit"))

  testthat::local_mocked_bindings(
    .get_posterior_samples = function(...) matrix(0, nrow = 2, ncol = 1),
    .thin_sample_rows = function(...) NULL,
    .is_PET = function(...) FALSE,
    .is_PEESE = function(...) FALSE,
    .is_weightfunction = function(...) calls$selection,
    .effect_direction = function(...) "positive",
    .zplot_predictive_components = function(..., extrapolate) {
      calls$predictive <- calls$predictive + 1L
      if (extrapolate) {
        stop("Invariant predictive components were recomputed.")
      }
      predictive
    },
    .zplot_selection_context = function(...) {
      if (calls$selection) list(selection = TRUE) else NULL
    },
    .zplot_density_vectorized = function(...) {
      calls$density <- calls$density + 1L
      matrix(1, nrow = 1, ncol = 1)
    },
    .zplot_selnorm_density_pair = function(...) {
      calls$paired <- calls$paired + 1L
      list(
        fitted       = matrix(2, nrow = 1, ncol = 1),
        extrapolated = matrix(3, nrow = 1, ncol = 1)
      )
    },
    .package = "RoBMA"
  )

  ordinary <- .zplot_density_pair(object, z_sequence = 0, max_samples = 10)
  expect_identical(calls$predictive, 1L)
  expect_identical(calls$density, 1L)
  expect_identical(calls$paired, 0L)
  expect_identical(ordinary$fitted, ordinary$extrapolated)

  calls$predictive <- 0L
  calls$density    <- 0L
  calls$paired     <- 0L
  calls$selection  <- TRUE
  selected <- .zplot_density_pair(object, z_sequence = 0, max_samples = 10)
  expect_identical(calls$predictive, 1L)
  expect_identical(calls$density, 0L)
  expect_identical(calls$paired, 1L)
  expect_false(identical(selected$fitted, selected$extrapolated))
})


test_that("zplot retains separate PET and PEESE predictive components", {

  calls <- new.env(parent = emptyenv())
  calls$predictive <- 0L
  calls$density    <- 0L

  object <- list(fit = structure(list(), class = "BayesTools_fit"))
  testthat::local_mocked_bindings(
    .get_posterior_samples = function(...) matrix(0, nrow = 2, ncol = 1),
    .thin_sample_rows = function(...) NULL,
    .is_PET = function(...) TRUE,
    .is_PEESE = function(...) FALSE,
    .is_weightfunction = function(...) FALSE,
    .effect_direction = function(...) "positive",
    .zplot_predictive_components = function(..., extrapolate) {
      calls$predictive <- calls$predictive + 1L
      list(
        mu         = matrix(as.numeric(extrapolate), nrow = 2, ncol = 1),
        tau_within = matrix(.1, nrow = 2, ncol = 1),
        sei        = .2
      )
    },
    .zplot_selection_context = function(...) NULL,
    .zplot_density_vectorized = function(..., mu_samples) {
      calls$density <- calls$density + 1L
      mu_samples
    },
    .package = "RoBMA"
  )

  result <- .zplot_density_pair(object, z_sequence = 0, max_samples = 10)
  expect_identical(calls$predictive, 2L)
  expect_identical(calls$density, 2L)
  expect_equal(result$fitted, matrix(0, nrow = 2, ncol = 1))
  expect_equal(result$extrapolated, matrix(1, nrow = 2, ncol = 1))
})


test_that("the chunked marginal assembler keeps only the quantities a chunk carries", {

  S    <- 4L
  rows <- 1:2
  result <- list(
    fitted       = matrix(0, S, 3L),
    extrapolated = matrix(0, S, 3L),
    weights      = numeric(S),
    EDR          = numeric(S)
  )

  full_chunk <- list(
    fitted       = matrix(1, length(rows), 3L),
    extrapolated = matrix(2, length(rows), 3L),
    weights      = c(.5, .25),
    EDR          = c(.8, .6)
  )
  assembled <- .zplot_marginal_chunk_assign(result, full_chunk, rows, TRUE)
  expect_equal(assembled$fitted[rows, ], full_chunk$fitted)
  expect_equal(assembled$extrapolated[rows, ], full_chunk$extrapolated)
  expect_equal(assembled$weights[rows], full_chunk$weights)
  expect_equal(assembled$EDR[rows], full_chunk$EDR)

  # A fitted-only chunk has no extrapolated curve, no inverse weights and no
  # EDR; writing them would fabricate values for rows that were never asked for.
  fitted_chunk <- list(fitted = matrix(1, length(rows), 3L))
  assembled <- .zplot_marginal_chunk_assign(result, fitted_chunk, rows, TRUE)
  expect_equal(assembled$fitted[rows, ], fitted_chunk$fitted)
  expect_null(assembled$extrapolated)
  expect_null(assembled$weights)
  expect_null(assembled$EDR)

  # a density request without probabilities leaves EDR untouched
  no_edr <- list(
    fitted       = matrix(1, length(rows), 3L),
    extrapolated = matrix(2, length(rows), 3L),
    weights      = c(.5, .25)
  )
  assembled <- .zplot_marginal_chunk_assign(result, no_edr, rows, FALSE)
  expect_equal(assembled$weights[rows], no_edr$weights)
  expect_equal(assembled$EDR, numeric(S))
})
