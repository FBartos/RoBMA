.factor_grid_dense_loglik <- function(outcome, mean, covariance) {

  root <- chol(covariance)
  residual <- outcome - mean
  solved <- backsolve(root, forwardsolve(t(root), residual))
  -0.5 * (
    length(outcome) * log(2 * pi) +
      2 * sum(log(diag(root))) +
      sum(residual * solved)
  )
}


test_that("compact factor and Markov grids match dense covariance likelihoods", {

  outcome  <- c(0.2, -0.4, 0.5, 0.1)
  means    <- rbind(c(0.1, -0.2, 0.3, 0), c(-0.1, 0.2, 0.1, 0.15))
  sampling <- matrix(c(
    0.50, 0.08, 0.00, 0.00,
    0.08, 0.60, 0.04, 0.00,
    0.00, 0.04, 0.55, 0.06,
    0.00, 0.00, 0.06, 0.65
  ), nrow = 4L, byrow = TRUE)
  plan_cache <- new.env(parent = emptyenv())
  scales <- rbind(
    c(0.35, 0.45, 0.55, 0.65),
    c(0.50, 0.40, 0.60, 0.30)
  )
  factor_plan <- list(
    type = "group",
    model_matrix = diag(4L),
    group_map = rep(1L, 4L),
    coefficient_structure = "markov"
  )
  baseline_states <- lapply(seq_len(nrow(scales)), function(draw) {
    list(list(
      coefficient_factor = diag(scales[draw, ]),
      coefficient_scale = scales[draw, ],
      markov_transition = rep(0, 3L),
      markov_innovation_variance = rep(1, 3L)
    ))
  })
  rho <- c(-0.25, 0.2, 0.65)
  correlation <- lapply(rho, function(value) {
    value^abs(outer(seq_len(4L), seq_len(4L), "-"))
  })
  cholesky <- array(NA_real_, dim = c(length(rho), 4L, 4L))
  transition <- innovation <- matrix(NA_real_, length(rho), 3L)
  for (grid in seq_along(rho)) {
    current <- t(chol(correlation[[grid]]))
    cholesky[grid, , ] <- current
    transition[grid, ] <-
      current[cbind(2:4, 1:3)] / diag(current)[1:3]
    innovation[grid, ] <- diag(current)[2:4]^2
  }
  update <- list(
    family = "markov",
    factor_index = 1L,
    coefficient_scale = scales,
    candidate_cholesky = cholesky,
    candidate_transition = transition,
    candidate_innovation_variance = innovation
  )

  observed <- .marglik_covariance_plan_factor_grid_loglik(
    cache = plan_cache,
    y = outcome,
    means = means,
    sampling_covariance = sampling,
    random_covariance_plans = list(time = factor_plan),
    random_covariance_states = baseline_states,
    block_indices = list(1:4),
    extra_variances = matrix(0, nrow(means), ncol(means)),
    update_grid = update
  )
  expected <- matrix(NA_real_, length(rho), nrow(means))
  for (grid in seq_along(rho)) {
    for (draw in seq_len(nrow(means))) {
      covariance <- sampling +
        correlation[[grid]] * tcrossprod(scales[draw, ])
      expected[grid, draw] <- .factor_grid_dense_loglik(
        outcome,
        means[draw, ],
        covariance
      )
    }
  }
  expect_equal(observed, expected, tolerance = 1e-12)

  correlation_draws <- array(NA_real_, dim = c(2L, 4L, 4L))
  for (draw in 1:2) {
    correlation_draws[draw, , ] <- t(chol(
      matrix(0.15 * draw, 4L, 4L) + diag(1 - 0.15 * draw, 4L)
    ))
  }
  candidate_scale <- rbind(c(0.2, 0.3), c(0.7, 0.8), c(1.1, 0.9))
  factor_update <- list(
    family = "factor",
    factor_index = 1L,
    component_index = 2L,
    coefficient_scale = scales,
    coefficient_cholesky = correlation_draws,
    candidate_scale = candidate_scale
  )
  factor_plan$coefficient_structure <- "dense"
  baseline_states <- lapply(seq_len(nrow(scales)), function(draw) {
    list(list(coefficient_factor = diag(scales[draw, ])))
  })
  observed <- .marglik_covariance_plan_factor_grid_loglik(
    cache = new.env(parent = emptyenv()),
    y = outcome,
    means = means,
    sampling_covariance = sampling,
    random_covariance_plans = list(time = factor_plan),
    random_covariance_states = baseline_states,
    block_indices = list(1:4),
    extra_variances = matrix(0, nrow(means), ncol(means)),
    update_grid = factor_update
  )
  for (grid in seq_len(nrow(candidate_scale))) {
    for (draw in seq_len(nrow(means))) {
      scale <- scales[draw, ]
      scale[[2L]] <- candidate_scale[grid, draw]
      root <- matrix(correlation_draws[draw, , ], 4L, 4L)
      covariance <- sampling + tcrossprod(root * scale)
      expected[grid, draw] <- .factor_grid_dense_loglik(
        outcome,
        means[draw, ],
        covariance
      )
    }
  }
  expect_equal(observed, expected, tolerance = 1e-12)
})


# A direct SD shared by several columns of a us() block drives all of them: the
# exchangeable unit of a three-level mean-difference factor, and the indicator
# columns of a logical predictor's interaction. The native factor grid writes
# the scale of one column, so their q-grids take the general route. The
# reference is the dense likelihood of the BayesTools covariance of the draws
# with the SD coordinate replaced, plus the declared SD prior; differences to
# the first grid value cancel the prior of the other coordinates.
test_that("known-V q-grids of an SD that drives several us() columns follow the replaced covariance", {

  dat <- data.frame(
    yi    = c(-.2, .1, .3, .5, 0, .2, .4, -.1, .25, .15, -.05, .35),
    group = factor(rep(c("a", "b", "c"), 4L)),
    flag  = rep(c(TRUE, FALSE, TRUE), 4L),
    x     = c(-1, .5, 2, 1, -.5, .3, 1.2, -1.4, .8, .1, -.7, 1.6),
    study = factor(rep(paste0("s", 1:4), each = 3L))
  )
  sd_prior <- BayesTools::prior("gamma", list(shape = 2, rate = 4))
  cases <- list(
    unit   = list(mods = ~ group, random = ~ us(1 + group | study), sd = "(mu) sd(group)", driven = 2:3),
    flag   = list(mods = ~ x, random = ~ us(1 + flag:x | study), sd = "(mu) sd(flag:x)", driven = 2:3),
    # the single-column SD of the same block keeps the native grid
    single = list(mods = ~ group, random = ~ us(1 + group | study), sd = "(mu) sd(intercept)", driven = NULL)
  )
  values <- c(.1, .25, .4, .7, 1.1)

  for (name in names(cases)) {
    case  <- cases[[name]]
    setup <- review_compile_object(brma.mv(
      yi = yi, V = diag(.04, 12L), mods = case$mods, random = case$random, data = dat,
      measure = "GEN", prior_unit_information_sd = 1,
      prior_heterogeneity = BayesTools::prior_random(sd = sd_prior),
      set_contrast_factor_predictors = "meandif", standardize_continuous_predictors = FALSE,
      only_priors = TRUE, silent = TRUE
    ))
    term       <- setup$compiled$mu$formula_design$random_effects[[1L]]
    primitives <- term$correlation$primitive_names
    sds        <- unique(term$sd_parameter_names)
    columns    <- unique(c(review_prior_columns(setup$compiled$mu$prior_list), primitives))
    # zero slopes: every estimate has the intercept as its mean
    samples <- matrix(0, 6L, length(columns), dimnames = list(NULL, columns))
    samples[, "mu_intercept"] <- .1
    samples[, sds]        <- seq(.2, .6, length.out = 6L * length(sds))
    samples[, primitives] <- seq(.3, .7, length.out = 6L * length(primitives))
    object <- review_attach_replay(setup, samples)

    catalog  <- BayesTools::parameter_catalog(object$fit)
    selected <- catalog$quantities[catalog$quantities$canonical_name == case$sd, , drop = FALSE]
    plan     <- BayesTools::random_effects_marginal_update_plan(
      object$fit, BayesTools::parameter_catalog_resolve(catalog, case$sd, "mu")
    )
    expect_identical(plan$family, "factor", info = name)
    expect_identical(plan$driven_columns, case$driven, info = name)

    target  <- .brma_random_parameter_density_target(
      object, .brma_random_parameter_io_labels(selected, "selector")
    )
    context <- .iwmde_context(object)
    spec    <- .iwmde_parameter_spec(context, target$parameter, target$parameter_spec)
    states  <- .iwmde_row_states(context, seq_len(nrow(samples)),
      parameter = target$parameter, parameter_spec = spec)
    replacement <- .iwmde_replacement_spec(context, target$parameter, spec)
    native <- .iwmde_log_q_grid_known_v_random_factor(
      context, target$parameter, values, states, replacement, states[[1L]]$active_setup
    )
    observed <- .iwmde_log_q_grid(context, target$parameter, values, states, replacement)

    expected <- matrix(NA_real_, length(values), nrow(samples))
    for (g in seq_along(values)) {
      candidate <- samples
      candidate[, plan$source_parameter] <- values[[g]]
      covariance <- BayesTools::random_effects_marginal_vcov(
        object$fit, parameter = "mu", posterior_samples = candidate
      )$samples
      for (s in seq_len(nrow(samples))) {
        expected[g, s] <- .factor_grid_dense_loglik(
          dat$yi, rep(.1, 12L), diag(.04, 12L) + covariance[s, , ]
        ) + stats::dgamma(values[[g]], 2, 4, log = TRUE)
      }
    }
    expect_equal(
      sweep(observed, 2L, observed[1L, ]), sweep(expected, 2L, expected[1L, ]),
      tolerance = 1e-10, info = name
    )
    if (is.null(case$driven)) {
      expect_true(is.matrix(native), info = name)
    } else {
      expect_null(native, info = name)
    }
  }
})
