context("Joint selection certified factor routing")
skip_on_cran()

.factor_selection_prior <- function(mode = "integrate", weight_rule = "product") {

  BayesTools::prior_weightfunction(
    "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
    model = BayesTools::selection_model(
      other_random_effects = mode, known_sampling_variance = mode,
      weight_rule = weight_rule, group = "study"
    )
  )
}

test_that("conditional selection accepts equivalent covariance representations without warnings", {

  # Exercise the real pre-fit boundary without running a sampler.
  testthat::local_mocked_bindings(
    .fit = function(object) list(),
    .stop_fit_errors = function(...) NULL,
    .object_summary = function(...) list(),
    .object_coefficients = function(...) list(),
    .refresh_selection_sensitivity_diagnostics = function(object) object,
    .autocompute_brma = function(object) object,
    .package = "RoBMA"
  )
  factor <- known_v_factor(rep(.03, 4L), cbind(
    c(.1, .1, 0, 0), c(0, 0, .1, .1)
  ))
  dense  <- diag(factor$diagonal) + tcrossprod(factor$loading)
  prior <- function(mode = "condition", weights = c(1, .5), prior_weights = 1) {

    BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_fixed(weights),
      prior_weights = prior_weights,
      model = BayesTools::selection_model(
        other_random_effects = mode, known_sampling_variance = mode, group = "study"
      )
    )
  }
  fit <- function(V = dense, prior_bias = prior(), silent = FALSE,
                  only_priors = FALSE, constructor = bselmodel.mv) {

    args <- list(
      yi = c(.1, .2, .3, .4), V = V, measure = "GEN",
      data = data.frame(study = c("a", "a", "b", "b")),
      prior_unit_information_sd = 1, silent = silent,
      only_priors = only_priors
    )
    if (identical(constructor, bselmodel.mv) || identical(constructor, RoBMA.mv)) {
      args$prior_bias <- prior_bias
    }
    warnings <- list()
    object <- withCallingHandlers(do.call(constructor, args), warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
    list(object = object, warnings = warnings)
  }
  plain <- fit()
  quiet <- fit(silent = TRUE)
  mixture <- fit(constructor = RoBMA.mv, prior_bias = list(prior(), prior()))
  for (result in list(plain, quiet, mixture)) {
    expect_length(result$warnings, 0L)
    expect_equal(.known_v_covariance_matrix(.data_known_v_data(result$object$data)), dense)
  }
  expect_identical(plain$object$data, quiet$object$data)
  expect_identical(plain$object$priors, quiet$object$priors)
  expect_length(fit(only_priors = TRUE)$warnings, 0L)
  expect_length(fit(V = diag(diag(dense)))$warnings, 0L)
  expect_length(fit(V = factor)$warnings, 0L)
  expect_length(fit(prior_bias = prior("integrate"))$warnings, 0L)
  expect_length(fit(prior_bias = prior(weights = c(1, 1)))$warnings, 0L)
  for (constructor in list(brma.mv, bPET.mv, bPEESE.mv)) {
    expect_length(fit(constructor = constructor)$warnings, 0L)
  }
})


test_that("ordinary sampling matrices and declared factors preserve vcalc covariance", {

  skip_if_not_installed("metafor")
  dat <- data.frame(
    yi    = seq(-.3, .4, length.out = 8L),
    vi    = seq(.02, .09, length.out = 8L),
    study = rep(c("a", "b"), each = 4L),
    type  = rep(c("x", "x", "y", "y"), 2L),
    obs   = rep(seq_len(4L), 2L)
  )
  V <- metafor::vcalc(
    vi, cluster = study, type = type, obs = obs,
    rho = c(.6, .3), data = dat
  )
  expected <- matrix(as.numeric(V), nrow = nrow(V))
  # Independent covariance identity: within-type correlation .6 and
  # between-type correlation .3, with independent residual variance .4*vi.
  type_root <- t(chol(matrix(c(.6, .3, .3, .6), 2L)))
  loading <- matrix(0, nrow(dat), 4L)
  for (study in seq_len(2L)) {
    rows <- which(dat$study == c("a", "b")[[study]])
    loading[rows, (2L * study - 1L):(2L * study)] <-
      sqrt(dat$vi[rows]) * type_root[match(dat$type[rows], c("x", "y")), ]
  }
  factor <- known_v_factor(.4 * dat$vi, loading)
  expect_equal(diag(factor$diagonal) + tcrossprod(factor$loading),
               expected, tolerance = 1e-14)

  object <- bselmodel.mv(
    yi = yi, V = V, data = dat, measure = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = .factor_selection_prior(weight_rule = "best"),
    only_priors = TRUE, silent = TRUE
  )
  known_V <- .data_known_v_data(object$data)
  expect_identical(.known_v_storage(known_V), "blocks")
  expect_identical(.known_v_selection_metadata(known_V)$origin, "matrix")
  expect_identical(.data_selection_model(object$data)$groups$provenance, "explicit")

  conditional_object <- function(V, constructor = bselmodel.mv, ...) {

    prior_bias <- if (identical(constructor, RoBMA.mv)) {
      .default_prior.bias_alt(
        model_type = "PSMA", measure = "SMD", data = dat,
        prior_unit_information_sd = 1,
        weightfunction_model = BayesTools::selection_model(
          known_sampling_variance = "condition", weight_rule = "best", group = "study")
      )
    } else {
      .factor_selection_prior("condition", weight_rule = "best")
    }
    constructor(
      yi = yi, V = V, data = dat, random = ~ 1 | study, measure = "SMD",
      prior_bias = prior_bias,
      prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE, ...
    )
  }
  for (constructor in list(bselmodel.mv, RoBMA.mv)) {
    for (input in list(V, factor)) {
      candidate <- conditional_object(input, constructor)
      known_V   <- .data_known_v_data(candidate$data)
      sampling  <- .selection_sampling_structure(candidate$data)
      expect_identical(sampling$residual_variance, numeric(nrow(dat)))
      expect_identical(sampling$source_ids, "sampling_error")
      retained <- matrix(0, nrow(dat), sampling$rank)
      for (block in sampling$latent_blocks) {
        retained[block$index, seq.int(block$z_start, block$z_end)] <- block$B
      }
      expect_equal(tcrossprod(retained), expected, tolerance = 1e-14)
      expect_equal(.known_v_covariance_matrix(known_V), expected)
      expect_identical(.data_selection_model(candidate$data)$groups$row_labels, dat$study)
    }
  }
  expect_error(
    conditional_object(V, known_v_residual_fraction = .2),
    "Unused argument in bselmodel.mv(): 'known_v_residual_fraction'",
    fixed = TRUE
  )
})

test_that("known_v_factor preserves exact provenance and structural routing", {

  K <- 6L
  dat <- data.frame(
    yi    = seq(-.3, .7, length.out = K),
    study = factor(rep(c("a", "b"), each = 3L))
  )
  loading <- cbind(
    c(.12, .08, .04, .10, .06, .03),
    c(-.03, .07, .11, .04, .09, .13)
  )
  V_factor <- known_v_factor(rep(.02, K), loading)
  covariance <- diag(.02, K) + tcrossprod(loading)

  factor_object <- bselmodel.mv(
    yi                        = yi,
    V                         = V_factor,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = .factor_selection_prior(),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  factor_setup <- .data_selection_execution_plan(factor_object[["data"]])
  factor_syntax <- .create_model_syntax(
    factor_object[["data"]], factor_object[["priors"]]
  )

  expect_identical(
    .known_v_storage(.data_known_v_data(factor_object[["data"]])),
    "factor"
  )
  expect_equal(
    .selection_joint_sampling_block(factor_setup[["sampling"]], seq_len(K)),
    covariance
  )
  expect_identical(factor_setup[["schema_version"]], 5L)
  expect_identical(factor_setup[["exactness"]], "EF")
  expect_identical(
    factor_setup[["factor_ranks"]],
    2L
  )
  expect_match(factor_syntax, "dselnorm_factor_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", factor_syntax, fixed = TRUE))

  dense_object <- bselmodel.mv(
    yi                        = yi,
    V                         = covariance,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = .factor_selection_prior(),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  dense_setup <- .data_selection_execution_plan(dense_object[["data"]])
  dense_syntax <- .create_model_syntax(
    dense_object[["data"]], dense_object[["priors"]]
  )
  expect_identical(dense_setup[["exactness"]], "E2")
  expect_match(dense_syntax, "dselnorm_mnorm_step", fixed = TRUE)
  expect_false(grepl("dselnorm_factor_step", dense_syntax, fixed = TRUE))
})


test_that("factor routing fails closed outside the certified kernel contract", {

  K <- 6L
  dat <- data.frame(yi = seq(-.3, .7, length.out = K), study = "a",
                    estimate = seq_len(K))
  high_rank <- matrix(seq_len(K * 5L) / 100, nrow = K, ncol = 5L)
  zero_residual <- matrix(c(
    .2, .1,
    .1, .2,
    .2, .1,
    .1, .2,
    .2, .1,
    .1, .2
  ), nrow = K, byrow = TRUE)

  make_object <- function(V) {
    bselmodel.mv(
      yi                        = yi,
      V                         = V,
      random                    = ~ 1 | estimate,
      data                      = dat,
      measure                   = "SMD",
      prior_unit_information_sd = 1,
      prior_bias = .factor_selection_prior(),
      only_priors               = TRUE,
      silent                    = TRUE
    )
  }
  high_rank_object <- make_object(known_v_factor(rep(.02, K), high_rank))
  expect_warning(
    zero_residual_object <- make_object(known_v_factor(
      rep(0, K), zero_residual
    )),
    "positive semidefinite"
  )

  for (object in list(high_rank_object, zero_residual_object)) {
    setup <- .data_selection_execution_plan(object[["data"]])
    syntax <- .create_model_syntax(object[["data"]], object[["priors"]])
    expect_identical(setup[["exactness"]], "E2")
    expect_match(syntax, "dselnorm_mnorm_step", fixed = TRUE)
    expect_false(grepl("dselnorm_factor_step", syntax, fixed = TRUE))
  }

  random_dat <- data.frame(
    yi    = dat[["yi"]],
    study = factor(rep("a", K))
  )
  random_high_rank_object <- bselmodel.mv(
    yi                        = yi,
    V                         = known_v_factor(rep(.02, K), high_rank),
    random                    = ~ 1 | study,
    data                      = random_dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = .factor_selection_prior(),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  random_setup    <- .data_selection_execution_plan(
    random_high_rank_object[["data"]]
  )
  random_fit_data <- .create_fit_data(
    random_high_rank_object[["data"]],
    random_high_rank_object[["priors"]]
  )
  random_syntax   <- .create_model_syntax(
    random_high_rank_object[["data"]],
    random_high_rank_object[["priors"]]
  )
  expect_identical(
    random_setup[["random_covariance"]][["representation"]],
    "diagonal_factor"
  )
  expect_identical(
    random_setup[["block_methods"]],
    "dense"
  )
  expect_true("sel_joint_block_1_diagonal" %in% names(random_fit_data))
  expect_match(random_syntax, "sel_joint_block_1_diagonal", fixed = TRUE)

  expect_error(known_v_factor(c(.1, -.1), matrix(0, 2L, 1L)),
               "non-negative")
  expect_error(known_v_factor(c(.1, .1), matrix(0, 3L, 1L)),
               "one row")
  expect_error(known_v_factor(c(.1, .1), matrix(c(0, NA), 2L, 1L)),
               "finite")
})


test_that("higher-rank factor controls retain budget semantics", {

  default_control <- set_selection_likelihood_control()
  odd_control <- set_selection_likelihood_control(
    points_per_scramble     = 9L,
    max_points_per_scramble = 9L,
    scrambles               = 2L
  )
  sampling <- structure(list(
    representation = "diagonal_factor",
    diagonal       = rep(1, 4L),
    loading        = matrix(0, nrow = 4L, ncol = 0L)
  ), class = c("RoBMA_selection_sampling_plan", "list"))
  factor_plan <- .selection_joint_execution_plan(
    row_blocks       = list(1:4),
    block_methods     = "factor",
    factor_ranks      = 3L,
    selection_control = odd_control,
    sampling          = sampling,
    sampling_factor_blocks = NULL,
    random_covariance = NULL
  )
  dense_plan <- .selection_joint_execution_plan(
    row_blocks       = list(1:4),
    block_methods     = "dense",
    factor_ranks      = NA_integer_,
    selection_control = odd_control,
    sampling          = sampling,
    sampling_factor_blocks = NULL,
    random_covariance = NULL
  )

  expect_identical(default_control[["points_per_scramble"]], 512L)
  expect_identical(default_control[["max_points_per_scramble"]], 8192L)
  expect_identical(default_control[["scrambles"]], 8L)
  expect_identical(default_control[["relative_tolerance"]], .005)
  expect_identical(default_control[["seed"]], 1L)
  expect_identical(factor_plan[["points_per_scramble"]], 9L)
  expect_identical(factor_plan[["max_points_per_scramble"]], 9L)
  expect_identical(factor_plan[["factor_points_per_proposal"]], 5L)
  expect_identical(factor_plan[["factor_max_points_per_proposal"]], 5L)
  expect_identical(
    dim(factor_plan[["designs"]][["factor_3"]]),
    c(2L, 5L, 6L)
  )
  expect_identical(
    dim(dense_plan[["designs"]][["4"]]),
    c(2L, 9L, 8L)
  )
})


test_that("post-fit likelihood retains certified sampling factors", {

  K <- 4L
  loading <- cbind(
    c(.12, .08, .04, .10),
    c(-.03, .07, .11, .04)
  )
  object <- bselmodel.mv(
    yi                        = seq(-.3, .6, length.out = K),
    data                      = data.frame(study = rep("a", K), estimate = seq_len(K)),
    V                         = known_v_factor(rep(.02, K), loading),
    random                    = ~ 1 | estimate,
    prior_heterogeneity       = BayesTools::prior("point", list(location = .15)),
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = .factor_selection_prior(),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  factor_call <- NULL
  testthat::local_mocked_bindings(
    .estimate_normal_covariance_target_location_from_setup = function(setup) {
      list(
        y     = seq(-.3, .6, length.out = K),
        means = matrix(0, nrow = 2L, ncol = K)
      )
    },
    .selection_joint_signed_context = function(setup, signed_yi) {
      list(obs_bin = rep(1L, K))
    },
    .selection_joint_random_covariance_samples = function(setup) NULL,
    .selection_joint_covariance_lower = function(...) {
      stop("Certified post-fit factors must not use dense covariance.")
    },
    .selection_joint_factor_loglik_block = function(
        yi, means, residual_sd, loading, sei, selection_context,
        execution_plan, block_index) {
      factor_call <<- list(
        residual_sd = residual_sd,
        loading     = loading,
        block_index = block_index
      )
      c(11, 12)
    },
    .package = "RoBMA"
  )

  observed <- .selection_joint_block_loglik_from_setup(list(
    data          = object[["data"]],
    priors        = object[["priors"]],
    fit           = structure(list(), formula_design = object$formula_design,
      prior_list = c(object$formula_design$mu$prior_list,
                     .create_fit_priors(object$data, object$priors))),
    posterior_samples = matrix(0, 2L, 1L, dimnames = list(NULL, "mu_intercept")),
    S             = 2L,
    K             = K,
    selection_sei = sqrt(.known_v_diagonal(
      .data_known_v_data(object[["data"]])
    )),
    tau_within    = matrix(0, nrow = 2L, ncol = K),
    tau_between   = matrix(0, nrow = 2L, ncol = K),
    is_multilevel = FALSE
  ))

  expect_equal(observed, matrix(c(11, 12), nrow = 2L))
  expect_identical(factor_call[["block_index"]], 1L)
  expect_equal(
    factor_call[["residual_sd"]],
    matrix(sqrt(.02 + .15^2), nrow = 2L, ncol = K),
    tolerance = 0
  )
  expect_equal(
    matrix(factor_call[["loading"]], nrow = 2L),
    matrix(array(
      rep(loading, each = 2L),
      dim = c(2L, K, ncol(loading))
    ), nrow = 2L),
    tolerance = 0
  )
})


test_that("exact random factors reuse repeated states without changing covariance", {

  dat <- data.frame(yi = c(-.2, .1), vi = c(.02, .03), study = "a")
  object <- bselmodel.mv(
    yi = yi, vi = vi, random = ~ 1 | study, data = dat, measure = "SMD",
    prior_unit_information_sd = 1, prior_bias = .factor_selection_prior(),
    only_priors = TRUE, silent = TRUE
  )
  factor_plans <- list(study = list(
    type = "group", model_matrix = matrix(1, 2L, 1L),
    group_map = c(1L, 1L), coefficient_structure = "diagonal"
  ))
  built_rows <- integer()
  testthat::local_mocked_bindings(
    .brma_mv_random_effects_marginal_factor_states = function(
        object, posterior_samples, blocks, row_blocks) {

      built_rows <<- c(built_rows, nrow(posterior_samples))
      structure(list(
        factor_plans = factor_plans,
        factor_states = lapply(posterior_samples[, "sd"], function(sd) {
          list(study = list(coefficient_factor = matrix(sd, 1L, 1L)))
        }),
        row_blocks = row_blocks,
        metadata = list(n_draws = nrow(posterior_samples), n_rows = 2L,
                        included_blocks = "study")
      ), class = c("BayesTools_random_effects_marginal_factor_states", "list"))
    },
    .package = "RoBMA"
  )
  setup <- list(data = object[["data"]], priors = object[["priors"]],
                S = 2L, K = 2L, posterior_samples = cbind(sd = c(.2, .5)))
  factors <- .selection_joint_random_factor_samples(setup)
  rows <- c(2L, 1L, 2L)
  factors[["diagonal"]] <- factors[["diagonal"]][rows, , drop = FALSE]
  factors[["loadings"]] <- lapply(factors[["loadings"]], function(loading) {
    loading[rows, , , drop = FALSE]
  })
  setup[["S"]] <- 3L
  setup[["posterior_samples"]] <- setup[["posterior_samples"]][rows, , drop = FALSE]
  setup[["selection_random_factor_samples"]] <- factors
  reused <- .selection_joint_random_factor_samples(setup)
  expect_identical(built_rows, 2L)
  for (i in seq_along(rows)) {
    covariance <- diag(reused[["diagonal"]][i, ]) +
      tcrossprod(matrix(reused[["loadings"]][[1L]][i, , ], nrow = 2L))
    expect_equal(covariance, matrix(c(.2, .5)[rows[[i]]]^2, 2L, 2L),
                 tolerance = 1e-14)
  }
  setup[["selection_random_factor_samples"]] <- NULL
  setup[["posterior_samples"]][, "sd"] <- c(.3, .4, .6)
  changed <- .selection_joint_random_factor_samples(setup)
  expect_identical(built_rows, c(2L, 3L))
  expect_equal(as.numeric(changed[["loadings"]][[1L]][, 1L, 1L]),
               c(.3, .4, .6), tolerance = 1e-14)
})


test_that("dense covariance assembly preserves diagonal and nonempty factor parts", {

  dat      <- data.frame(yi = c(-.2, .1, .3), study = "a", esid = seq_len(3L))
  sampling <- diag(c(.02, .03, .04)) + tcrossprod(c(.04, .02, .03))
  diagonal <- rbind(rep(.04, 3L), rep(.09, 3L))
  for (rank in 0:1) {
    object <- bselmodel.mv(
      yi                        = yi,
      V                         = sampling,
      random                    = ~ 1 | study/esid,
      data                      = dat,
      measure                   = "SMD",
      prior_unit_information_sd = 1,
      selection = BayesTools::selection_model(
        other_random_effects = if (rank == 0L) "condition" else "integrate",
        group                = "study"
      ),
      only_priors = TRUE, silent = TRUE
    )
    plan <- .data_selection_execution_plan(object[["data"]])
    # Dense JAGS storage has no factor columns; post-fit factors may still
    # contain a nonempty loading component, independently of that storage.
    expect_identical(plan[["random_covariance"]][["representation"]], "dense")
    expect_identical(plan[["random_covariance"]][["loading_ranks"]][[1L]], 0L)
    loading <- array(0, dim = c(2L, 3L, rank))
    if (rank > 0L) {
      loading[1L, , 1L] <- .3
      loading[2L, , 1L] <- .4
    }
    factors <- list(
      diagonal         = diagonal,
      loadings         = list(loading),
      ranks            = rank,
      loading_supports = list(matrix(TRUE, 3L, rank)),
      row_blocks       = list(seq_len(3L))
    )
    observed <- .selection_joint_covariance_lower(
      setup                 = list(data = object[["data"]], S = 2L),
      block_index           = 1L,
      random_factor_samples = factors
    )
    expected <- t(vapply(seq_len(2L), function(draw) {

      covariance <- sampling + diag(diagonal[draw, ])
      if (rank > 0L) {
        covariance <- covariance + tcrossprod(rep(c(.3, .4)[draw], 3L))
      }
      covariance[lower.tri(covariance, diag = TRUE)]
    }, numeric(6L)))
    expect_equal(observed, expected, tolerance = 1e-14)
  }
})

test_that("bridge factor states retain the certified covariance exactly", {

  dat <- data.frame(
    yi    = c(-.20, .05, .35),
    vi    = c(.012, .018, .025),
    study = factor(rep("a", 3L)),
    time  = 1:3
  )
  object <- bselmodel.mv(
    yi                        = yi,
    vi                        = vi,
    random                    = ~ har(time | study),
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = .factor_selection_prior(),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  exact_setup <- .data_selection_execution_plan(object[["data"]])
  random_term <- object[["formula_design"]][["mu"]][["random_effects"]][[1L]]
  coefficient_factor <- matrix(
    c(.25, 0, 0, .08, .22, 0, -.03, .07, .18),
    nrow = 3L,
    byrow = TRUE
  )
  bridge_factor <- list(
    representation = "factor_state",
    row_blocks     = exact_setup[["row_blocks"]],
    factor_plans   = list(study = list(
      type                  = "group",
      model_matrix          = random_term[["model_matrix"]],
      group_map             = random_term[["group_map"]],
      coefficient_structure = "markov"
    )),
    factor_states  = list(study = list(
      coefficient_factor = coefficient_factor
    ))
  )
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  bridge_data <- .marglik_add_selection_bridge_data(
    fit_data         = fit_data,
    priors           = object[["priors"]],
    effect_direction = .data_effect_direction(object[["data"]]),
    model_data       = object[["data"]]
  )
  factor_call <- NULL
  testthat::local_mocked_bindings(
    .marglik_bridge_random_covariance = function(...) bridge_factor,
    .marglik_selection_context = function(parameters, data) {
      list(obs_bin = rep(1L, 3L))
    },
    .selection_joint_factor_loglik_block = function(
        yi, means, residual_sd, loading, sei, selection_context,
        execution_plan, block_index) {
      factor_call <<- list(
        residual_sd = residual_sd,
        loading     = loading,
        block_index = block_index
      )
      17
    },
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    random_effects_marginal_factor_vcov = function(...) {
      stop("Certified bridge factors must not be materialized as dense.")
    },
    .package = "BayesTools"
  )

  observed <- .marglik_joint_selection_log_lik(
    parameters            = list(),
    data                  = bridge_data,
    model_data            = object[["data"]],
    bridge_context        = list(),
    covariance_plan_cache = NULL,
    mu_samples            = matrix(0, nrow = 1L, ncol = 3L),
    tau_within_samples    = matrix(0, nrow = 1L, ncol = 3L),
    tau_between_samples   = matrix(0, nrow = 1L, ncol = 3L),
    is_random             = TRUE,
    is_multilevel         = FALSE,
    fixed_zero_random     = FALSE,
    K                     = 3L
  )

  rank <- exact_setup[["factor_ranks"]][[1L]]
  loading <- matrix(factor_call[["loading"]], nrow = 3L, ncol = rank)
  reconstructed <- diag(as.numeric(factor_call[["residual_sd"]])^2) +
    tcrossprod(loading)
  basis <- random_term[["model_matrix"]] %*% coefficient_factor
  expected <- diag(dat$vi) + tcrossprod(basis)

  expect_identical(observed, 17)
  expect_identical(factor_call[["block_index"]], 1L)
  expect_identical(rank, 2L)
  expect_equal(unname(reconstructed), unname(expected), tolerance = 1e-14)
})


test_that("bridge dense fallback reconstructs the exact required block", {

  dat <- data.frame(
    yi    = c(-.20, .05, .35),
    study = factor(rep("a", 3L))
  )
  sampling <- diag(c(.02, .03, .04)) + tcrossprod(c(.03, .02, .01))
  object <- bselmodel.mv(
    yi                        = yi,
    V                         = sampling,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = .factor_selection_prior(),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  exact_setup <- .data_selection_execution_plan(object[["data"]])
  random_term <- object[["formula_design"]][["mu"]][["random_effects"]][[1L]]
  bridge_factor <- list(
    representation = "factor_state",
    row_blocks     = exact_setup[["row_blocks"]],
    factor_plans   = list(study = list(
      type                  = "group",
      model_matrix          = random_term[["model_matrix"]],
      group_map             = random_term[["group_map"]],
      coefficient_structure = "diagonal"
    )),
    factor_states  = list(study = list(
      coefficient_factor = matrix(.2, 1L, 1L)
    ))
  )
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  bridge_data <- .marglik_add_selection_bridge_data(
    fit_data         = fit_data,
    priors           = object[["priors"]],
    effect_direction = .data_effect_direction(object[["data"]]),
    model_data       = object[["data"]]
  )
  observed_lower <- NULL
  testthat::local_mocked_bindings(
    .marglik_bridge_random_covariance = function(...) bridge_factor,
    .marglik_selection_context = function(parameters, data) {
      list(obs_bin = rep(1L, 3L))
    },
    .selection_joint_dense_loglik_block = function(
        yi, means, covariance_lower, sei, selection_context,
        execution_plan, block_size, normalizer_grid = NULL) {
      expect_null(normalizer_grid)
      observed_lower <<- covariance_lower
      23
    },
    .package = "RoBMA"
  )

  observed <- .marglik_joint_selection_log_lik(
    parameters            = list(),
    data                  = bridge_data,
    model_data            = object[["data"]],
    bridge_context        = list(),
    covariance_plan_cache = NULL,
    mu_samples            = matrix(0, nrow = 1L, ncol = 3L),
    tau_within_samples    = matrix(0, nrow = 1L, ncol = 3L),
    tau_between_samples   = matrix(0, nrow = 1L, ncol = 3L),
    is_random             = TRUE,
    is_multilevel         = FALSE,
    fixed_zero_random     = FALSE,
    K                     = 3L
  )

  expected <- sampling + matrix(.04, 3L, 3L)
  pairs <- .selection_joint_lower_pairs(
    exact_setup,
    1:3
  )
  expected_lower <- matrix(
    expected[cbind(pairs[["row_1"]], pairs[["row_2"]])],
    nrow = 1L
  )
  expect_identical(exact_setup[["exactness"]], "E2")
  expect_identical(observed, 23)
  expect_equal(observed_lower, expected_lower, tolerance = 1e-14)
})


test_that("native factor likelihood agrees with an independent MVN oracle", {

  reference_loglik <- function(y, mu, covariance, sei, omega, z_lower,
                               z_upper, obs_bin, sign) {
    assignments <- expand.grid(rep(list(seq_along(omega)), length(y)))
    terms <- apply(assignments, 1L, function(bins) {
      lower <- if (sign == 1L) {
        z_lower[bins] * sei
      } else {
        -z_upper[bins] * sei
      }
      upper <- if (sign == 1L) {
        z_upper[bins] * sei
      } else {
        -z_lower[bins] * sei
      }
      probability <- suppressWarnings(as.numeric(mvtnorm::pmvnorm(
        lower     = lower,
        upper     = upper,
        mean      = mu,
        sigma     = covariance,
        algorithm = mvtnorm::Miwa(steps = 4096L)
      )))
      prod(omega[bins]) * probability
    })
    mvtnorm::dmvnorm(y, mu, covariance, log = TRUE) +
      sum(log(omega[obs_bin])) - log(sum(terms))
  }

  set.seed(927)
  errors <- diagnostics <- normal_errors <- numeric(3L)
  for (rank in 2:4) {
    K <- 4L
    residual_sd <- seq(.28, .48, length.out = K)
    loading <- matrix(stats::rnorm(K * rank, sd = .22), K, rank)
    loading[, 1L] <- loading[, 1L] + seq(.15, .35, length.out = K)
    mu  <- seq(-.35, .45, length.out = K)
    y   <- c(-.42, -.05, .35, .91)
    sei <- seq(.20, .36, length.out = K)
    sign     <- if (rank == 3L) -1L else 1L
    omega    <- c(.18, .45, 1)
    z_lower  <- stats::qnorm(c(.025, .10, 1), lower.tail = FALSE)
    z_upper  <- stats::qnorm(c(0, .025, .10), lower.tail = FALSE)
    score    <- sign * y / sei
    obs_bin  <- vapply(score, function(value) {
      which(value >= z_lower)[[1L]]
    }, integer(1L))
    qmc <- BayesTools::selection_qmc_design(
      dimensions = 2L * rank,
      points     = 4096L,
      scrambles  = 16L,
      seed       = 417L
    )
    quadrature <- .selection_joint_factor_quadrature_rules(rank)[[
      as.character(rank)
    ]]
    observed <- .Call(
      "RoBMA_selnorm_factor_step_loglik_batch",
      as.double(y),
      matrix(mu, nrow = 1L),
      matrix(residual_sd, nrow = 1L),
      matrix(as.double(loading), nrow = 1L),
      as.double(sei),
      matrix(omega, nrow = 1L),
      as.double(z_lower),
      as.double(z_upper),
      as.integer(obs_bin),
      sign,
      TRUE,
      1L,
      quadrature[["nodes"]],
      quadrature[["log_weights"]],
      as.double(quadrature[["orders"]]),
      as.double(quadrature[["rule_counts"]]),
      as.double(qmc),
      4096L,
      4096L,
      16L,
      .005,
      FALSE,
      0L, PACKAGE = "RoBMA"
    )
    covariance <- diag(residual_sd^2) + tcrossprod(loading)
    normal_observed <- .Call(
      "RoBMA_selnorm_factor_step_loglik_batch",
      as.double(y), matrix(mu, nrow = 1L),
      matrix(residual_sd, nrow = 1L),
      matrix(as.double(loading), nrow = 1L),
      as.double(sei), matrix(omega, nrow = 1L),
      as.double(z_lower), as.double(z_upper), as.integer(obs_bin),
      sign, TRUE, 0L, quadrature[["nodes"]],
      quadrature[["log_weights"]], as.double(quadrature[["orders"]]),
      as.double(quadrature[["rule_counts"]]),
      as.double(qmc), 4096L, 4096L, 16L, .005, FALSE,
      0L, PACKAGE = "RoBMA"
    )
    reference <- reference_loglik(
      y, mu, covariance, sei, omega, z_lower, z_upper, obs_bin, sign
    )
    errors[[rank - 1L]] <- observed[["log_density"]] - reference
    normal_errors[[rank - 1L]] <- normal_observed[["log_density"]] -
      mvtnorm::dmvnorm(y, mu, covariance, log = TRUE)
    diagnostics[[rank - 1L]] <- max(
      observed[["relative_mcse"]], observed[["relative_change"]]
    )
  }

  expect_true(all(abs(errors) < .005))
  expect_equal(normal_errors, rep(0, 3L), tolerance = 1e-11)
  expect_true(all(diagnostics < .005))

  K <- 4L
  residual_sd <- c(.28, .32, .36, .40)
  loading <- cbind(
    c(.25, .20, .15, .10),
    c(.12, .08, .10, .05),
    c(0, .09, .07, .06),
    c(0, 0, 0, .08)
  )
  mu      <- c(-.20, -.05, .15, .35)
  y       <- c(-.31, .02, .29, .61)
  sei     <- c(.20, .24, .28, .32)
  omega   <- c(.35, 1)
  z_lower <- c(stats::qnorm(.025, lower.tail = FALSE), -Inf)
  z_upper <- c(Inf, stats::qnorm(.025, lower.tail = FALSE))
  obs_bin <- ifelse(y / sei >= z_lower[[1L]], 1L, 2L)
  quadrature <- .selection_joint_factor_quadrature_rules(4L)[["4"]]
  qmc <- BayesTools::selection_qmc_design(
    dimensions = 8L,
    points     = 8L,
    scrambles  = 2L,
    seed       = 419L
  )
  observed <- .Call(
    "RoBMA_selnorm_factor_step_loglik_batch",
    y, matrix(mu, nrow = 1L), matrix(residual_sd, nrow = 1L),
    matrix(as.double(loading), nrow = 1L), sei, matrix(omega, nrow = 1L),
    z_lower, z_upper, as.integer(obs_bin), 1L, TRUE, 1L,
    quadrature[["nodes"]], quadrature[["log_weights"]],
    as.double(quadrature[["orders"]]),
    as.double(quadrature[["rule_counts"]]), as.double(qmc),
    8L, 8L, 2L, .005, FALSE,
    0L, PACKAGE = "RoBMA"
  )
  covariance <- diag(residual_sd^2) + tcrossprod(loading)
  reference <- reference_loglik(
    y, mu, covariance, sei, omega, z_lower, z_upper, obs_bin, 1L
  )

  expect_equal(observed[["relative_mcse"]], 0, tolerance = 0)
  expect_lte(observed[["relative_change"]], .005)
  expect_lt(abs(observed[["log_density"]] - reference), .005)
})


test_that("factor quadrature resolves an Assink state before QMC fallback", {

  y <- c(
    2.2844, 2.1771, 1.7777, 1.5480, 1.4855,
    1.4836, 1.2777, 1.0311, .9409, .6263
  )
  vi <- c(
    .3325, .3073, .2697, .4533, .1167,
    .1706, .1538, .3132, .1487, .2139
  )
  mean <- rep(-.462179, length(y))
  residual_sd <- sqrt(.3 * vi + .315042^2)
  loading <- cbind(sqrt(.7 * vi), rep(.245830, length(y)))
  sei <- sqrt(vi)
  omega <- c(1, .562858)
  z_lower <- c(stats::qnorm(.025, lower.tail = FALSE), -Inf)
  z_upper <- c(Inf, stats::qnorm(.025, lower.tail = FALSE))
  obs_bin <- ifelse(y / sei >= z_lower[[1L]], 1L, 2L)
  quadrature <- .selection_joint_factor_quadrature_rules(2L)[["2"]]
  qmc <- BayesTools::selection_qmc_design(
    dimensions = 4L,
    points     = 8L,
    scrambles  = 2L,
    seed       = 1L
  )
  observed <- .Call(
    "RoBMA_selnorm_factor_step_loglik_batch",
    y, matrix(mean, nrow = 1L), matrix(residual_sd, nrow = 1L),
    matrix(as.double(loading), nrow = 1L), sei, matrix(omega, nrow = 1L),
    z_lower, z_upper, as.integer(obs_bin), 1L, TRUE, 1L,
    quadrature[["nodes"]], quadrature[["log_weights"]],
    as.double(quadrature[["orders"]]),
    as.double(quadrature[["rule_counts"]]), as.double(qmc),
    8L, 8L, 2L, .005, FALSE,
    0L, PACKAGE = "RoBMA"
  )
  covariance <- diag(residual_sd^2) + tcrossprod(loading)
  normal_log_density <- mvtnorm::dmvnorm(
    y,
    mean       = mean,
    sigma      = covariance,
    log        = TRUE
  )
  observed_log_normalizer <- normal_log_density +
    sum(log(omega[obs_bin])) - observed[["log_density"]]

  expect_equal(observed[["relative_mcse"]], 0, tolerance = 0)
  expect_lte(observed[["relative_change"]], .005)
  expect_equal(
    observed_log_normalizer,
    -5.5163303949665297,
    tolerance = 5e-5
  )
})


test_that("factor integration retains inactive and opposing selection modes", {

  # Captured valid Assink product-space states: three sampling factors and an
  # inactive random factor. A rank-four budget previously rejected both states.
  y <- c(.7156, .7067, .6475, .6428, .6271, .6238, .6025, .5763,
         .5171, -.3797, -.4228, -.4245, -.4671, -.5230, -.5675, -.7586)
  vi <- c(.0914, .0875, .0330, .0861, .0400, .0680, .1287, .0332,
          .0517, .0390, .0664, .0809, .0667, .0988, .0340, .0437)
  loading <- cbind(
    c(.252942681254074, .176776695296637, .108562029668362,
      .175356779167502, .119522860933439, .155838744494796,
      .300149962518738, .108890508572340, .135883353337654,
      .118019368870416, .153994434036707, .169978990298381,
      .154341920978808, .262982889177224, .110194633003868,
      .124928551008738),
    c(0, .173205080756888, .0443202630213959, .171813852759316,
      .0487950036474267, .0636209010280352, 0, .0444543639723985,
      .133137952086226, .115634893399132, .150883114647446,
      .166544717289810, .151223580927617, 0, .107968249301092,
      .122404481710668),
    c(0, 0, .0966953980290686, 0, .106458129484475,
      .138804418757713, 0, .0969879717628257, rep(0, 8L)),
    0
  )
  weights <- rbind(
    c(1, .288785085821257, .067972879147966,
      .067972879147966, .288785085821257, 1),
    c(1, .557347309764094, .483225646479775,
      .483225646479775, .557347309764094, 1)
  )
  z_lower <- stats::qnorm(c(.025, .05, .5, .95, .975, 1),
                           lower.tail = FALSE)
  z_upper <- c(Inf, head(z_lower, -1L))
  obs_bin <- vapply(y / sqrt(vi), function(z) {
    which(z >= z_lower)[[1L]]
  }, integer(1L))
  quadrature <- .selection_joint_factor_quadrature_rules(3:4)
  evaluate <- function(columns, initial_points = 256L, max_points = 4096L) {
    rank <- length(columns)
    rules <- quadrature[[as.character(rank)]]
    qmc <- BayesTools::selection_qmc_design(
      dimensions = 2L * rank, points = max_points, scrambles = 8L, seed = 1L
    )
    .Call(
      "RoBMA_selnorm_factor_step_loglik_batch",
      y, matrix(0, 2L, length(y)),
      matrix(rep(sqrt(.3 * vi), each = 2L), 2L),
      matrix(rep(as.double(loading[, columns]), each = 2L), 2L),
      sqrt(vi), weights, z_lower, z_upper, obs_bin, 1L, TRUE, 1L,
      rules[["nodes"]], rules[["log_weights"]],
      as.double(rules[["orders"]]), as.double(rules[["rule_counts"]]),
      as.double(qmc), initial_points, max_points, 8L, .005, TRUE,
      0L, PACKAGE = "RoBMA"
    )
  }
  actual <- evaluate(1:4)
  reduced <- evaluate(1:3)

  expect_equal(actual, reduced, tolerance = 1e-12)
  expect_identical(actual[["relative_mcse"]], c(0, 0))
  expect_true(all(actual[["relative_change"]] <= .005))
  # Independent full tensor integration of base-R interval probabilities at
  # orders 95 and 127 agreed within 8e-13; it did not use the native reduction.
  expect_lt(max(abs(actual[["log_normalizer"]] -
                      c(-8.2838610492094045, -7.0139713230657499))),
            3e-5)

  # A small active fourth factor still has both opposing selection modes.
  # Independent base-R interval quadrature (95 sampling nodes, 5 study nodes)
  # gave -8.2656437531115206; a single-mode approximation loses half the mass.
  loading[, 4L] <- .01
  active <- evaluate(1:4, initial_points = 257L, max_points = 4103L)
  expect_true(all(active[["relative_mcse"]] > 0))
  expect_true(all(active[["relative_mcse"]] <= .005))
  expect_true(all(active[["relative_change"]] <= .005))
  expect_lt(abs(active[["log_normalizer"]][[1L]] + 8.2656437531115206), .005)

  # Four coarse points cannot cover all seven proposals. This valid minimum
  # R-interface budget must produce an estimate and honest diagnostics.
  limited <- evaluate(1:4, initial_points = 4L, max_points = 4L)
  expect_true(all(is.finite(limited[["log_normalizer"]])))
  expect_true(all(limited[["relative_mcse"]] > .005))
})


test_that("JAGS instantiates the certified factor distribution", {

  skip_if_not_installed("rjags")
  qmc <- BayesTools::selection_qmc_design(
    dimensions = 4L,
    points     = 8L,
    scrambles  = 2L,
    seed       = 19L
  )
  quadrature <- .selection_joint_factor_quadrature_rules(2L)[["2"]]
  model_text <- paste0(
    "model{\n",
    "  y[1:3] ~ dselnorm_factor_step(",
    "mu[1:3],residual_sd[1:3],loading[1:3,1:2],",
    "sei[1:3],omega,z_lower,z_upper,obs_bin[1:3],",
    "1,1,1,nodes,log_weights,orders,rule_counts,",
    "qmc[1:2,1:8,1:4],8,8,2,0.005,0)\n",
    "}\n"
  )
  connection <- textConnection(model_text)
  on.exit(close(connection), add = TRUE)

  model <- rjags::jags.model(
    file = connection,
    data = list(
      y           = c(-.2, .1, .4),
      mu          = c(0, .05, .1),
      residual_sd = c(.2, .25, .3),
      loading     = matrix(c(.1, .05, .03, -.04, .02, .08), 3L, 2L),
      sei         = c(.15, .2, .25),
      omega       = 1,
      z_lower     = -1e300,
      z_upper     = 1e300,
      obs_bin     = rep(1L, 3L),
      nodes       = quadrature[["nodes"]],
      log_weights = quadrature[["log_weights"]],
      orders      = quadrature[["orders"]],
      rule_counts = quadrature[["rule_counts"]],
      qmc         = qmc
    ),
    n.chains = 1L,
    n.adapt  = 0L,
    quiet    = TRUE
  )

  expect_s3_class(model, "jags")
})


test_that("the early dense rule stays outside rank-one and sampling-conditioned plans", {

  testthat::local_mocked_bindings(
    .fit = function(object) list(),
    .stop_fit_errors = function(...) NULL,
    .object_summary = function(...) list(),
    .object_coefficients = function(...) list(),
    .autocompute_brma = function(object) object,
    .package = "RoBMA"
  )
  data("dat.assink2016", package = "metadat", envir = environment())
  dat <- dat.assink2016
  V <- metafor::vcalc(vi, cluster = study, type = deltype, obs = esid,
    rho = c(.7, .5), data = dat)
  make_plan <- function(sampling = "integrate") {

    prior <- BayesTools::prior_weightfunction("one-sided", steps = .025,
      weights = BayesTools::wf_fixed(c(1, .5)),
      model = BayesTools::selection_model(known_sampling_variance = sampling,
        group = "study"))
    object <- bselmodel.mv(yi = yi, V = V, random = ~ 1 | study / esid,
      data = dat, measure = "SMD", prior_unit_information_sd = 1,
      prior_bias = prior, only_priors = TRUE, silent = TRUE)
    .data_selection_execution_plan(object[["data"]])
  }
  original <- c(15L, 31L, 63L, 127L, 255L, 511L, 1023L)
  expect_identical(SELNORM_CLUSTER_QUADRATURE_ORDERS, original)
  dense <- make_plan()
  expect_true(any(dense[["block_methods"]] == "dense"))
  expect_false(any(dense[["block_methods"]] == "rank_one"))
  expect_identical(dense[["quadrature"]][["orders"]], c(7L, original))
  conditioned <- make_plan("condition")
  expect_identical(conditioned[["statistical_target"]], "whole_sampling_error_selection")
  expect_identical(conditioned[["quadrature"]][["orders"]], original)

  # The compiled rank-one route also owns mixed-sign loadings. A shared plan
  # containing that route must not inherit a dense-only experimental schedule.
  for (methods in list("rank_one", c("dense", "rank_one"))) {
    plan <- .selection_joint_execution_plan(
      row_blocks = lapply(seq_along(methods), function(index) (3L * index - 2L):(3L * index)),
      block_methods = methods,
      factor_ranks = ifelse(methods == "rank_one", 1L, NA_integer_),
      selection_control = set_selection_likelihood_control(
        points_per_scramble = 8L, max_points_per_scramble = 8L, scrambles = 2L),
      sampling = NULL, sampling_factor_blocks = NULL, random_covariance = NULL)
    expect_identical(plan[["quadrature"]][["orders"]], original)
  }
})

# Append to the existing test-00-selection-factor-routing.R after root review.
# Constructor-only regressions; no new fitted-model cache or sampling budget.
test_that("factor cancellation preserves publication and prediction partitions", {

  construct <- function(V, dat, group = "paper", rule = "best") {

    bselmodel.mv(
      yi = dat$yi, V = V, data = dat, measure = "GEN",
      prior_bias = BayesTools::prior_weightfunction(
        "one-sided", steps = .025, weights = BayesTools::wf_fixed(c(1, .5)),
        model = selection_model(group = group, weight_rule = rule)),
      prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  }

  # The columns overlap every row, but their cross-publication covariance
  # cancels exactly. These are positive-definite, correlated 2x2 blocks.
  loading <- rbind(c(1, 1), c(2, 2), c(1, -1), c(2, -2))
  factor <- known_v_factor(c(1, 2, 3, 4), loading)
  covariance <- matrix(c(
    3, 4, 0, 0,
    4, 10, 0, 0,
    0, 0, 5, 4,
    0, 0, 4, 12), 4L, byrow = TRUE)
  expect_identical(diag(factor$diagonal) + tcrossprod(factor$loading), covariance)
  dat <- data.frame(yi = c(.1, -.2, .3, -.4), paper = c("a", "a", "b", "b"))
  dense_fit <- construct(covariance, dat)
  factor_fit <- construct(factor, dat)
  expected <- list(1:2, 3:4)
  for (fit in list(dense_fit, factor_fit)) {
    known <- .data_known_v_data(fit$data)
    plan <- .data_selection_execution_plan(fit$data)
    expect_identical(lapply(.known_v_blocks(known), `[[`, "index"), expected)
    expect_identical(.known_v_covariance_matrix(known), covariance)
    expect_identical(plan$row_blocks, expected)
    expect_identical(plan$sampling$representation, "dense")
    expect_identical(plan$sampling$covariance, covariance)
  }

  # Default singleton publication binding must also agree for a diagonal V
  # expressed with cancelling columns, both at fitting and for explicit V_new.
  diagonal_factor <- known_v_factor(c(1, 1), rbind(c(1, 1), c(1, -1)))
  diagonal <- diag(3, 2L)
  dat <- data.frame(yi = c(.1, -.2))
  for (input in list(diagonal, diagonal_factor)) {
    fit <- construct(input, dat, group = NULL, rule = "product")
    known <- .data_known_v_data(fit$data)
    expect_length(.known_v_correlated_blocks(known), 0L)
    expect_identical(.data_selection_model(fit$data)$groups$row_blocks, list(1L, 2L))
    for (new_input in list(diagonal, diagonal_factor)) {
      known_new <- .known_v_newdata_prepare(new_input, 2L)
      expect_length(.known_v_correlated_blocks(known_new), 0L)
      expect_identical(.known_v_covariance_matrix(known_new), diagonal)
      context <- list(object = fit, same_data = FALSE, known_V_new = known_new,
        K = 2L, raw_newdata = dat, outcome_data = dat)
      expect_identical(.predict_joint_selection_groups(context), list(1L, 2L))
    }
  }
})

test_that("latent fit data distinguish independent rows from rows without loadings", {

  loading <- rbind(c(1, 1), c(2, 2), c(1, -1), c(0, 0))
  factor <- known_v_factor(rep(1, 4L), loading)
  dat <- data.frame(yi = c(.1, -.2, .3, -.4))
  fit <- brma.mv(yi = dat$yi, V = factor, data = dat, measure = "GEN",
    known_v_parameterization = "latent", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE)
  known <- .data_known_v_data(fit$data)
  blocks <- .known_v_backend_blocks(known, "latent")
  expected_zero <- setdiff(seq_len(4L), unlist(lapply(blocks, `[[`, "index")))
  expect_identical(expected_zero, 4L)
  expect_identical(.known_v_independent_indices(known), 3:4)
  graph_data <- .create_fit_data(fit$data, fit$priors)
  expect_identical(graph_data$known_v_independent_index, expected_zero)
  expect_identical(graph_data$known_v_independent_n, 1L)
  expect_identical(.known_v_rank(known), 4L)
  expect_identical(as.integer(attr(.create_fit_priors(fit$data, fit$priors)$sampling_z,
    "levels")), 4L)
  expect_identical(unlist(lapply(blocks, function(block) {
    seq.int(block$z_start, block$z_end)
  })), 1:4)
})
