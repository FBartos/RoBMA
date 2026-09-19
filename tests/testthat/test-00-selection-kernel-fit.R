context("Selection kernel")
skip_on_cran()

test_that("selection constructors reject the superseded target argument", {

  expect_error(
    bselmodel(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      measure                   = "SMD",
      prior_unit_information_sd = 1,
      selection_likelihood      = "exact",
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "Unused argument in bselmodel(): 'selection_likelihood'",
    fixed = TRUE
  )
  expect_error(
    RoBMA(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      measure                   = "SMD",
      prior_unit_information_sd = 1,
      selection_likelihood      = "approximate",
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "Unused argument in RoBMA(): 'selection_likelihood'",
    fixed = TRUE
  )
})

test_that("selection model titles omit the likelihood implementation", {

  args <- list(
    yi                        = c(.1, .2),
    sei                       = c(.1, .1),
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  objects <- lapply(c("integrate", "condition"), function(mode) {

    prior_bias <- BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
      model = BayesTools::selection_model(
        other_random_effects = mode, known_sampling_variance = mode
      )
    )
    do.call(bselmodel, c(args, list(prior_bias = prior_bias)))
  })
  marginal    <- objects[[1L]]
  conditional <- objects[[2L]]

  expect_identical(
    .summary.brma_model_names(marginal),
    "Bayesian Random-Effects Selection Model (k = 2)"
  )
  expect_identical(
    .summary.brma_model_names(conditional),
    "Bayesian Random-Effects Selection Model (k = 2)"
  )
  joint_plan <- .data_selection_execution_plan(marginal[["data"]])
  joint_data <- .create_fit_data(marginal[["data"]], marginal[["priors"]])
  expect_identical(
    .data_selection_model(marginal[["data"]])[["applicability"]],
    list(estimate_random_effects = TRUE, other_random_effects = FALSE, known_sampling_variance = TRUE)
  )
  expect_identical(
    .data_selection_execution_plan(conditional[["data"]])[["exactness"]], "E0"
  )
  expect_s3_class(joint_plan, "RoBMA_selection_execution_plan")
  expect_identical(joint_plan[["exactness"]], "E0")
  expect_identical(
    joint_plan[["block_methods"]],
    rep("singleton", length(marginal[["data"]][["outcome"]][["yi"]]))
  )
  expect_length(joint_plan[["designs"]], 0L)
  expect_null(joint_data[["sel_joint_cluster_nodes"]])
  expect_identical(joint_data[["sel_joint_singleton_n"]], 2L)
  expect_equal(joint_data[["sel_joint_singleton_sampling_variance"]],
               rep(.01, 2L))
  expect_false(any(grepl("^sel_joint_block_", names(joint_data))))
})

test_that("context-only conditional selection has no residual random variance node", {

  skip_if_not_installed("rjags")
  dat <- data.frame(
    yi = c(.1, .2, -.1, .05), vi = c(.01, .02, .03, .04),
    study = c("a", "a", "b", "b")
  )
  object <- bselmodel.mv(
    yi = yi, vi = vi, random = ~ 1 | study, data = dat,
    measure = "GEN", prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_fixed(c(1, .5)),
      model = BayesTools::selection_model(group = "study")
    ),
    only_priors = TRUE, silent = TRUE
  )
  plan <- .data_selection_execution_plan(object[["data"]])
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax <- BayesTools::JAGS_add_priors(
    .create_model_syntax(object[["data"]], object[["priors"]]),
    .create_fit_priors(object[["data"]], object[["priors"]])
  )

  expect_identical(plan[["block_methods"]], rep("singleton", nrow(dat)))
  expect_null(plan[["random_covariance"]])
  expect_equal(fit_data[["sel_joint_singleton_sampling_variance"]], dat$vi,
               tolerance = 0)
  expect_match(syntax, paste0(
    "sel_joint_singleton_variance[s] = ",
    "sel_joint_singleton_sampling_variance[s]\n"
  ), fixed = TRUE)
  expect_false(grepl("sel_joint_singleton_random_variance", syntax, fixed = TRUE))

  # Hold the study contexts fixed while compiling the actual likelihood syntax.
  fit_data$mu <- c(.1, .1, -.1, -.1)
  fit_data <- fit_data[vapply(names(fit_data), function(name) {
    grepl(name, syntax, fixed = TRUE)
  }, logical(1))]
  connection <- textConnection(syntax)
  on.exit(close(connection), add = TRUE)
  model <- rjags::jags.model(connection, data = fit_data, n.chains = 1L,
                             n.adapt = 0L, quiet = TRUE)
  expect_s3_class(model, "jags")
})


test_that("selection reference weights are structural convergence parameters", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_cumulative(c(1, 1, 1))
  )

  object <- bselmodel(
    yi          = c(.1, .2, .3),
    sei         = c(.1, .1, .1),
    measure     = "SMD",
    prior_bias  = prior_bias,
    only_priors = TRUE,
    silent      = TRUE
  )

  expect_identical(
    .convergence_structural_parameters(object[["priors"]]),
    c("omega[1]", "omega[0,0.025]")
  )

  object[["priors"]][["outcome"]][["bias"]] <- NULL
  expect_length(.convergence_structural_parameters(object[["priors"]]), 0L)
})

test_that("selection model fit data and syntax use only the selected-normal kernel", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )

  object <- bselmodel(
    yi                        = c(.1, .2, .3),
    sei                       = c(.1, .1, .1),
    measure                   = "SMD",
    prior_bias                = prior_bias,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax   <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_true(all(c(
    "sel_z_lower", "sel_z_upper", "sel_joint_singleton_obs_bin", "sel_sign"
  ) %in% names(fit_data)))
  expect_false(any(grepl("sel_phack|phack_z|sel_segment|sel_kernel_mode", names(fit_data))))
  expect_match(syntax, "dselnorm_step", fixed = TRUE)
  expect_match(syntax, "sqrt(sel_joint_singleton_variance[s])", fixed = TRUE)
  expect_false(grepl("dselnorm_kernel|sel_phack|phack_z|sel_segment|sel_kernel_mode", syntax))
})

test_that("mixed normal-step bias syntax uses scalar step switch", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )

  object <- RoBMA(
    yi                        = c(.1, .2, .3),
    sei                       = c(.1, .1, .1),
    measure                   = "SMD",
    prior_bias                = prior_bias,
    prior_bias_null           = BayesTools::prior_none(),
    prior_effect              = BayesTools::prior("normal", parameters = list(0, 1)),
    prior_effect_null         = NULL,
    prior_heterogeneity       = BayesTools::prior("invgamma", parameters = list(1, .15)),
    prior_heterogeneity_null  = NULL,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax   <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_true(all(c(
    "sel_z_lower", "sel_z_upper", "sel_joint_singleton_obs_bin", "sel_sign"
  ) %in% names(fit_data)))
  expect_false(any(grepl("sel_phack|phack_z|sel_segment", names(fit_data))))
  expect_match(syntax, "sel_kernel_mode_active", fixed = TRUE)
  expect_match(syntax, "dselnorm_step_switch", fixed = TRUE)
  expect_false(grepl("dselnorm_kernel", syntax, fixed = TRUE))
})


test_that("RoBMA marginal product selection integrates cluster effects", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = .025,
    weights = BayesTools::wf_fixed(c(1, .5)),
    model = BayesTools::selection_model(
      other_random_effects = "integrate", known_sampling_variance = "integrate"
    )
  )
  args <- list(
    yi                        = c(.10, .20, .05, .15),
    sei                       = rep(.10, 4L),
    cluster                   = c("a", "a", "b", "b"),
    measure                   = "SMD",
    prior_bias                = prior_bias,
    prior_bias_null           = BayesTools::prior_none(),
    prior_effect              = BayesTools::prior("normal", parameters = list(0, 1)),
    prior_effect_null         = NULL,
    prior_heterogeneity       = BayesTools::prior("invgamma", parameters = list(1, .15)),
    prior_heterogeneity_null  = NULL,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  marginal <- do.call(RoBMA, args)

  expect_true(.is_data_joint_selection(marginal[["data"]]))
  expect_identical(
    .data_selection_model(marginal[["data"]])[c("other_random_effects", "known_sampling_variance")],
    list(other_random_effects = "integrate", known_sampling_variance = "integrate")
  )
  expect_identical(
    .data_selection_execution_plan(marginal[["data"]])[["row_blocks"]],
    list(1:2, 3:4)
  )
  joint_plan <- .data_selection_execution_plan(marginal[["data"]])
  expect_s3_class(joint_plan, "RoBMA_selection_execution_plan")
  expect_identical(joint_plan[["exactness"]], "E1")
  expect_identical(joint_plan[["block_methods"]], rep("rank_one", 2L))
  expect_identical(
    joint_plan[["quadrature"]][["orders"]],
    c(15L, 31L, 63L, 127L, 255L, 511L, 1023L)
  )
  expect_named(joint_plan[["designs"]], "factor_1")
  joint_priors <- .create_fit_priors(marginal[["data"]], marginal[["priors"]])
  joint_data   <- .create_fit_data(marginal[["data"]], marginal[["priors"]])
  joint_syntax <- .create_model_syntax(marginal[["data"]], marginal[["priors"]])
  expect_null(joint_priors[["gamma"]])
  expect_identical(grep("^sel_joint_qmc_", names(joint_data), value = TRUE),
                   "sel_joint_qmc_factor_1")
  expect_length(joint_data[["sel_joint_cluster_nodes"]], 2025L)
  expect_length(joint_data[["sel_joint_cluster_log_weights"]], 2025L)
  expect_true(all(is.finite(joint_data[["sel_joint_cluster_nodes"]])))
  expect_true(all(is.finite(
    joint_data[["sel_joint_cluster_log_weights"]]
  )))
  expect_identical(
    joint_data[["sel_joint_cluster_orders"]],
    c(15L, 31L, 63L, 127L, 255L, 511L, 1023L)
  )
  expect_false(grepl("sel_joint_qmc_2", joint_syntax, fixed = TRUE))
  expect_match(joint_syntax, "dselnorm_cluster_step", fixed = TRUE)
  expect_match(joint_syntax, "sel_kernel_mode_active", fixed = TRUE)
  expect_false(grepl("gamma[", joint_syntax, fixed = TRUE))
  expect_false(grepl("dselnorm_step_switch", joint_syntax, fixed = TRUE))

  args[["prior_bias"]] <- BayesTools::prior_weightfunction(
    side = "one-sided", steps = .025, weights = BayesTools::wf_fixed(c(1, .5))
  )
  conditional <- do.call(RoBMA, args)
  conditional_syntax <- .create_model_syntax(
    conditional[["data"]],
    conditional[["priors"]]
  )
  expect_true(.is_data_joint_selection(conditional[["data"]]))
  expect_identical(
    .data_selection_model(conditional[["data"]])[c("other_random_effects", "known_sampling_variance")],
    list(other_random_effects = "condition", known_sampling_variance = "integrate")
  )
  expect_true(.selection_retains_other_random(conditional[["data"]]))
  expect_false(.selection_retains_sampling(conditional[["data"]]))
  expect_match(conditional_syntax, "dselnorm_step_switch", fixed = TRUE)
  expect_match(conditional_syntax, "gamma[cluster[i]]", fixed = TRUE)
})


test_that("powered selection weights require integrated sampling and independent product kernels", {

  weighted_args <- list(
    yi                        = c(.10, .20, .05),
    sei                       = rep(.10, 3L),
    weights                   = c(1, .5, 1),
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  for (mode in c("condition", "integrate")) {
    weighted_args[["prior_bias"]] <- BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
      model = BayesTools::selection_model(
        other_random_effects = mode, known_sampling_variance = mode
      )
    )
    if (mode == "condition") {
      expect_error(do.call(RoBMA, weighted_args), paste0(
        "Non-unit 'weights' are unavailable with 'known_sampling_variance = \"condition\"'. ",
        "Omit 'weights' or set 'known_sampling_variance = \"integrate\"' in 'selection_model()'."
      ), fixed = TRUE)
      next
    }
    independent <- do.call(RoBMA, weighted_args)
    expect_identical(
      .data_selection_execution_plan(independent[["data"]])[["row_blocks"]],
      as.list(seq_len(3L))
    )
    expect_equal(independent[["data"]][["outcome"]][["weights"]], c(1, .5, 1))
  }
  expect_error(
    do.call(RoBMA, c(weighted_args, list(cluster = c("a", "a", "b")))),
    paste0(
      "'weights' are unavailable for jointly normalized selection vectors. ",
      "Omit 'weights'."
    ),
    fixed = TRUE
  )
  weighted_args[["prior_bias"]] <- BayesTools::prior_weightfunction(
    "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
    model = BayesTools::selection_model(weight_rule = "best")
  )
  expect_error(
    do.call(RoBMA, weighted_args),
    paste0(
      "'weights' are unavailable for jointly normalized selection vectors. ",
      "Omit 'weights'."
    ),
    fixed = TRUE
  )

  no_selection <- RoBMA(
    yi                        = c(.10, .20, .05),
    sei                       = rep(.10, 3L),
    measure                   = "SMD",
    model_type                = "PP",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  expect_false(.is_priors_weightfunction(no_selection[["priors"]]))
  expect_null(.data_selection_model(no_selection[["data"]]))
  expect_null(.data_selection_execution_plan(no_selection[["data"]]))
})


test_that("exact independent selection blocks use the scalar kernel", {

  object <- bselmodel(
    yi                        = c(.10, .20, .05),
    sei                       = rep(.10, 3L),
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])

  expect_identical(
    .data_selection_execution_plan(object[["data"]])[["row_blocks"]],
    as.list(seq_len(3L))
  )
  expect_match(syntax, "dselnorm_step_switch", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", syntax, fixed = TRUE))
  expect_false(grepl("_qmc", syntax, fixed = TRUE))
  expect_false(any(grepl("_qmc$", names(fit_data))))
})

test_that("marginal cluster log likelihood routes mixed blocks by size", {

  object <- bselmodel(
    yi                        = c(.10, .20, .05),
    sei                       = rep(.10, 3L),
    cluster                   = c("a", "a", "b"),
    measure                   = "SMD",
    selection = BayesTools::selection_model(other_random_effects = "integrate", known_sampling_variance = "integrate"),
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  expect_identical(
    .data_selection_execution_plan(object[["data"]])[["row_blocks"]],
    list(1:2, 3L)
  )
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  expect_true("sel_joint_block_1_sampling_variance" %in% names(fit_data))
  expect_false("sel_joint_block_2_sampling_variance" %in% names(fit_data))

  cluster_rows   <- NULL
  singleton_rows <- NULL
  testthat::local_mocked_bindings(
    .estimate_normal_covariance_target_location_from_setup = function(setup){
      list(y = c(.10, .20, .05), means = matrix(0, nrow = 2L, ncol = 3L))
    },
    .selection_joint_signed_context = function(setup, signed_yi){
      list(obs_bin = rep(1L, 3L))
    },
    .selection_joint_random_covariance_samples = function(setup) NULL,
    .selection_joint_singleton_variances = function(
        setup, rows, block_indices, random_covariance_samples,
        random_factor_samples, static = NULL){
      singleton_rows <<- rows
      matrix(.04, nrow = 2L, ncol = 1L)
    },
    .selection_joint_singleton_loglik_matrix = function(
        yi, means, variances, sei, selection_context){
      matrix(c(31, 32), nrow = 2L, ncol = 1L)
    },
    .selection_joint_cluster_loglik_block = function(
        yi, means, residual_sd, loading, sei, selection_context,
        execution_plan, normalizer_grid = NULL, plan_native = NULL,
        block_native = NULL){
      expect_null(normalizer_grid)
      cluster_rows <<- length(yi)
      c(11, 12)
    },
    .package = "RoBMA"
  )

  observed <- .selection_joint_block_loglik_from_setup(list(
    data          = object[["data"]],
    S             = 2L,
    selection_sei = rep(.10, 3L),
    tau_within    = matrix(.15, nrow = 2L, ncol = 3L),
    tau_between   = matrix(.20, nrow = 2L, ncol = 3L)
  ))

  expect_identical(singleton_rows, 3L)
  expect_identical(cluster_rows, 2L)
  expect_equal(observed, matrix(c(11, 12, 31, 32), nrow = 2L))
})

test_that("exact singleton variances preserve covariance representations", {

  ordinary <- bselmodel(
    yi                        = c(.10, .20),
    vi                        = c(.01, .02),
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection = BayesTools::selection_model(other_random_effects = "integrate", known_sampling_variance = "integrate"),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  tau_within <- rbind(c(.10, .20), c(.30, .40))
  ordinary_setup <- list(
    data          = ordinary[["data"]],
    S             = 2L,
    tau_within    = tau_within,
    tau_between   = matrix(0, nrow = 2L, ncol = 2L),
    is_multilevel = FALSE
  )
  sampling <- matrix(c(.01, .02, .01, .02), nrow = 2L, byrow = TRUE)
  expect_equal(
    .selection_joint_singleton_variances(
      setup         = ordinary_setup,
      rows          = 1:2,
      block_indices = 1:2
    ),
    sampling + tau_within^2,
    tolerance = 1e-15
  )

  dat <- data.frame(
    yi    = c(.10, .20),
    vi    = c(.01, .02),
    study = c("a", "b")
  )
  random <- bselmodel.mv(
    yi                        = yi,
    vi                        = vi,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
      model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "integrate", group = "study"
      )
    ),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  random_setup <- list(
    data          = random[["data"]],
    S             = 2L,
    tau_within    = matrix(0, nrow = 2L, ncol = 2L),
    tau_between   = matrix(0, nrow = 2L, ncol = 2L),
    is_multilevel = FALSE
  )
  random_covariance <- array(0, dim = c(2L, 2L, 2L))
  random_covariance[1L, , ] <- matrix(c(.04, .01, .01, .09), 2L)
  random_covariance[2L, , ] <- matrix(c(.16, .02, .02, .25), 2L)
  expect_equal(
    .selection_joint_singleton_variances(
      setup                     = random_setup,
      rows                      = 1:2,
      block_indices             = 1:2,
      random_covariance_samples = random_covariance
    ),
    sampling + rbind(c(.04, .09), c(.16, .25)),
    tolerance = 1e-15
  )

  random_factor <- list(
    diagonal = rbind(c(.01, .02), c(.03, .04)),
    loadings = list(
      array(c(.20, .30), dim = c(2L, 1L, 1L)),
      array(c(.40, .50), dim = c(2L, 1L, 1L))
    ),
    ranks = c(1L, 1L)
  )
  expect_equal(
    .selection_joint_singleton_variances(
      setup                 = random_setup,
      rows                  = 1:2,
      block_indices         = 1:2,
      random_factor_samples = random_factor
    ),
    sampling + random_factor[["diagonal"]] +
      rbind(c(.20^2, .40^2), c(.30^2, .50^2)),
    tolerance = 1e-15
  )
})

test_that("independent marginal and conditional bridge densities are identical", {

  make_object <- function(mode) {

    bselmodel(
      yi                        = c(-.20, .05, .35),
      vi                        = c(.012, .018, .025),
      measure                   = "SMD",
      prior_unit_information_sd = 1,
      selection = BayesTools::selection_model(other_random_effects = mode, known_sampling_variance = mode),
      only_priors               = TRUE,
      silent                    = TRUE
    )
  }
  evaluate <- function(object) {

    data     <- object[["data"]]
    priors   <- object[["priors"]]
    fit_data <- .create_fit_data(data, priors)
    fit_data <- .marglik_add_selection_bridge_data(
      fit_data         = fit_data,
      priors           = priors,
      effect_direction = .data_effect_direction(data),
      model_data       = data
    )
    .log_posterior(
      parameters             = list(mu = .10, tau = .15, omega = c(.60, 1)),
      data                   = fit_data,
      is_scale              = FALSE,
      is_multilevel         = FALSE,
      is_weights            = FALSE,
      is_known_v            = FALSE,
      is_PET                = FALSE,
      is_PEESE              = FALSE,
      is_weightfunction     = TRUE,
      effect_direction      = "positive",
      outcome_type          = "norm",
      model_data            = data,
      is_random             = FALSE,
      joint_selection = FALSE
    )
  }

  marginal       <- make_object("integrate")
  conditional <- make_object("condition")
  expect_identical(
    .data_selection_execution_plan(marginal[["data"]])[["exactness"]],
    "E0"
  )
  expect_equal(evaluate(marginal), evaluate(conditional), tolerance = 0)

  joint_data <- .create_fit_data(marginal[["data"]], marginal[["priors"]])
  joint_data <- .marglik_add_selection_bridge_data(
    fit_data         = joint_data,
    priors           = marginal[["priors"]],
    effect_direction = "positive",
    model_data       = marginal[["data"]]
  )
  expect_false(any(grepl("^sel_", names(joint_data))))
  first_context  <- .marglik_selection_context(list(omega = c(.60, 1)), joint_data)
  second_context <- .marglik_selection_context(list(omega = c(.40, 1)), joint_data)
  expect_identical(
    first_context[["native_cache"]],
    second_context[["native_cache"]]
  )
  expect_equal(second_context[["omega"]], matrix(c(.40, 1), nrow = 1L))
})

test_that("marginal selection bridge routes cluster plans through quadrature", {

  object <- bselmodel(
    yi                        = c(.10, .20, .05),
    sei                       = rep(.10, 3L),
    cluster                   = c("a", "a", "b"),
    measure                   = "SMD",
    selection = BayesTools::selection_model(other_random_effects = "integrate", known_sampling_variance = "integrate"),
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  bridge_data <- .marglik_add_selection_bridge_data(
    fit_data         = fit_data,
    priors           = object[["priors"]],
    effect_direction = .data_effect_direction(object[["data"]]),
    model_data       = object[["data"]]
  )

  cluster_call   <- NULL
  singleton_call <- NULL
  joint_calls    <- 0L
  testthat::local_mocked_bindings(
    .marglik_selection_context = function(parameters, data) {
      list(obs_bin = rep(1L, 3L))
    },
    .selection_joint_singleton_loglik_matrix = function(
        yi, means, variances, sei, selection_context) {
      singleton_call <<- list(
        yi        = yi,
        means     = means,
        variances = variances,
        sei       = sei
      )
      matrix(7, nrow = 1L, ncol = 1L)
    },
    .selection_joint_cluster_loglik_block = function(
        yi, means, residual_sd, loading, sei, selection_context,
        execution_plan, normalizer_grid = NULL, plan_native = NULL,
        block_native = NULL) {
      expect_null(normalizer_grid)
      cluster_call <<- list(
        yi               = yi,
        means            = means,
        residual_sd      = residual_sd,
        loading          = loading,
        sei              = sei,
        execution_plan = execution_plan
      )
      11
    },
    .selection_joint_dense_loglik_block = function(...) {
      joint_calls <<- joint_calls + 1L
      stop("The QMC bridge path must not be used for an E1 cluster plan.")
    },
    .package = "RoBMA"
  )

  mu          <- matrix(c(.01, .02, .03), nrow = 1L)
  tau_within  <- matrix(c(.15, .16, .17), nrow = 1L)
  tau_between <- matrix(c(.20, .21, .22), nrow = 1L)
  observed <- .marglik_joint_selection_log_lik(
    parameters            = list(),
    data                  = bridge_data,
    model_data            = object[["data"]],
    bridge_context        = NULL,
    covariance_plan_cache = NULL,
    mu_samples            = mu,
    tau_within_samples    = tau_within,
    tau_between_samples   = tau_between,
    is_random             = FALSE,
    is_multilevel         = TRUE,
    fixed_zero_random     = FALSE,
    K                     = 3L
  )

  expect_identical(observed, 18)
  expect_identical(joint_calls, 0L)
  expect_equal(cluster_call[["yi"]], c(.10, .20))
  expect_equal(cluster_call[["means"]], mu[, 1:2, drop = FALSE])
  expect_equal(
    cluster_call[["residual_sd"]],
    sqrt(tau_within[, 1:2, drop = FALSE]^2 + .10^2)
  )
  expect_equal(cluster_call[["loading"]], tau_between[, 1:2, drop = FALSE])
  expect_s3_class(
    cluster_call[["execution_plan"]],
    "RoBMA_selection_execution_plan"
  )
  expect_equal(singleton_call[["yi"]], .05)
  expect_equal(singleton_call[["means"]], mu[, 3L, drop = FALSE])
  expect_equal(
    singleton_call[["variances"]],
    matrix(.10^2 + .17^2 + .22^2, nrow = 1L)
  )
})

test_that("selection spec sets probability telescoping flag once", {

  yi  <- c(.1, .2, .3)
  sei <- c(.1, .1, .1)
  safe_prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = .05,
    weights = BayesTools::wf_fixed(c(1, .5))
  )
  safe_spec <- .selection_spec(
    priors           = list(outcome = list(bias = safe_prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive",
    signed_data      = TRUE
  )

  expect_true(safe_spec[["telescope_probabilities"]])
  expect_identical(safe_spec[["jags_data"]][["sel_telescope_probabilities"]], 1L)

  wrapped_safe_prior <- BayesTools::prior_bias(selection = safe_prior)
  expect_silent(
    wrapped_safe_spec <- .selection_spec(
      priors           = list(outcome = list(bias = wrapped_safe_prior)),
      yi               = yi,
      sei              = sei,
      effect_direction = "positive",
      signed_data      = TRUE
    )
  )
  expect_true(wrapped_safe_spec[["telescope_probabilities"]])
  expect_identical(wrapped_safe_spec[["jags_data"]][["sel_telescope_probabilities"]], 1L)

  wide_prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = .05,
    weights = BayesTools::wf_fixed(c(1, 200))
  )
  expect_warning(
    wide_spec <- .selection_spec(
      priors           = list(outcome = list(bias = wide_prior)),
      yi               = yi,
      sei              = sei,
      effect_direction = "positive",
      signed_data      = TRUE
    ),
    "probability telescoping disabled"
  )

  expect_false(wide_spec[["telescope_probabilities"]])
  expect_identical(wide_spec[["jags_data"]][["sel_telescope_probabilities"]], 0L)
})

test_that("branch kernel modes follow BayesTools phack and combined labels", {

  selection <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025),
    weights = BayesTools::wf_fixed(c(1, .5))
  )
  phacking <- BayesTools::prior_phacking(form = "linear")
  bias <- BayesTools::prior_mixture(list(
    BayesTools::prior_none(),
    selection,
    phacking,
    BayesTools::prior_bias(selection, phacking)
  ))
  backend <- BayesTools::selection_backend_spec(bias)

  expect_equal(
    backend[["branch_type"]],
    c("none", "weightfunction", "phack", "combined")
  )
  expect_equal(
    backend[["branch_kernel_mode"]],
    c(
      SELKERNEL_NORMAL,
      SELKERNEL_STEP,
      SELKERNEL_PHACK_POWER,
      SELKERNEL_STEP_PHACK_POWER
    )
  )
})

test_that("single two-sided bselmodel uses active full-grid omega in JAGS", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "two-sided",
    steps   = c(.05, .10),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )

  object <- bselmodel(
    yi                        = c(.24, .31, -.18, .05),
    sei                       = c(.10, .12, .09, .20),
    measure                   = "SMD",
    prior_bias                = prior_bias,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  priors   <- .create_fit_priors(object[["data"]], object[["priors"]])
  syntax   <- BayesTools::JAGS_add_priors(
    .create_model_syntax(object[["data"]], object[["priors"]]),
    priors
  )

  expect_equal(length(fit_data[["sel_z_lower"]]), 5L)
  expect_equal(length(fit_data[["sel_z_upper"]]), 5L)
  expect_match(syntax, "omega_local\\[2\\] <- 0.5")
  expect_match(syntax, "omega\\[4\\] <- omega_local\\[2\\]")
  expect_match(syntax, "omega\\[5\\] <- omega_local\\[1\\]")
  expect_false(grepl("omega\\[6\\]", syntax))
  expect_true("omega" %in% BayesTools::JAGS_to_monitor(priors))

})

test_that("JAGS permits fixed zero weights for empty p-value bins", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, 0, .5))
  )
  sei <- c(.10, .10)
  yi  <- stats::qnorm(c(.01, .20), lower.tail = FALSE) * sei

  object <- bselmodel(
    yi                        = yi,
    sei                       = sei,
    measure                   = "SMD",
    prior_bias                = prior_bias,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  priors   <- .create_fit_priors(object[["data"]], object[["priors"]])
  syntax   <- BayesTools::JAGS_add_priors(
    .create_model_syntax(object[["data"]], object[["priors"]]),
    priors
  )

  expect_equal(fit_data[["sel_joint_singleton_obs_bin"]], c(1L, 3L))
  expect_match(syntax, "omega[2] <- 0", fixed = TRUE)

  expect_false(grepl("omega[2] ~", syntax, fixed = TRUE))
})

test_that("selection p-value bins are upper-closed at step boundaries", {

  p_cuts <- c(0, .025, .50, 1)
  sei    <- rep(1, 5)
  p_val  <- c(.025, .0250001, .50, .5000001, .0249999)
  yi     <- stats::qnorm(p_val, lower.tail = FALSE) * sei

  expect_equal(
    .selection_obs_bin(yi, sei, p_cuts, sign = 1L),
    c(1L, 2L, 2L, 3L, 1L)
  )
  expect_equal(.selection_step_bin_from_z(0, p_cuts), 2L)
  expect_equal(.selection_step_bin_from_z(stats::qnorm(.025, lower.tail = FALSE), p_cuts), 1L)
})

test_that("selection omega extraction orders indexed posterior columns numerically", {

  posterior_samples <- matrix(
    seq_len(30),
    nrow = 3,
    dimnames = list(NULL, c("omega[10]", "mu", "omega[2]", "omega[1]", "omega[3]",
                            "omega[4]", "omega[5]", "omega[6]", "omega[7]", "omega[8]"))
  )
  posterior_samples <- cbind(posterior_samples, "omega[9]" = 31:33)

  selection_spec <- list(jags_omega = "omega", n_bins = 10L)
  omega          <- .extract_selection_omega_samples(posterior_samples, selection_spec)

  expect_equal(colnames(omega), paste0("omega[", 1:10, "]"))
  expect_equal(omega[, 1], posterior_samples[, "omega[1]"])
  expect_equal(omega[, 10], posterior_samples[, "omega[10]"])

  custom_samples <- posterior_samples
  colnames(custom_samples) <- sub("^omega", "custom.omega+beta", colnames(custom_samples))
  custom_samples <- cbind(
    custom_samples,
    "omega[1]" = 101:103,
    "omega[2]" = 201:203
  )
  selection_spec <- list(jags_omega = "custom.omega+beta", n_bins = 10L)
  custom_omega   <- .extract_selection_omega_samples(custom_samples, selection_spec)

  expect_equal(colnames(custom_omega), paste0("custom.omega+beta[", 1:10, "]"))
  expect_equal(custom_omega[, 2], custom_samples[, "custom.omega+beta[2]"])
  expect_false(any(custom_omega[, 1] == custom_samples[, "omega[1]"]))

  missing_custom <- custom_samples[, grepl("^omega\\[", colnames(custom_samples)), drop = FALSE]
  expect_error(
    .extract_selection_omega_samples(missing_custom, selection_spec),
    "custom.omega\\+beta"
  )
})

test_that("marginal selection constructors integrate Gaussian dependence", {

  cluster_object <- bselmodel(
    yi                        = c(.10, .20, .05, .15),
    sei                       = rep(.10, 4L),
    cluster                   = c("a", "a", "b", "b"),
    measure                   = "SMD",
    selection = BayesTools::selection_model(other_random_effects = "integrate", known_sampling_variance = "integrate"),
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  cluster_fit_priors <- .create_fit_priors(
    cluster_object[["data"]],
    cluster_object[["priors"]]
  )
  cluster_syntax <- .create_model_syntax(
    cluster_object[["data"]],
    cluster_object[["priors"]]
  )

  expect_true(.is_data_joint_selection(cluster_object[["data"]]))
  expect_identical(
    .data_selection_model(cluster_object[["data"]])[["other_random_effects"]],
    "integrate"
  )
  expect_identical(
    .data_selection_execution_plan(cluster_object[["data"]])[["row_blocks"]],
    list(1:2, 3:4)
  )
  expect_null(cluster_fit_priors[["gamma"]])
  expect_match(cluster_syntax, "dselnorm_cluster_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", cluster_syntax, fixed = TRUE))
  expect_false(grepl("gamma[", cluster_syntax, fixed = TRUE))

  data <- data.frame(study = factor(c("a", "a", "a", "b")))
  # No exact low-rank representation reproduces Markov correlations inside the
  # evidence bound, so the block declines recovery and keeps the dense
  # multivariate syntax.
  V <- .dense_route_block_diagonal(
    c(.010, .014, .011, .012), data[["study"]]
  )
  mv_object <- bselmodel.mv(
    yi                        = c(.10, .20, .05, .15),
    V                         = V,
    random                    = ~ diag(1 | study),
    data                      = data,
    measure                   = "SMD",
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
      model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "integrate", group = "study"
      )
    ),
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  random_terms <- mv_object[["formula_design"]][["mu"]][["random_effects"]]
  mv_fit_data  <- .create_fit_data(mv_object[["data"]], mv_object[["priors"]])
  mv_syntax    <- .create_model_syntax(
    mv_object[["data"]],
    mv_object[["priors"]]
  )

  expect_s3_class(mv_object, "bselmodel.mv")
  expect_true(.is_data_joint_selection(mv_object[["data"]]))
  expect_true(all(vapply(
    random_terms,
    function(term) identical(term[["compile_mode"]], "marginalized"),
    logical(1L)
  )))
  expect_match(mv_syntax, "sel_joint_random_block_1_lower", fixed = TRUE)
  expect_length(
    grep(
      "^sel_joint_block_[0-9]+_diagonal$",
      names(mv_fit_data),
      value = TRUE
    ),
    0L
  )

  conditional_prior <- BayesTools::prior_weightfunction(
    "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
    model = BayesTools::selection_model(group = "study")
  )

  conditional_auto <- bselmodel.mv(
    yi                        = c(.10, .20, .05, .15),
    V                         = V,
    random                    = ~ diag(1 | study),
    data                      = data,
    measure                   = "SMD",
    prior_bias                = conditional_prior,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  known_V <- .data_known_v_data(conditional_auto[["data"]])
  expect_identical(.known_v_requested_parameterization(known_V), "auto")
  for (backend in c("whitened", "block_mvn")) {
    forced <- bselmodel.mv(
      yi = c(.10, .20, .05, .15), V = V, random = ~ diag(1 | study), data = data,
      measure = "SMD", prior_bias = conditional_prior,
      prior_unit_information_sd = 1, known_v_parameterization = backend,
      only_priors = TRUE, silent = TRUE
    )
    expect_false(.selection_retains_sampling(forced[["data"]]))
    expect_null(.data_known_v_data(forced[["data"]])[["selection_structure"]])
    expect_equal(.selection_joint_sampling_block(
      .selection_joint_sampling_plan(forced[["data"]]), seq_len(nrow(V))
    ), V, tolerance = 0)
  }

  variance_plan <- .marglik_marginalized_variance_plan(
    conditional_auto[["data"]]
  )
  expect_length(variance_plan[["terms"]], 0L)
  expect_null(.marglik_variance_plan_node_names(variance_plan))
})


test_that("selection prediction partitions covariance without known-V metadata", {

  prediction_data <- list(outcome = list(sei = c(.10, .15, .20)))
  random_terms <- list(list(
    block_name = "study",
    group_map  = c(1L, 1L, 2L)
  ))

  expect_identical(
    .selection_joint_dependency_blocks(prediction_data, random_terms),
    list(1:2, 3L)
  )
})


test_that("marginal selection covariance batches preserve multilevel algebra", {

  object <- bselmodel(
    yi                        = c(.10, .20),
    sei                       = c(.10, .15),
    cluster                   = c("a", "a"),
    measure                   = "SMD",
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_fixed(c(1, .5)),
      model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "integrate"
      )
    ),
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  tau_within <- rbind(c(.20, .30), c(.25, .35))
  tau_between <- rbind(c(.40, .50), c(.45, .55))
  setup <- list(
    data          = object[["data"]],
    S             = 2L,
    tau_within    = tau_within,
    tau_between   = tau_between,
    is_multilevel = TRUE
  )

  evaluated_setup <- .log_lik_evaluated_setup(
    fit                  = object[["fit"]],
    data                 = object[["data"]],
    priors               = object[["priors"]],
    unit                 = "cluster",
    data_hash            = NULL,
    mu_samples           = matrix(0, nrow = 2L, ncol = 2L),
    tau_within_samples   = tau_within,
    tau_between_samples  = tau_between,
    posterior_samples    = NULL
  )
  factor_setup <- .selection_joint_factor_block_samples(
    evaluated_setup,
    block_index = 1L
  )

  expect_true(evaluated_setup[["is_multilevel"]])
  expect_identical(dim(factor_setup[["loading"]]), c(2L, 2L))

  observed <- .selection_joint_covariance_lower(setup, block_index = 1L)
  expected_covariance <- array(NA_real_, dim = c(2L, 2L, 2L))
  for (draw in seq_len(2L)) {
    expected_covariance[draw, , ] <- diag(c(.10, .15)^2) +
      diag(tau_within[draw, ]^2) +
      outer(tau_between[draw, ], tau_between[draw, ])
  }
  expected <- t(vapply(seq_len(2L), function(draw) {
    expected_covariance[draw, , ][
      cbind(c(1L, 2L, 2L), c(1L, 1L, 2L))
    ]
  }, numeric(3L)))

  expect_equal(observed, expected)

  response_setup <- .predict_joint_selection_response_setup(
    context        = list(
      object             = object,
      posterior_samples  = matrix(0, nrow = 2L, ncol = 1L),
      conditioning_depth = "marginal",
      new_data            = object[["data"]],
      known_V_new         = NULL,
      K                   = 2L,
      is_known_v          = FALSE,
      outcome_data        = object[["data"]][["outcome"]],
      random_mv           = FALSE,
      is_multilevel       = TRUE,
      same_data           = TRUE
    ),
    location_state = list(fixed_mu = matrix(0, nrow = 2L, ncol = 2L)),
    scale_state    = list(within = tau_within, between = tau_between)
  )
  for (draw in seq_len(2L)) {
    expect_equal(.block_covariance_dense(response_setup[["covariance"]], draw),
                 expected_covariance[draw, , ])
  }
})

test_that("exact bivariate selection kernel matches rectangle integration", {

  y     <- c(.30, .50)
  mu    <- c(.10, .15)
  sigma <- matrix(c(.09, .03, .03, .16), 2L, 2L)
  sei   <- c(.20, .25)
  z     <- stats::qnorm(.025, lower.tail = FALSE)
  omega <- c(.4, 1)
  design <- BayesTools::selection_qmc_design(
    dimensions = 4L,
    points     = 16384L,
    scrambles  = 8L,
    seed       = 5L
  )
  lower <- matrix(
    sigma[cbind(c(1L, 2L, 2L), c(1L, 1L, 2L))],
    nrow = 1L
  )
  actual <- .Call(
    "RoBMA_selnorm_mnorm_step_loglik_batch",
    y, matrix(mu, nrow = 1L), lower, sei, matrix(omega, nrow = 1L),
    c(z, -Inf), c(Inf, z), c(2L, 1L), 1L, TRUE, 1L,
    as.double(design), 16384L, 8L, .01, FALSE,
    0L, NULL, PACKAGE = "RoBMA"
  )
  reflected <- .Call(
    "RoBMA_selnorm_mnorm_step_loglik_batch",
    -y, matrix(-mu, nrow = 1L), lower, sei, matrix(omega, nrow = 1L),
    c(z, -Inf), c(Inf, z), c(2L, 1L), -1L, TRUE, 1L,
    as.double(design), 16384L, 8L, .01, FALSE,
    0L, NULL, PACKAGE = "RoBMA"
  )

  threshold <- z * sei
  rectangles <- list(
    list(lower = threshold, upper = c(Inf, Inf), weight = .4^2),
    list(lower = c(threshold[[1L]], -Inf),
         upper = c(Inf, threshold[[2L]]), weight = .4),
    list(lower = c(-Inf, threshold[[2L]]),
         upper = c(threshold[[1L]], Inf), weight = .4),
    list(lower = c(-Inf, -Inf), upper = threshold, weight = 1)
  )
  normalizer <- sum(vapply(rectangles, function(rectangle) {
    rectangle[["weight"]] * as.numeric(mvtnorm::pmvnorm(
      lower = rectangle[["lower"]],
      upper = rectangle[["upper"]],
      mean  = mu,
      sigma = sigma
    ))
  }, numeric(1L)))
  expected <- mvtnorm::dmvnorm(y, mean = mu, sigma = sigma, log = TRUE) +
    log(prod(omega)) - log(normalizer)

  expect_equal(actual[["log_density"]], expected, tolerance = 5e-4)
  expect_equal(reflected[["log_density"]], actual[["log_density"]],
               tolerance = 1e-12)
  expect_lt(actual[["relative_mcse"]], 5e-4)
})


test_that("dense selection supports an excluded middle interval", {

  y      <- c(.60, -.10)
  mu     <- c(.10, .15)
  sigma  <- matrix(c(.09, .03, .03, .16), 2L, 2L)
  sei    <- c(.20, .25)
  omega  <- c(1, 0, .4)
  bounds <- c(Inf, stats::qnorm(.05, lower.tail = FALSE), 0, -Inf)
  lower  <- tail(bounds, -1L)
  upper  <- head(bounds, -1L)
  # Sum the four nonzero weighted rectangles, independently of the GHK path.
  rectangles <- expand.grid(first = c(1L, 3L), second = c(1L, 3L))
  normalizer <- sum(vapply(seq_len(nrow(rectangles)), function(i) {

    bins      <- as.integer(rectangles[i, ])
    tail_sign <- ifelse(bins == 1L, 1, -1)
    prod(omega[bins]) * as.numeric(mvtnorm::pmvnorm(
      lower     = ifelse(bins == 1L, lower[bins], -upper[bins]) * sei,
      upper     = c(Inf, Inf),
      mean      = tail_sign * mu,
      sigma     = tcrossprod(tail_sign) * sigma,
      algorithm = mvtnorm::TVPACK(abseps = 1e-10)
    ))
  }, numeric(1L)))
  expected <- mvtnorm::dmvnorm(y, mu, sigma, log = TRUE) +
    log(.4) - log(normalizer)
  design <- BayesTools::selection_qmc_design(
    dimensions = 4L, points = 16384L, scrambles = 8L, seed = 5L
  )
  for (sign in c(-1L, 1L)) {
    for (telescope in c(FALSE, TRUE)) {
      actual <- .Call(
        "RoBMA_selnorm_mnorm_step_loglik_batch",
        sign * y, matrix(sign * mu, nrow = 1L),
        matrix(sigma[lower.tri(sigma, diag = TRUE)], nrow = 1L),
        sei, matrix(omega, nrow = 1L), lower, upper, c(1L, 3L), sign,
        telescope, SELKERNEL_STEP, as.double(design), 16384L, 8L, .01, FALSE,
        0L, NULL, PACKAGE = "RoBMA"
      )
      expect_equal(actual[["log_density"]], expected, tolerance = 5e-4)
      expect_lt(actual[["relative_mcse"]], 5e-4)
    }
  }
})


test_that("exact cluster reduction matches bivariate rectangle integration", {

  y           <- c(.30, .50)
  mu          <- c(.10, .15)
  residual_sd <- c(.20, .25)
  loading     <- c(.30, .35)
  sei         <- c(.20, .25)
  z           <- stats::qnorm(.025, lower.tail = FALSE)
  omega       <- c(.4, 1)
  quadrature  <- .selection_joint_cluster_quadrature_rules(
    SELNORM_CLUSTER_QUADRATURE_ORDERS
  )
  qmc <- as.double(BayesTools::selection_qmc_design(2L, 4096L, 8L, 1L))
  actual <- .Call(
    "RoBMA_selnorm_cluster_step_loglik_batch",
    y, matrix(mu, nrow = 1L), matrix(residual_sd, nrow = 1L),
    matrix(loading, nrow = 1L), sei, matrix(omega, nrow = 1L),
    c(z, -Inf), c(Inf, z), c(2L, 1L), 1L, TRUE, SELKERNEL_STEP,
    quadrature[["nodes"]], quadrature[["log_weights"]],
    as.numeric(quadrature[["orders"]]), qmc, 256L, 4096L, 8L, .005, FALSE,
    0L, PACKAGE = "RoBMA"
  )

  covariance <- diag(residual_sd^2) + outer(loading, loading)
  threshold  <- z * sei
  rectangles <- list(
    list(lower = threshold, upper = c(Inf, Inf), weight = .4^2),
    list(
      lower = c(threshold[[1L]], -Inf),
      upper = c(Inf, threshold[[2L]]),
      weight = .4
    ),
    list(
      lower = c(-Inf, threshold[[2L]]),
      upper = c(threshold[[1L]], Inf),
      weight = .4
    ),
    list(lower = c(-Inf, -Inf), upper = threshold, weight = 1)
  )
  normalizer <- sum(vapply(rectangles, function(rectangle) {
    rectangle[["weight"]] * as.numeric(mvtnorm::pmvnorm(
      lower = rectangle[["lower"]],
      upper = rectangle[["upper"]],
      mean  = mu,
      sigma = covariance
    ))
  }, numeric(1L)))
  expected <- mvtnorm::dmvnorm(
    y,
    mean  = mu,
    sigma = covariance,
    log   = TRUE
  ) + log(prod(omega)) - log(normalizer)

  expect_equal(actual[["log_density"]], expected, tolerance = 5e-8)
  expect_lt(actual[["relative_change"]], 5e-4)
})


test_that("rank-one selection quadrature preserves direction and log fallback", {

  y           <- c(.3, -.1, .5)
  mu          <- c(.1, .2, -.1)
  residual_sd <- c(.2, .3, .25)
  loading     <- c(.3, -.2, 0)
  sei         <- c(.2, .15, .25)
  z           <- stats::qnorm(.025, lower.tail = FALSE)
  omega       <- c(1, .4)
  quadrature  <- .selection_joint_cluster_quadrature_rules(
    SELNORM_CLUSTER_QUADRATURE_ORDERS
  )
  qmc <- as.double(BayesTools::selection_qmc_design(2L, 4096L, 8L, 1L))
  normalizer <- stats::integrate(function(gamma) {
    vapply(gamma, function(value) {
      probability <- stats::pnorm(
        z * sei, mu + loading * value, residual_sd, lower.tail = FALSE
      )
      stats::dnorm(value) * prod(omega[1L] * probability +
                                omega[2L] * (1 - probability))
    }, numeric(1L))
  }, -Inf, Inf, rel.tol = 1e-10)$value
  covariance <- diag(residual_sd^2) + tcrossprod(loading)
  bins <- as.integer(ifelse(y >= z * sei, 1L, 2L))
  expected <- mvtnorm::dmvnorm(y, mu, covariance, log = TRUE) +
    sum(log(omega[bins])) - log(normalizer)

  for(sign in c(-1L, 1L)){
    for(telescope in c(FALSE, TRUE)){
      actual <- .Call(
        "RoBMA_selnorm_cluster_step_loglik_batch",
        sign * y, matrix(sign * mu, 1L), matrix(residual_sd, 1L),
        matrix(sign * loading, 1L), sei, matrix(omega, 1L),
        c(z, -Inf), c(Inf, z), bins, sign, telescope, SELKERNEL_STEP,
        quadrature$nodes, quadrature$log_weights,
        as.numeric(quadrature$orders), qmc, 256L, 4096L, 8L,
        1e-8, FALSE, 0L, PACKAGE = "RoBMA"
      )
      expect_equal(actual$log_density, expected, tolerance = 1e-8)
      expect_lte(actual$relative_change, 1e-8)
    }
  }

  # Constant selection weights cancel even when their product underflows:
  # the defined likelihood is the ordinary multivariate normal density.
  y <- rep(y, 10L)
  mu <- rep(mu, 10L)
  residual_sd <- rep(residual_sd, 10L)
  loading <- rep(loading, 10L)
  actual <- .Call(
    "RoBMA_selnorm_cluster_step_loglik_batch",
    y, matrix(mu, 1L), matrix(residual_sd, 1L), matrix(loading, 1L),
    rep(sei, 10L), matrix(rep(1e-200, 2L), 1L),
    c(z, -Inf), c(Inf, z), rep(bins, 10L), 1L, TRUE, SELKERNEL_STEP,
    quadrature$nodes, quadrature$log_weights,
    as.numeric(quadrature$orders), qmc, 256L, 4096L, 8L,
    1e-8, FALSE, 0L, PACKAGE = "RoBMA"
  )
  expected <- mvtnorm::dmvnorm(
    y, mu, diag(residual_sd^2) + tcrossprod(loading), log = TRUE
  )
  expect_equal(actual$log_density, expected, tolerance = 1e-10)
})


test_that("rank-one batches keep posterior-row weights and quadrature state separate", {

  y       <- c(.3, .5)
  sei     <- c(.2, .25)
  mu      <- rbind(c(.1, .15), c(-.2, .3), c(.4, -.1))
  sd      <- rbind(c(.2, .25), c(.3, .2), c(.15, .3))
  loading <- rbind(c(.3, .35), c(0, 0), c(.5, -.2))
  omega   <- rbind(c(1, .4), c(.2, 1), c(1, .1))
  modes   <- c(SELKERNEL_STEP, SELKERNEL_NORMAL, SELKERNEL_STEP)
  z       <- stats::qnorm(.025, lower.tail = FALSE)
  bins    <- as.integer(ifelse(y >= z * sei, 1L, 2L))
  quadrature <- .selection_joint_cluster_quadrature_rules(
    SELNORM_CLUSTER_QUADRATURE_ORDERS
  )
  qmc <- as.double(BayesTools::selection_qmc_design(2L, 4096L, 8L, 1L))
  actual <- .Call(
    "RoBMA_selnorm_cluster_step_loglik_batch",
    y, mu, sd, loading, sei, omega, c(z, -Inf), c(Inf, z), bins,
    1L, TRUE, modes, quadrature$nodes, quadrature$log_weights,
    as.numeric(quadrature$orders), qmc, 256L, 4096L, 8L,
    1e-8, FALSE, 0L, PACKAGE = "RoBMA"
  )
  expected <- vapply(seq_len(nrow(mu)), function(row) {
    value <- mvtnorm::dmvnorm(
      y, mu[row, ], diag(sd[row, ]^2) + tcrossprod(loading[row, ]), log = TRUE
    )
    if (modes[row] == SELKERNEL_NORMAL) return(value)
    normalizer <- stats::integrate(function(gamma) {
      probability1 <- stats::pnorm(z * sei[1], mu[row, 1] + loading[row, 1] * gamma,
                                   sd[row, 1], lower.tail = FALSE)
      probability2 <- stats::pnorm(z * sei[2], mu[row, 2] + loading[row, 2] * gamma,
                                   sd[row, 2], lower.tail = FALSE)
      weight <- function(p) omega[row, 1] * p + omega[row, 2] * (1 - p)
      stats::dnorm(gamma) * weight(probability1) * weight(probability2)
    }, -Inf, Inf, rel.tol = 1e-10)$value
    value + sum(log(omega[row, bins])) - log(normalizer)
  }, numeric(1L))
  expect_equal(actual$log_density, expected, tolerance = 1e-8)
  expect_true(all(actual$relative_change <= 1e-8))
})


test_that("exact diagonal selection kernel reduces to analytic row factors", {

  y     <- c(.30, .50)
  mu    <- c(.10, .15)
  sd    <- c(.30, .40)
  sei   <- c(.20, .25)
  z     <- stats::qnorm(.025, lower.tail = FALSE)
  omega <- c(.4, 1)
  design <- BayesTools::selection_qmc_design(
    dimensions = 4L,
    points     = 8L,
    scrambles  = 2L,
    seed       = 9L
  )
  actual <- .Call(
    "RoBMA_selnorm_mnorm_step_loglik_batch",
    y, matrix(mu, nrow = 1L), matrix(c(sd[[1L]]^2, 0, sd[[2L]]^2), nrow = 1L),
    sei, matrix(omega, nrow = 1L), c(z, -Inf), c(Inf, z), c(2L, 1L),
    1L, TRUE, 1L, as.double(design), 8L, 2L, .01, FALSE,
    0L, NULL, PACKAGE = "RoBMA"
  )

  threshold       <- z * sei
  significant     <- y >= threshold
  observed_weight <- ifelse(significant, omega[[1L]], omega[[2L]])
  normalizer      <- omega[[1L]] * stats::pnorm(
    threshold,
    mean       = mu,
    sd         = sd,
    lower.tail = FALSE
  ) + omega[[2L]] * stats::pnorm(threshold, mean = mu, sd = sd)
  expected <- sum(
    stats::dnorm(y, mean = mu, sd = sd, log = TRUE) +
      log(observed_weight) - log(normalizer)
  )

  expect_equal(actual[["log_density"]], expected, tolerance = 1e-13)
  expect_identical(actual[["relative_mcse"]], 0)
})


test_that("exact singleton selection kernel reduces to the scalar density", {

  y     <- .30
  mu    <- .10
  sd    <- .30
  sei   <- .20
  z     <- stats::qnorm(.025, lower.tail = FALSE)
  omega <- c(.4, 1)
  design <- BayesTools::selection_qmc_design(
    dimensions = 2L,
    points     = 8L,
    scrambles  = 2L,
    seed       = 13L
  )
  actual <- .Call(
    "RoBMA_selnorm_mnorm_step_loglik_batch",
    y, matrix(mu, nrow = 1L), matrix(sd^2, nrow = 1L), sei,
    matrix(omega, nrow = 1L), c(z, -Inf), c(Inf, z), 2L,
    1L, TRUE, SELKERNEL_STEP, as.double(design),
    8L, 2L, .01, FALSE,
    0L, NULL, PACKAGE = "RoBMA"
  )
  threshold  <- z * sei
  normalizer <- omega[[1L]] * stats::pnorm(
    threshold,
    mean       = mu,
    sd         = sd,
    lower.tail = FALSE
  ) + omega[[2L]] * stats::pnorm(threshold, mean = mu, sd = sd)
  expected <- stats::dnorm(y, mean = mu, sd = sd, log = TRUE) -
    log(normalizer)

  expect_equal(actual[["log_density"]], expected, tolerance = 1e-13)
  expect_identical(actual[["relative_mcse"]], 0)
})


test_that("zero observed weights bypass unrequested exact normalization", {

  y      <- c(0, 0)
  mean   <- matrix(0, 1L, 2L)
  sei    <- c(.2, .25)
  sd     <- matrix(c(.2, .25), 1L)
  weight <- matrix(c(1, 0), 1L)
  cutoff <- stats::qnorm(.975)
  common <- list(sei, weight, c(cutoff, -Inf), c(Inf, cutoff),
                 c(2L, 2L), 1L, TRUE, SELKERNEL_STEP)
  qmc <- as.double(BayesTools::selection_qmc_design(
    dimensions = 4L, points = 8L, scrambles = 2L, seed = 1L
  ))
  cluster_qmc <- as.double(BayesTools::selection_qmc_design(
    dimensions = 2L, points = 8L, scrambles = 2L, seed = 1L
  ))
  quadrature <- .selection_joint_cluster_quadrature_rules(c(1L, 3L, 5L))
  rules <- list(quadrature$nodes, quadrature$log_weights,
                as.double(quadrature$orders))
  inputs <- list(
    mnorm = c(list(y, mean, matrix(c(.13, .08, .185), 1L)), common,
              list(qmc, 8L, 2L, 1e-12, FALSE)),
    cluster = c(list(y, mean, sd, matrix(c(.3, .35), 1L)), common,
                rules, list(cluster_qmc, 8L, 8L, 2L, 1e-12, FALSE)),
    factor = c(list(y, mean, sd, matrix(c(.3, .35, .1, -.15), 1L)),
               common, rules,
               list(qmc, 8L, 8L, 2L, 1e-12, FALSE))
  )
  for (method in names(inputs)) {
    tail_arguments <- if (method == "mnorm") list(0L, NULL) else list(0L)
    result <- do.call(.Call, c(list(
      paste0("RoBMA_selnorm_", method, "_step_loglik_batch")
    ), inputs[[method]], tail_arguments, list(PACKAGE = "RoBMA")))
    expect_identical(result$log_density, -Inf)
    expect_identical(result$log_normalizer, NA_real_)
    diagnostics <- result[intersect(c("relative_mcse", "relative_change"),
                                    names(result))]
    expect_identical(unlist(diagnostics, use.names = FALSE),
                     rep(0, length(diagnostics)))
  }
})


test_that("exact singleton blocks share the compiled scalar batch", {

  y     <- c(.30, .50, -.10)
  means <- rbind(c(.10, .15, 0), c(.15, .20, -.05))
  sd    <- rbind(c(.30, .40, .25), c(.35, .45, .30))
  sei   <- c(.20, .25, .15)
  z     <- stats::qnorm(.025, lower.tail = FALSE)
  omega <- rbind(c(.4, 1), c(.6, 1))
  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = .025,
    weights = BayesTools::wf_fixed(c(1, .4))
  )
  context <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = y,
    sei              = sei,
    effect_direction = "positive",
    signed_data      = TRUE
  )
  context[["omega"]]       <- omega
  context[["kernel_mode"]] <- rep(SELKERNEL_STEP, 2L)

  actual <- .selection_joint_singleton_loglik_matrix(
    yi                = y,
    means             = means,
    variances         = sd^2,
    sei               = sei,
    selection_context = context
  )
  expected <- matrix(NA_real_, nrow = 2L, ncol = 3L)
  for (sample_index in seq_len(2L)) {
    threshold <- z * sei
    normalizer <- omega[sample_index, 1L] * stats::pnorm(
      threshold,
      mean       = means[sample_index, ],
      sd         = sd[sample_index, ],
      lower.tail = FALSE
    ) + omega[sample_index, 2L] * stats::pnorm(
      threshold,
      mean = means[sample_index, ],
      sd   = sd[sample_index, ]
    )
    observed_weight <- omega[
      sample_index,
      ifelse(y >= threshold, 1L, 2L)
    ]
    expected[sample_index, ] <- stats::dnorm(
      y,
      mean = means[sample_index, ],
      sd   = sd[sample_index, ],
      log  = TRUE
    ) + log(observed_weight) - log(normalizer)
  }

  expect_equal(actual, expected, tolerance = 1e-13)
})


test_that("exact selected multivariate response RNG matches region masses", {

  set.seed(41)
  S     <- 6000L
  mu    <- c(.10, .15)
  sigma <- matrix(c(.09, .03, .03, .16), 2L, 2L)
  sei   <- c(.20, .25)
  omega <- c(40, 100)
  context <- list(
    omega       = matrix(omega, nrow = S, ncol = 2L, byrow = TRUE),
    kernel_mode = rep(SELKERNEL_STEP, S),
    vector_rule = 0L,
    use_normal  = rep(FALSE, S),
    p_cuts      = c(0, .025, 1),
    sign        = 1L
  )
  covariance <- array(
    rep(sigma, each = S),
    dim = c(S, 2L, 2L)
  )
  draws <- .outcome_rng.selnorm_mvn(
    mu_samples         = matrix(mu, nrow = S, ncol = 2L, byrow = TRUE),
    covariance_samples = covariance,
    sei                = sei,
    selection_context  = context,
    dependency_blocks  = list(1:2)
  )

  threshold <- stats::qnorm(.025, lower.tail = FALSE) * sei
  region <- 1L + (draws[, 1L] >= threshold[[1L]]) +
    2L * (draws[, 2L] >= threshold[[2L]])
  rectangles <- list(
    list(lower = c(-Inf, -Inf), upper = threshold, weight = 1),
    list(lower = c(threshold[[1L]], -Inf),
         upper = c(Inf, threshold[[2L]]), weight = .4),
    list(lower = c(-Inf, threshold[[2L]]),
         upper = c(threshold[[1L]], Inf), weight = .4),
    list(lower = threshold, upper = c(Inf, Inf), weight = .4^2)
  )
  expected <- vapply(rectangles, function(rectangle) {
    rectangle[["weight"]] * as.numeric(mvtnorm::pmvnorm(
      lower = rectangle[["lower"]],
      upper = rectangle[["upper"]],
      mean  = mu,
      sigma = sigma
    ))
  }, numeric(1L))
  expected <- expected / sum(expected)

  expect_equal(as.numeric(prop.table(table(factor(region, levels = 1:4)))),
               expected, tolerance = .025)
})


test_that("selected Gaussian proposals retain rare product-event moments and direction", {

  S <- 2000L
  rho <- .8
  se <- c(1, 1.5)
  sigma <- outer(se, se) * matrix(c(1, rho, rho, 1), 2L)
  cutoff <- 5
  probability <- stats::pnorm(cutoff, lower.tail = FALSE)
  # Independent one-dimensional Gaussian conditioning integral for the
  # bivariate tail event, not the rejection proposal's own normalizer.
  mass <- stats::integrate(function(z) {
    stats::dnorm(z) * stats::pnorm((cutoff - rho * z) / sqrt(1 - rho^2),
                                  lower.tail = FALSE)
  }, cutoff, Inf, rel.tol = 1e-10, abs.tol = 1e-20)$value
  first <- stats::integrate(function(z) {
    z * stats::dnorm(z) * stats::pnorm((cutoff - rho * z) / sqrt(1 - rho^2),
                                      lower.tail = FALSE)
  }, cutoff, Inf, rel.tol = 1e-10, abs.tol = 1e-20)$value / mass
  withr::local_seed(364)
  for (direction in c(1L, -1L)) {
    context <- list(
      omega = matrix(c(1e-250, 0), S, 2L, byrow = TRUE),
      kernel_mode = rep(SELKERNEL_STEP, S), vector_rule = 0L,
      use_normal = rep(FALSE, S), p_cuts = c(0, probability, 1), sign = direction
    )
    draws <- .outcome_rng.selnorm_mvn(
      matrix(0, S, 2L), array(rep(sigma, each = S), c(S, 2L, 2L)),
      se, context, list(1:2)
    )
    z <- sweep(direction * draws, 2L, se, `/`)
    expect_true(all(z >= cutoff))
    expect_equal(colMeans(z), rep(first, 2L), tolerance = .035)
  }
})


test_that("exact selected response RNG preserves singular covariance support", {

  S <- 4L
  loading <- c(.2, -.3)
  covariance_matrix <- tcrossprod(loading)
  covariance <- array(
    rep(covariance_matrix, each = S),
    dim = c(S, 2L, 2L)
  )
  context <- list(
    omega       = matrix(1, nrow = S, ncol = 1L),
    kernel_mode = rep(SELKERNEL_NORMAL, S),
    vector_rule = 0L,
    use_normal  = rep(TRUE, S),
    p_cuts      = c(0, 1),
    sign        = 1L
  )
  mean <- matrix(c(.1, -.2), nrow = S, ncol = 2L, byrow = TRUE)

  set.seed(842)
  actual <- .outcome_rng.selnorm_mvn(
    mu_samples         = mean,
    covariance_samples = covariance,
    sei                = c(.1, .1),
    selection_context  = context,
    dependency_blocks  = list(1:2)
  )
  set.seed(842)
  factor <- .covariance_sampling_factor(
    .covariance_factorization(covariance_matrix)
  )
  expected <- matrix(NA_real_, nrow = S, ncol = 2L)
  for (draw in seq_len(S)) {
    expected[draw, ] <- mean[draw, ] +
      as.vector(stats::rnorm(2L) %*% factor)
  }

  expect_equal(actual, expected, tolerance = 1e-15)

  indefinite <- array(
    rep(matrix(c(1, 2, 2, 1), nrow = 2L), each = S),
    dim = c(S, 2L, 2L)
  )
  expect_error(
    .outcome_rng.selnorm_mvn(
      mu_samples         = mean,
      covariance_samples = indefinite,
      sei                = c(.1, .1),
      selection_context  = context,
      dependency_blocks  = list(1:2)
    ),
    "Selected response covariance is not positive semidefinite\\."
  )
})


test_that("exact selected response RNG rejects asymmetric covariance", {

  context <- list(
    omega       = matrix(1, nrow = 1L, ncol = 1L),
    kernel_mode = SELKERNEL_NORMAL,
    vector_rule = 0L,
    use_normal  = TRUE,
    p_cuts      = c(0, 1),
    sign        = 1L
  )
  covariance <- array(c(1, .1, .2, 1), dim = c(1L, 2L, 2L))

  expect_error(
    .outcome_rng.selnorm_mvn(
      mu_samples         = matrix(0, nrow = 1L, ncol = 2L),
      covariance_samples = covariance,
      sei                = c(1, 1),
      selection_context  = context,
      dependency_blocks  = list(1:2)
    ),
    "Selected response covariance must be symmetric.",
    fixed = TRUE
  )
})


test_that("selected response RNG distinguishes impossible contexts from exhausted searches", {

  context <- list(
    omega       = matrix(c(0, 1), nrow = 1L),
    kernel_mode = SELKERNEL_STEP,
    vector_rule = 0L,
    use_normal  = FALSE,
    p_cuts      = c(0, .025, 1),
    sign        = 1L
  )
  expect_error(
    .outcome_rng.selnorm_mvn(
      mu_samples         = matrix(1, nrow = 1L),
      covariance_samples = array(0, dim = c(1L, 1L, 1L)),
      sei                = .1,
      selection_context  = context,
      dependency_blocks  = list(1L),
      max_attempts       = 3L
    ),
    paste0(
      "Selected response simulation is unavailable because a fully retained outcome ",
      "has zero acceptance probability. Use strictly positive selection weights ",
      "or integrate an outcome-generating source."
    ),
    fixed = TRUE
  )
  # A positive scalar selected law is sampled directly, including this tail
  # that previously exhausted a short unselected-Gaussian rejection search.
  set.seed(842)
  uniforms <- stats::runif(2L)
  cutoff <- .1 * stats::qnorm(.025, lower.tail = FALSE)
  expected <- stats::qnorm(log(uniforms[2L]) +
    stats::pnorm(cutoff, mean = 1, sd = .2, log.p = TRUE),
    mean = 1, sd = .2, log.p = TRUE)
  set.seed(842)
  scalar <- .outcome_rng.selnorm_mvn(
    mu_samples = matrix(1, nrow = 1L),
    covariance_samples = array(.04, dim = c(1L, 1L, 1L)),
    sei = .1, selection_context = context, dependency_blocks = list(1L),
    max_attempts = 3L
  )
  expect_true(is.finite(scalar[1L, 1L]))
  expect_lte(scalar[1L, 1L], cutoff)
  expect_equal(as.numeric(scalar), expected, tolerance = 1e-12)

  # Retain the exhausted-search diagnostic on an SPD, strongly negatively
  # correlated joint event. Both rows below the cutoff remain mathematically
  # possible, but cannot be reached by this reproducible three-attempt search.
  set.seed(842)
  expect_error(
    .outcome_rng.selnorm_mvn(
      mu_samples         = matrix(1, nrow = 1L, ncol = 2L),
      covariance_samples = array(c(.04, -.0396, -.0396, .04), dim = c(1L, 2L, 2L)),
      sei                = c(.1, .1),
      selection_context  = context,
      dependency_blocks  = list(1:2),
      max_attempts       = 3L
    ),
    paste0(
    "Selected response RNG was rejected by diagnostics: no ",
      "proposal was accepted in 3 attempts for a dependency block of size ",
      "2. Use 'bias_adjusted = TRUE' to draw responses before selection."
    ),
    fixed = TRUE
  )
})


test_that("fixed selection mixtures require exact branch indicators", {

  selection_spec <- list(
    fixed_omega = matrix(c(1, 0.5), ncol = 1L),
    n_bins      = 1L,
    jags_omega  = "omega"
  )
  samples <- matrix(
    c(1, 2),
    ncol     = 1L,
    dimnames = list(NULL, "bias_indicator")
  )

  expect_equal(
    .extract_selection_fixed_omega_samples(samples, selection_spec),
    matrix(c(1, 0.5), ncol = 1L,
           dimnames = list(NULL, "omega[1]"))
  )

  samples[1L, 1L] <- 1 + .Machine$double.eps
  expect_error(
    .extract_selection_fixed_omega_samples(samples, selection_spec),
    "integer-valued"
  )
})
