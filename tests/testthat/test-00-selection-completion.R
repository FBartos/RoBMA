test_that("fitting warns when all applicable selection sources are conditioned", {

  make <- function(estimate = "condition", sampling = "condition",
                   other = "condition", cluster = NULL) {

    bselmodel(yi = c(.2, -.1), sei = c(.4, .5), cluster = cluster,
      measure = "GEN", effect_direction = "positive", prior_unit_information_sd = 1,
      selection = selection_model(estimate_random_effects = estimate,
        known_sampling_variance = sampling, other_random_effects = other),
      only_priors = TRUE, silent = TRUE)
  }
  expect_no_warning(conditioned <- make())
  estimate_integrated <- make(estimate = "integrate")
  sampling_integrated <- make(sampling = "integrate")
  cluster_integrated <- make(other = "integrate", cluster = c("a", "a"))
  inapplicable_other <- make(other = "integrate")
  fixed <- bselmodel.mv(yi = c(.2, -.1), V = diag(c(.16, .25)), random = NULL,
    data = data.frame(paper = c("a", "a")), measure = "GEN",
    effect_direction = "positive", prior_unit_information_sd = 1,
    selection = selection_model(known_sampling_variance = "condition", group = "paper"),
    only_priors = TRUE, silent = TRUE)

  # Stop before backend preparation: warning coverage requires no sampling.
  testthat::local_mocked_bindings(
    .create_fit_priors = function(...) stop("Stopped before fitting.", call. = FALSE),
    .package = "RoBMA"
  )
  capture <- function(object, extend = FALSE) {

    messages <- character()
    calls <- list()
    expect_error(withCallingHandlers(.fit(object, extend = extend), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      calls <<- c(calls, list(conditionCall(w)))
      invokeRestart("muffleWarning")
    }), "Stopped before fitting.", fixed = TRUE)
    list(messages = messages, calls = calls)
  }
  expected <- paste0(
    "All applicable variation sources are set to 'condition'. ",
    "Selection weights cancel, so the model uses the ordinary Gaussian likelihood ",
    "without adjustment by the weight function."
  )
  for (object in list(conditioned, inapplicable_other, fixed)) {
    warning <- capture(object)
    expect_identical(warning$messages, expected)
    expect_identical(warning$calls, list(NULL))
  }
  for (object in list(estimate_integrated, sampling_integrated, cluster_integrated)) {
    expect_identical(capture(object)$messages, character())
  }
  conditioned$fit <- list(model = "existing fit")
  expect_identical(capture(conditioned, extend = TRUE)$messages, character())
})


test_that("complete sampling conditioning requires positive acceptance", {

  make <- function(estimate = "integrate", weights = c(1, 0), tau = .3,
                   observation_weights = NULL, sampling = "condition") {

    bselmodel(
      yi = c(.2, -.1), sei = c(.4, .5), measure = "GEN",
      effect_direction = "positive", weights = observation_weights,
      prior_effect = prior("point", list(location = 0)),
      prior_heterogeneity = prior("point", list(location = tau)),
      prior_bias = prior_weightfunction("one-sided", .025, wf_fixed(weights),
        model = selection_model(estimate_random_effects = estimate,
          known_sampling_variance = sampling)),
      only_priors = TRUE, silent = TRUE
    )
  }
  expect_s3_class(make(), "bselmodel")
  expect_s3_class(make("condition", c(1, .2)), "bselmodel")
  message <- paste0(
    "Selection with zero weights is unavailable for retained rows {1}: ",
    "positive acceptance is not guaranteed for every retained context under the supplied random-effect priors. ",
    "Use 'weights = wf_cumulative()' in 'prior_weightfunction()' to assign positive weights."
  )
  expect_error(make("condition"), message, fixed = TRUE)
  expect_error(make(tau = 0), message, fixed = TRUE)
  expect_s3_class(make(observation_weights = c(1, 1)), "bselmodel")
  expect_error(make(observation_weights = c(1, 2)), paste0(
    "Non-unit 'weights' are unavailable with 'known_sampling_variance = \"condition\"'. ",
    "Omit 'weights' or set 'known_sampling_variance = \"integrate\"' in 'selection_model()'."
  ), fixed = TRUE)
  expect_s3_class(make(observation_weights = c(1, 2), sampling = "integrate"), "bselmodel")
})


test_that("all-conditioned fitting leaves auxiliary and weight priors disconnected", {

  skip_if_not_installed("rjags")
  cases <- list(
    list(cluster = NULL, tau = 0, direction = "positive"),
    list(cluster = NULL, tau = .3, direction = "negative"),
    list(cluster = c("a", "a"), tau = .3, direction = "positive")
  )
  for (case in cases) {
    allocation <- if (!is.null(case$cluster)) {
      list(prior_heterogeneity_allocation = prior("point", list(location = .4)))
    } else list()
    object <- do.call(bselmodel, c(list(
      yi = c(.2, -.1), sei = c(.4, .5), cluster = case$cluster, measure = "GEN",
      effect_direction = case$direction, prior_unit_information_sd = 1,
      prior_effect = prior("point", list(location = .1)),
      prior_heterogeneity = prior("point", list(location = case$tau)),
      selection = selection_model(estimate_random_effects = "condition",
        other_random_effects = "condition", known_sampling_variance = "condition"),
      only_priors = TRUE, silent = TRUE
    ), allocation))
    expect_true(.selection_all_sources_conditioned(object$data))
    syntax <- .create_model_syntax(object$data, object$priors)
    expect_false(grepl("theta|sampling_z|omega|gamma|dselnorm", syntax))
    fit_priors <- .create_fit_priors(object$data, object$priors)
    expect_true(all(c("theta", "sampling_z", "omega") %in% names(fit_priors)))
    syntax <- BayesTools::JAGS_add_priors(syntax, fit_priors)
    fit_data <- .create_fit_data(object$data, object$priors)
    fit_data <- fit_data[vapply(names(fit_data), function(name) {
      grepl(name, syntax, fixed = TRUE)
    }, logical(1L))]
    connection <- textConnection(syntax)
    model <- tryCatch(rjags::jags.model(
      connection, data = fit_data, n.chains = 1L, n.adapt = 0L, quiet = TRUE,
      inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 191L)
    ), finally = close(connection))
    # Fixed mu/tau/rho leave only prior-only auxiliaries and omega. JAGS has no
    # posterior sampler for any of them because the likelihood has no path to
    # those nodes; they remain available as ordinary forward prior draws.
    expect_length(rjags::list.samplers(model), 0L)
    blocks <- .data_selection_execution_plan(object$data)$row_blocks
    covariance_names <- paste0("sel_joint_block_", seq_along(blocks), "_retained_covariance")
    samples <- as.matrix(rjags::coda.samples(model,
      c("sel_joint_mu", covariance_names, "theta", "sampling_z", "omega"),
      n.iter = 2L, progress.bar = "none"))
    expect_true(all(is.finite(samples)))
    expect_true(all(c("theta[1]", "theta[2]", "sampling_z[1]", "sampling_z[2]",
      "omega[1]", "omega[2]") %in% colnames(samples)))
    direction <- if (case$direction == "positive") 1 else -1
    expect_equal(unname(samples[, c("sel_joint_mu[1]", "sel_joint_mu[2]")]),
      matrix(direction * .1, 2L, 2L), tolerance = 0)
    covariance <- diag(c(.4, .5)^2 + case$tau^2, 2L)
    if (!is.null(case$cluster)) covariance[1L, 2L] <- covariance[2L, 1L] <- .4 * case$tau^2
    for (i in seq_along(blocks)) {
      expected <- covariance[blocks[[i]], blocks[[i]], drop = FALSE]
      expected <- expected[lower.tri(expected, diag = TRUE)]
      columns <- if (length(expected) == 1L) covariance_names[i] else {
        paste0(covariance_names[i], "[", seq_along(expected), "]")
      }
      expect_equal(unname(samples[, columns, drop = FALSE]),
        matrix(expected, 2L, length(expected), byrow = TRUE), tolerance = 1e-14)
    }
  }
})


test_that("all-conditioned fixed-effect fitting preserves known sampling covariance", {

  skip_if_not_installed("rjags")
  V <- matrix(c(.16, .03, .03, .25), 2L)
  object <- bselmodel.mv(
    yi = c(.2, -.1), V = V, random = NULL, data = data.frame(paper = c("a", "a")),
    measure = "GEN", effect_direction = "positive", prior_unit_information_sd = 1,
    prior_effect = prior("point", list(location = .1)),
    selection = selection_model(known_sampling_variance = "condition",
      weight_rule = "best", group = "paper"),
    only_priors = TRUE, silent = TRUE
  )
  expect_true(.selection_all_sources_conditioned(object$data))
  expect_length(.data_selection_model(object$data)$sources$random, 0L)
  fit_priors <- .create_fit_priors(object$data, object$priors)
  expect_false("theta" %in% names(fit_priors))
  syntax <- BayesTools::JAGS_add_priors(
    .create_model_syntax(object$data, object$priors), fit_priors)
  connection <- textConnection(syntax)
  model <- tryCatch(rjags::jags.model(
    connection, data = .create_fit_data(object$data, object$priors),
    n.chains = 1L, n.adapt = 0L, quiet = TRUE,
    inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 391L)
  ), finally = close(connection))
  expect_length(rjags::list.samplers(model), 0L)
  samples <- as.matrix(rjags::coda.samples(model,
    "sel_joint_block_1_retained_covariance", n.iter = 2L, progress.bar = "none"))
  expect_equal(unname(samples),
    matrix(V[lower.tri(V, diag = TRUE)], 2L, 3L, byrow = TRUE), tolerance = 0)
})


test_that("hard selection respects the direction of singular candidate variation", {

  make <- function(x, rule = "product", side = "one-sided") {

    dat <- data.frame(yi = c(.2, .1), study = c("a", "a"), x = x)
    bselmodel.mv(
      yi = yi, V = diag(c(.16, .25)), data = dat, random = ~ 0 + x | study,
      measure = "GEN", effect_direction = "positive",
      standardize_continuous_predictors = FALSE,
      prior_effect = prior("point", list(location = 0)),
      prior_heterogeneity = prior("point", list(location = .3)),
      prior_bias = prior_weightfunction(side, .025, wf_fixed(c(1, 0)),
        model = selection_model(other_random_effects = "integrate",
          known_sampling_variance = "condition", weight_rule = rule, group = study)),
      only_priors = TRUE, silent = TRUE
    )
  }
  # A common positive direction can reach both upper tails for every context.
  expect_s3_class(make(c(1, 1)), "bselmodel.mv")
  # Opposed effects cannot overcome an arbitrarily negative sum of retained errors.
  expect_error(make(c(1, -1)), paste0(
    "Selection with zero weights is unavailable for retained rows {1, 2}: ",
    "positive acceptance is not guaranteed for every retained context under the supplied random-effect priors. ",
    "Use 'weights = wf_cumulative()' in 'prior_weightfunction()' to assign positive weights."
  ), fixed = TRUE)
  # Either outcome suffices for best selection; both tails suffice when two-sided.
  expect_s3_class(make(c(1, -1), "best"), "bselmodel.mv")
  expect_s3_class(make(c(1, -1), side = "two-sided"), "bselmodel.mv")
})


test_that("singular observed covariance has a sampling-condition remedy", {

  expect_warning(
    expect_error(
      bselmodel.mv(
        yi = yi, V = matrix(.04, 2, 2), measure = "GEN",
        data = data.frame(yi = c(.2, .1), paper = "paper"),
        prior_unit_information_sd = 1,
        prior_bias = prior_weightfunction("one-sided", .025, wf_fixed(c(1, .2)),
          model = selection_model(known_sampling_variance = "condition",
                                  group = paper)),
        only_priors = TRUE, silent = TRUE
      ),
      paste0(
        "The sampling-conditioned selection likelihood is unavailable for retained row block(s) {1, 2}: ",
        "the full sampling and random-effect covariance is not guaranteed to be positive definite. ",
        "Supply a positive-definite 'V' or strictly positive random-effect variation covering every singular direction."
      ),
      fixed = TRUE
    ),
    paste0(
      "The 'V' argument is positive semidefinite, not positive definite, ",
      "because at least one dependency block has a rank-deficient correlation structure."
    ),
    fixed = TRUE
  )
})
