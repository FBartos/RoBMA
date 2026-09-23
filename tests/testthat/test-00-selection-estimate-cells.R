.selection_estimate_cell_fixture <- function(
    model, group_covariance = NULL, covariance_input = "factor",
    selection_control = set_selection_likelihood_control(relative_tolerance = 1e-6)) {

  dat <- data.frame(yi = c(-.4, .25), study = factor(c("s1", "s1")),
                    esid = factor(c("e1", "e2")))
  residual <- c(.64, 1)
  loading <- matrix(c(.3, .15), 2L, 1L)
  sampling_covariance <- diag(residual) + tcrossprod(loading)
  V <- if (covariance_input == "factor") known_v_factor(residual, loading) else
    sampling_covariance
  random <- if (is.null(group_covariance)) ~ 1 | study / esid else ~ 1 | esid
  object <- bselmodel.mv(
    yi = yi, V = V, random = random,
    R = if (is.null(group_covariance)) NULL else list(esid = group_covariance),
    data = dat, measure = "GEN", effect_direction = "positive",
    prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    prior_heterogeneity = BayesTools::prior_random(
      sd = BayesTools::prior("point", list(location = .25))
    ),
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .5, BayesTools::wf_fixed(c(1, .35)), model = model
    ),
    selection_control = selection_control,
    only_priors = TRUE, silent = TRUE
  )
  fit <- structure(list(), formula_design = object$formula_design,
    prior_list = c(object$formula_design$mu$prior_list,
                   .create_fit_priors(object$data, object$priors)))
  samples <- c(mu_intercept = .1)
  terms <- .formula_design_random_effects_by_mode(object$formula_design$mu, "sampled")
  for (term in terms) {
    values <- if (term$n_groups == 1L) .35 else c(-.15, .2)
    names(values) <- paste0(term$parameter_stem, "_xRE_COEFx[",
                            seq_len(term$n_groups), ",1]")
    samples <- c(samples, values)
  }
  sampling_auxiliary <- c(-.18, .12)
  if (model$known_sampling_variance == "condition") {
    structure <- .selection_sampling_structure(object$data)
    basis <- matrix(0, 2L, structure$rank)
    for (block in structure$latent_blocks) {
      basis[block$index, seq.int(block$z_start, block$z_end)] <- block$B
    }
    # Encode a fixed physical error in whichever computational basis the input
    # uses. The reference below depends only on that error and the stated V.
    z <- as.vector(t(basis) %*% solve(sampling_covariance, sampling_auxiliary))
    names(z) <- paste0("sampling_z[", seq_along(z), "]")
    samples <- c(samples, z)
  }
  samples <- matrix(samples, 1L, dimnames = list(NULL, names(samples)))
  list(object = object, fit = fit, samples = samples)
}


.selection_estimate_cell_setup <- function(fixture) {

  .log_lik_posterior_setup(
    fit = fixture$fit, posterior_samples = fixture$samples,
    data = fixture$object$data, priors = fixture$object$priors,
    unit = "estimate", data_hash = NULL
  )
}


.selection_estimate_cell_reference <- function(
    y, mean, covariance, rule, normalizer_mean = mean,
    normalizer_covariance = covariance) {

  # Independent bivariate Gaussian calculation. The only integral is P(Y1<=0,
  # Y2<=0), evaluated by conditioning Y2 on Y1 with ordinary normal probabilities.
  sd <- sqrt(diag(normalizer_covariance))
  negative <- stats::pnorm(0, normalizer_mean, sd)
  if (all(normalizer_covariance == 0)) {
    joint_negative <- list(value = as.numeric(all(normalizer_mean <= 0)), abs.error = 0)
  } else if (all(normalizer_covariance == normalizer_covariance[1L, 1L])) {
    # The integrated study intercept is one shared scalar normal variable.
    joint_negative <- list(value = stats::pnorm(min(-normalizer_mean / sd)), abs.error = 0)
  } else {
    conditional_sd <- sqrt(normalizer_covariance[2L, 2L] -
      normalizer_covariance[1L, 2L]^2 / normalizer_covariance[1L, 1L])
    joint_negative <- stats::integrate(function(x) {
      stats::dnorm(x, normalizer_mean[1L], sd[1L]) * stats::pnorm(0,
        normalizer_mean[2L] + normalizer_covariance[1L, 2L] /
          normalizer_covariance[1L, 1L] * (x - normalizer_mean[1L]), conditional_sd)
    }, -Inf, 0, rel.tol = 1e-11, abs.tol = 1e-12)
  }
  weight <- .35
  normalizer <- if (rule == "product") {
    1 - (1 - weight) * sum(negative) + (1 - weight)^2 * joint_negative$value
  } else {
    1 - (1 - weight) * joint_negative$value
  }
  observed_weight <- if (rule == "product") {
    prod(ifelse(y > 0, 1, weight))
  } else if (any(y > 0)) 1 else weight
  root <- chol(covariance)
  standardized <- forwardsolve(t(root), y - mean)
  gaussian <- -log(2 * pi) - sum(log(diag(root))) - sum(standardized^2) / 2
  list(log_density = gaussian + log(observed_weight) - log(normalizer),
       log_normalizer = log(normalizer), integration_error = joint_negative$abs.error)
}


test_that("all eight source cells preserve full-error selection and their source roles", {

  modes <- c("integrate", "condition")
  default_density <- NULL
  for (estimate in modes) for (other in modes) for (sampling in modes) {
    # Direct covariance construction from the stated toy model, not from the
    # compiled factor/covariance implementation under test.
    V <- diag(c(.64, 1)) + tcrossprod(c(.3, .15))
    C_estimate <- diag(.25^2, 2L)
    C_other <- matrix(.25^2, 2L, 2L)
    integrated <- (if (estimate == "integrate") C_estimate else 0) +
      (if (other == "integrate") C_other else 0) + matrix(0, 2L, 2L)
    retained <- c(0, 0)
    if (estimate == "condition") retained <- retained + c(-.15, .2)
    if (other == "condition") retained <- retained + .35
    if (sampling == "integrate") {
      expected_mean <- .1 + retained
      expected_covariance <- V + integrated
      candidate_mean <- expected_mean
      candidate_covariance <- expected_covariance
    } else {
      # Gaussian conditioning in physical source coordinates is independent of
      # the fitting factorization. All sources receive their covariance-owned
      # correction, even when selection preserves their population law.
      expected_mean <- rep(.1, 2L)
      expected_covariance <- V + C_estimate + C_other
      delta <- as.vector(solve(expected_covariance,
        c(-.4, .25) - .1 - c(-.18, .12) - c(-.15, .2) - .35))
      expected_sampling <- c(-.18, .12) + as.vector(V %*% delta)
      expected_estimate <- c(-.15, .2) + as.vector(C_estimate %*% delta)
      expected_other <- .35 + as.vector(C_other %*% delta)
      candidate_mean <- .1 + expected_sampling
      if (estimate == "condition") candidate_mean <- candidate_mean + expected_estimate
      if (other == "condition") candidate_mean <- candidate_mean + expected_other
      candidate_covariance <- integrated
    }

    for (rule in c("product", "best")) {
      label <- paste(estimate, other, sampling, rule, sep = "/")
      model <- BayesTools::selection_model(estimate, other, sampling, rule,
                                           group = "study")
      fixture <- .selection_estimate_cell_fixture(model)
      object <- fixture$object
      bound <- .data_selection_model(object$data)
      sources <- bound$sources$random
      source_names <- vapply(sources, `[[`, character(1L), "name")
      expected_names <- c("study", "esid_study")
      expect_setequal(source_names, expected_names)
      roles <- stats::setNames(vapply(sources, `[[`, character(1L), "role"), source_names)
      expect_identical(unname(roles[expected_names]), c("other", "estimate"), info = label)
      expected_retained <- stats::setNames(
        c(other == "condition", estimate == "condition"), expected_names
      )
      retained_sources <- unname(expected_retained[source_names])
      expect_identical(unname(vapply(sources, `[[`, logical(1L), "retained")),
                       retained_sources, info = label)
      for (mode in c("sampled", "marginalized")) {
        terms <- .formula_design_random_effects_by_mode(object$formula_design$mu, mode)
        backend_sampled <- if (sampling == "condition") rep(TRUE, length(sources)) else retained_sources
        expected <- source_names[if (mode == "sampled") backend_sampled else !backend_sampled]
        expect_setequal(vapply(terms, .random_effect_term_block_name, character(1L)),
                        expected)
      }
      expect_identical(.selection_retains_sampling(object$data),
                       sampling == "condition", info = label)
      expect_equal(.outcome_data_sei(object)^2, c(.64, 1) + c(.3, .15)^2,
                   tolerance = 1e-12, info = label)

      setup <- .selection_estimate_cell_setup(fixture)
      expect_equal(as.numeric(setup$mu), expected_mean, tolerance = 1e-12, info = label)
      plan <- .data_selection_execution_plan(object$data)
      if (sampling == "integrate") {
        covariance_samples <- .selection_joint_random_covariance_samples(setup)
        for (block in seq_along(plan$row_blocks)) {
          rows <- plan$row_blocks[[block]]
          expected <- expected_covariance[rows, rows, drop = FALSE]
          expect_equal(as.numeric(.selection_joint_covariance_lower(
            setup, block, random_covariance_samples = covariance_samples
          )), expected[lower.tri(expected, diag = TRUE)], tolerance = 1e-12, info = label)
        }
      } else {
        state <- .selection_conditioned_sampling_state(setup)
        expect_equal(as.numeric(state$e), expected_sampling, tolerance = 1e-12, info = label)
        expect_equal(as.numeric(state$baseline_mu + state$e), candidate_mean,
          tolerance = 1e-12, info = label)
        expect_equal(unname(.selection_covariance_draw(state$total_covariance, 1L)), expected_covariance,
          tolerance = 1e-12, info = label)
        expect_equal(unname(.selection_covariance_draw(state$integrated_covariance, 1L)), candidate_covariance,
          tolerance = 1e-12, info = label)
        posterior <- .selection_random_source_posterior(setup, state)
        expect_equal(as.numeric(posterior$esid_study), expected_estimate,
          tolerance = 1e-12, info = label)
        expect_equal(as.numeric(posterior$study), expected_other,
          tolerance = 1e-12, info = label)
      }
      reference <- .selection_estimate_cell_reference(
        c(-.4, .25), expected_mean, expected_covariance, rule,
        candidate_mean, candidate_covariance
      )
      expect_lt(reference$integration_error, 1e-9)
      actual <- if (sampling == "condition") as.numeric(state$log_lik) else
        as.numeric(.selection_joint_loglik_from_setup(setup))
      expect_equal(actual, reference$log_density, tolerance = 1e-6, info = label)
      if (estimate == "integrate" && other == "condition" &&
          sampling == "integrate" && rule == "product") default_density <- actual
    }
  }
  default <- BayesTools::selection_model(group = "study")
  fixture <- .selection_estimate_cell_fixture(default)
  expect_identical(default$estimate_random_effects, "integrate")
  expect_equal(as.numeric(.selection_joint_loglik_from_setup(
    .selection_estimate_cell_setup(fixture)
  )), default_density, tolerance = 1e-12)
})


test_that("an estimate-level term preserves correlated known group covariance", {

  known_R <- matrix(c(1, .4, .4, 1), 2L,
                     dimnames = list(c("e1", "e2"), c("e1", "e2")))
  # Correlated known R takes the dense route. Increase its explicit point budget,
  # preserving the same diagnostic and independent-reference accuracy criteria.
  control <- set_selection_likelihood_control(
    points_per_scramble = 8192L, relative_tolerance = 1e-6
  )
  conditioned_control <- set_selection_likelihood_control(
    points_per_scramble = 8192L, max_points_per_scramble = 262144L,
    relative_tolerance = 1e-6
  )
  V <- diag(c(.64, 1)) + tcrossprod(c(.3, .15))
  C <- unname(.25^2 * known_R)
  covariance <- V + C
  for (sampling in c("integrate", "condition")) for (rule in c("product", "best")) {
    model <- BayesTools::selection_model(known_sampling_variance = sampling,
      weight_rule = rule, group = "study")
    cell_control <- if (sampling == "condition") conditioned_control else control
    fixture <- .selection_estimate_cell_fixture(model, known_R, selection_control = cell_control)
    object <- fixture$object
    sources <- .data_selection_model(object$data)$sources$random
    expect_identical(unname(vapply(sources, `[[`, character(1L), "role")), "estimate")
    expect_false(sources[[1L]]$retained)
    setup <- .selection_estimate_cell_setup(fixture)
    plan <- .data_selection_execution_plan(object$data)
    expect_identical(plan$row_blocks, list(1:2))
    expect_equal(as.numeric(setup$mu), rep(.1, 2L), tolerance = 1e-12)
    if (sampling == "integrate") {
      candidate_mean <- rep(.1, 2L)
      candidate_covariance <- covariance
      expect_equal(as.numeric(.selection_joint_covariance_lower(
        setup, 1L, random_covariance_samples = .selection_joint_random_covariance_samples(setup)
      )), covariance[lower.tri(covariance, diag = TRUE)], tolerance = 1e-12)
    } else {
      delta <- solve(covariance, c(-.4, .25) - .1 - c(-.18, .12) - c(-.15, .2))
      candidate_mean <- .1 + c(-.18, .12) + as.vector(V %*% delta)
      candidate_covariance <- C
      state <- .selection_conditioned_sampling_state(setup)
      expect_equal(.selection_covariance_draw(state$integrated_covariance, 1L), C, tolerance = 1e-12)
      expect_equal(as.numeric(state$baseline_mu + state$e), candidate_mean,
        tolerance = 1e-12)
    }
    reference <- .selection_estimate_cell_reference(c(-.4, .25), rep(.1, 2L),
      covariance, rule, candidate_mean, candidate_covariance)
    expect_lt(reference$integration_error, 1e-9)
    actual <- if (sampling == "condition") as.numeric(state$log_lik) else
      as.numeric(.selection_joint_loglik_from_setup(setup))
    expect_equal(actual, reference$log_density, tolerance = 1e-6)
    if (sampling == "condition") {
      dense <- .selection_estimate_cell_fixture(model, known_R,
        covariance_input = "dense", selection_control = cell_control)
      expect_equal(as.numeric(.selection_joint_loglik_from_setup(
        .selection_estimate_cell_setup(dense)
      )), actual, tolerance = 1e-12)
    }
  }
})
