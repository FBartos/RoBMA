test_that("post-fit factors preserve the selected law independently of fitted covariance syntax", {

  dat <- data.frame(yi = c(-.2, .1), study = c("a", "a"), esid = 1:2)
  sampling_diagonal <- c(.04, .09)
  sampling_loading <- matrix(c(.1, .2), 2L, 1L)
  V <- diag(sampling_diagonal) + tcrossprod(sampling_loading)
  estimate_sd <- c(.2, .4)
  study_sd <- c(.3, .5)
  means <- c(-.1, .25)
  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .5)), model = selection_model(
      other_random_effects = "integrate", known_sampling_variance = "integrate", group = study))
  fits <- list()
  # The general covariance helper is a fallback only. A supported post-fit
  # factor law must not construct an S x K x K cube, even for dense JAGS syntax.
  testthat::local_mocked_bindings(
    .selection_joint_random_covariance_samples = function(...) stop("Unexpected dense covariance cube."),
    .package = "RoBMA"
  )
  for (representation in c("dense", "diagonal_factor")) {
    object <- bselmodel.mv(yi = yi, random = ~ 1 | study/esid, data = dat,
      V = if (representation == "dense") V else known_v_factor(sampling_diagonal, sampling_loading),
      measure = "GEN", prior_unit_information_sd = 1, prior_bias = prior,
      effect_direction = "positive", only_priors = TRUE, silent = TRUE,
      selection_control = set_selection_likelihood_control(relative_tolerance = 1e-8))
    plan <- .data_selection_execution_plan(object$data)
    expect_identical(plan$random_covariance$representation, representation)
    fit <- structure(list(), formula_design = object$formula_design,
      prior_list = c(object$formula_design$mu$prior_list,
        .create_fit_priors(object$data, object$priors)))
    samples <- cbind(mu_intercept = means)
    for (source in .data_selection_model(object$data)$sources$random) {
      terms <- object$formula_design$mu$random_effects
      term <- terms[[which(vapply(terms, function(x) identical(x$block_name, source$name), logical(1L)))]]
      name <- unique(term$sd_parameter_names)
      stopifnot(length(name) == 1L)
      samples <- cbind(samples, if (source$role == "estimate") estimate_sd else study_sd)
      colnames(samples)[ncol(samples)] <- name
    }
    setup <- .estimate_likelihood_setup_from_parts(fit, object$data, object$priors, samples)
    factors <- .selection_joint_random_factor_samples(setup)
    expect_false(is.null(factors))
    expect_identical(factors$ranks, 1L)
    expected_random <- lapply(seq_along(means), function(draw) {
      diag(estimate_sd[draw]^2, 2L) + matrix(study_sd[draw]^2, 2L, 2L)
    })
    packed <- t(vapply(expected_random, function(random) {
      covariance <- V + random
      covariance[lower.tri(covariance, diag = TRUE)]
    }, numeric(3L)))
    expect_equal(.selection_joint_covariance_lower(setup, 1L,
      random_factor_samples = factors), packed, tolerance = 1e-14)

    # Independent bivariate conditional-normal integration defines the
    # selected normalizer; it does not reuse any selection-kernel calculation.
    expected <- vapply(seq_along(means), function(draw) {
      covariance <- V + expected_random[[draw]]
      mean <- rep(means[draw], 2L)
      sd <- sqrt(diag(covariance))
      cut <- sqrt(diag(V)) * stats::qnorm(.975)
      both <- stats::integrate(function(x) {
        stats::dnorm(x, mean[1L], sd[1L]) * stats::pnorm(cut[2L],
          mean[2L] + covariance[2L, 1L] / covariance[1L, 1L] * (x - mean[1L]),
          sqrt(covariance[2L, 2L] - covariance[2L, 1L]^2 / covariance[1L, 1L]),
          lower.tail = FALSE)
      }, cut[1L], Inf, rel.tol = 1e-11, abs.tol = 1e-13)$value
      normalizer <- .25 + .25 * sum(stats::pnorm(cut, mean, sd, lower.tail = FALSE)) + .25 * both
      mvtnorm::dmvnorm(dat$yi, mean, covariance, log = TRUE) +
        sum(log(ifelse(dat$yi >= cut, 1, .5))) - log(normalizer)
    }, numeric(1L))
    observed <- as.numeric(.selection_joint_block_loglik_from_setup(setup))
    expect_equal(observed, expected, tolerance = 1e-7)
    fits[[representation]] <- list(factors = factors, log_lik = observed)
  }
  expect_equal(fits$dense$factors, fits$diagonal_factor$factors, tolerance = 1e-14)
  expect_equal(fits$dense$log_lik, fits$diagonal_factor$log_lik, tolerance = 1e-7)
})

test_that("only structural factor unavailability uses the generic fallback", {

  plan <- structure(list(schema_version = 5L,
    random_covariance = list(representation = "dense", term_names = "study"),
    row_blocks = list(1:2)), class = c("RoBMA_selection_execution_plan", "list"))
  data <- structure(list(), random = TRUE, selection_execution_plan = plan)
  setup <- list(data = data, S = 1L, K = 2L, posterior_samples = matrix(.2))
  testthat::local_mocked_bindings(
    .brma_mv_random_effects_marginal_factor_states = function(...) list(),
    .package = "RoBMA"
  )
  structural <- TRUE
  testthat::local_mocked_bindings(
    random_effects_marginal_diagonal_factor = function(...) {
      if (!structural) stop("invalid numerical factor state")
      stop(structure(list(message = "unsupported factor structure", call = NULL),
        class = c("BayesTools_random_effects_marginal_factor_unavailable", "error", "condition")))
    },
    .package = "BayesTools"
  )
  expect_null(.selection_joint_random_factor_samples(setup))
  structural <- FALSE
  expect_error(.selection_joint_random_factor_samples(setup),
    "invalid numerical factor state", fixed = TRUE)
})
