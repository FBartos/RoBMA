test_that("post-fit factors preserve the selected law independently of fitted covariance syntax", {

  dat <- data.frame(yi = c(-.2, .1, .25), study = "a", esid = 1:3)
  sampling_diagonal <- c(.04, .09, .06)
  sampling_loading <- matrix(c(.1, -.2, .15), 3L, 1L)
  # A declared factor and a plain matrix that declines exact recovery exercise
  # the factor and the dense fitted syntax. They cannot be the same matrix any
  # more: a matrix an exact representation reproduces takes the factor route
  # however it is supplied, so each leg carries its own covariance and its own
  # oracle. The marginal variances agree, so both legs share their standard
  # errors and their selection thresholds.
  V_factor <- diag(sampling_diagonal) + tcrossprod(sampling_loading)
  V_dense  <- .dense_route_negative_covariance(diag(V_factor))
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
    V <- if (representation == "dense") V_dense else V_factor
    object <- bselmodel.mv(yi = yi, random = ~ 1 | study/esid, data = dat,
      V = if (representation == "dense") V_dense else
        known_v_factor(sampling_diagonal, sampling_loading),
      measure = "GEN", prior_unit_information_sd = 1, prior_bias = prior,
      effect_direction = "positive", only_priors = TRUE, silent = TRUE,
      # A three-row block on the general route is integrated by quasi-Monte
      # Carlo, so its budget and its agreement are stated on that scale; the
      # factor route stays deterministic.
      selection_control = if (representation == "dense") {
        set_selection_likelihood_control(
          relative_tolerance = 1e-4, points_per_scramble = 16384L,
          max_points_per_scramble = 65536L, scrambles = 16L)
      } else set_selection_likelihood_control(relative_tolerance = 1e-8))
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
      diag(estimate_sd[draw]^2, 3L) + matrix(study_sd[draw]^2, 3L, 3L)
    })
    packed <- t(vapply(expected_random, function(random) {
      covariance <- V + random
      covariance[lower.tri(covariance, diag = TRUE)]
    }, numeric(6L)))
    expect_equal(.selection_joint_covariance_lower(setup, 1L,
      random_factor_samples = factors), packed, tolerance = 1e-14)

    # Independent orthant integration defines the selected normalizer; it does
    # not reuse any selection-kernel calculation. The step weights expand the
    # product of one-sided weights into one trivariate probability per subset
    # of significant rows.
    cut <- sqrt(diag(V)) * stats::qnorm(.975)
    expected <- vapply(seq_along(means), function(draw) {
      covariance <- V + expected_random[[draw]]
      mean <- rep(means[draw], 3L)
      normalizer <- sum(apply(
        expand.grid(rep(list(c(FALSE, TRUE)), 3L)), 1L, function(significant) {
          .5^sum(!significant) * as.numeric(mvtnorm::pmvnorm(
            lower = ifelse(significant, cut, -Inf),
            upper = ifelse(significant, Inf, cut),
            mean  = mean, sigma = covariance,
            algorithm = mvtnorm::Miwa(steps = 4096L)
          ))
        }
      ))
      mvtnorm::dmvnorm(dat$yi, mean, covariance, log = TRUE) +
        sum(log(ifelse(dat$yi >= cut, 1, .5))) - log(normalizer)
    }, numeric(1L))
    observed <- as.numeric(.selection_joint_block_loglik_from_setup(setup))
    expect_equal(observed, expected,
                 tolerance = if (representation == "dense") 1e-3 else 1e-7)
    fits[[representation]] <- list(factors = factors, log_lik = observed)
  }
  # The random factors are a property of the random structure, which both legs
  # share; their selected laws belong to their own covariances.
  expect_equal(fits$dense$factors, fits$diagonal_factor$factors, tolerance = 1e-14)
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
