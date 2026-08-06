context("BayesTools formula coefficient density contract")


test_that("transformed coefficient hypotheses use exact structural weights", {

  transform <- list(
    schema_version             = 1L,
    formula_design_version     = 3L,
    parameter_registry_version = 3L,
    parameter                  = "mu",
    target_scale               = "original",
    source_names               = c("mu_intercept", "mu_x"),
    target_names               = c("mu_intercept", "mu_x"),
    matrix = rbind(
      mu_intercept = c(mu_intercept = 1, mu_x = -2),
      mu_x         = c(mu_intercept = 0, mu_x = 0.5)
    ),
    source_transforms = c(mu_intercept = "identity", mu_x = "identity"),
    output_transforms = c(mu_intercept = "identity", mu_x = "identity"),
    dependencies = data.frame(
      target      = c("mu_intercept", "mu_intercept", "mu_x"),
      source      = c("mu_intercept", "mu_x", "mu_x"),
      coefficient = c(1, -2, 0.5),
      stringsAsFactors = FALSE
    ),
    sources = data.frame(),
    targets = data.frame(
      target            = c("mu_intercept", "mu_x"),
      structural_status = "dependent",
      fixed_value       = NA_real_,
      reason            = "",
      stringsAsFactors  = FALSE
    )
  )
  class(transform) <- c("BayesTools_formula_coefficient_transform", "list")
  density <- structure(list(sentinel = TRUE), class = "prior_linear_density")
  fit <- structure(list(sentinel = TRUE), class = "BayesTools_fit")
  object <- list(fit = fit)
  selected <- list(
    parameter = "mu_intercept",
    component = "mods",
    entry     = list(formula_parameter = "mu")
  )
  context <- structure(list(sentinel = TRUE), class = "prior_density_context")
  samples <- structure(list(), prior_density_context = context)
  observed_context <- NULL
  observed_values  <- numeric()
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    JAGS_formula_prior_density = function(..., context) {
      observed_context <<- context
      density
    },
    prior_density_ordinate = function(density, value) {
      observed_values <<- c(observed_values, value)
      list(behavior = "regular", exact = TRUE)
    },
    .package = "BayesTools"
  )

  target <- .hypothesis_brma_formula_coefficient_target(
    object                    = object,
    selected                  = selected,
    standardized_coefficients = FALSE
  )
  target <- .hypothesis_brma_formula_prior_target(
    object      = object,
    samples     = samples,
    hypothesis  = BayesTools::hypothesis_parse("mu_intercept = 3"),
    target_info = target
  )

  expect_identical(observed_context, context)
  expect_identical(observed_values, 3)
  expect_identical(target[["prior_density"]], density)
  expect_identical(target[["parameter_spec"]][["type"]], "linear")
  expect_identical(
    target[["parameter_spec"]][["weights"]],
    c(mu_intercept = 1, mu_x = -2)
  )
})


test_that("nonlinear joint coefficient transforms fail qCMDE/IWMDE closed", {

  target <- list(
    formula_parameter = "log_tau",
    target            = "log_tau_intercept",
    target_i          = 1L,
    transform = list(
      target_scale = "original",
      matrix = matrix(
        c(1, -2),
        nrow = 1L,
        dimnames = list(
          "log_tau_intercept",
          c("log_tau_intercept", "log_tau_x")
        )
      ),
      source_transforms = c(
        log_tau_intercept = "log",
        log_tau_x         = "identity"
      ),
      output_transforms = c(log_tau_intercept = "exp")
    )
  )
  class(target[["transform"]]) <- c(
    "BayesTools_formula_coefficient_transform",
    "list"
  )
  density <- structure(list(sentinel = TRUE), class = "prior_linear_density")
  testthat::local_mocked_bindings(
    JAGS_formula_prior_density = function(...) density,
    prior_density_ordinate = function(...) list(exact = TRUE),
    .package = "BayesTools"
  )

  observed <- .hypothesis_brma_formula_prior_target(
    object      = list(fit = structure(list(), class = "BayesTools_fit")),
    samples     = list(),
    hypothesis  = BayesTools::hypothesis_parse("log_tau_intercept = 1"),
    target_info = target
  )

  expect_identical(
    observed[["parameter_spec"]][["type"]],
    "unsupported_formula_transform"
  )
  expect_match(observed[["parameter_spec"]][["reason"]], "nonlinear joint")
})


test_that("exp-affine KDE requires continuous unconditional structure", {

  transform <- list(
    target_scale = "original",
    matrix = matrix(
      c(1, -2),
      nrow = 1L,
      dimnames = list(
        "log_tau_intercept",
        c("log_tau_intercept", "log_tau_x")
      )
    ),
    source_transforms = c(
      log_tau_intercept = "log",
      log_tau_x         = "identity"
    ),
    output_transforms = c(log_tau_intercept = "exp")
  )
  class(transform) <- c(
    "BayesTools_formula_coefficient_transform",
    "list"
  )
  target <- list(
    target    = "log_tau_intercept",
    target_i  = 1L,
    transform = transform
  )
  target[["route"]] <- .hypothesis_brma_formula_transform_route(target)
  condition_event <- structure(list(
    conditional      = character(),
    conditional_rule = "AND",
    families         = list(),
    condition_key    = "<averaged>"
  ), class = "BayesTools_condition_event")
  atoms <- structure(list(
    declared               = TRUE,
    locations              = matrix(
      numeric(), nrow = 0L, ncol = 1L,
      dimnames = list(NULL, "log_tau_intercept")
    ),
    mass                   = numeric(),
    source                 = "single_model_structure",
    component_probabilities = 1
  ), class = c("BayesTools_posterior_atoms", "list"))
  posterior <- structure(
    c(.15, .25, .35),
    class                    = c("mixed_posteriors.simple", "mixed_posteriors"),
    conditional              = character(),
    condition_key            = "<averaged>",
    resolved_condition_event = condition_event,
    posterior_atoms          = atoms
  )
  prior_density <- structure(list(
    points = data.frame(x = numeric(), p = numeric())
  ), class = c("prior_linear_density", "prior_density"))
  samples <- structure(
    list(log_tau_intercept = posterior),
    prior_densities = list(log_tau_intercept = prior_density)
  )
  captured <- NULL
  testthat::local_mocked_bindings(
    transform_prior_samples = function(...) {
      cbind(log_tau_intercept = c(.1, .2, .4))
    },
    hypothesis_BF = function(posterior, prior, hypothesis, ...) {
      captured <<- list(
        posterior  = posterior,
        prior      = prior,
        hypothesis = hypothesis
      )
      parsed <- list(list(
        input = "log_tau_intercept = -1.6094379124341003",
        left  = list(
          type = "point", label = "log_tau_intercept = -1.6094379124341003",
          expr = as.name("log_tau_intercept"),
          expression = as.name("log_tau_intercept"),
          value = log(.2)
        ),
        right = list(
          type = "not_point", label = "log_tau_intercept != -1.6094379124341003",
          expr = as.name("log_tau_intercept"),
          expression = as.name("log_tau_intercept"),
          value = log(.2)
        ),
        explicit = FALSE
      ))
      structure(data.frame(
        Alternative = "log_tau_intercept = -1.6094379124341003",
        Null        = "log_tau_intercept != -1.6094379124341003"
      ), hypothesis_ast = hypothesis, parsed = parsed)
    },
    .package = "BayesTools"
  )

  result <- .hypothesis_brma_exp_affine_kde(
    object      = list(fit = structure(list(), class = "BayesTools_fit")),
    samples     = samples,
    hypothesis  = BayesTools::hypothesis_parse("log_tau_intercept = 0.2"),
    parameter   = "log_tau_intercept",
    target_info = target,
    conditional = FALSE,
    logBF       = FALSE,
    BF01        = FALSE,
    seed        = 1,
    n_samples   = 3L,
    columns     = "default"
  )

  expect_equal(
    captured[["posterior"]][["log_tau_intercept"]],
    log(c(.15, .25, .35))
  )
  expect_equal(
    captured[["prior"]][["log_tau_intercept"]],
    log(c(.1, .2, .4))
  )
  expect_equal(
    captured[["hypothesis"]][["statements"]][[1L]][["left"]][["value"]],
    log(.2)
  )
  expect_identical(result[["Alternative"]], "log_tau_intercept = 0.2")
  expect_identical(result[["Null"]], "log_tau_intercept != 0.2")
  expect_identical(attr(result, "hypothesis_ast"),
                   BayesTools::hypothesis_parse("log_tau_intercept = 0.2"))
  expect_equal(attr(result, "parsed")[[1L]][["left"]][["value"]], .2)
  expect_identical(
    attr(result, "parsed")[[1L]][["left"]][["label"]],
    "log_tau_intercept = 0.2"
  )
  expect_identical(target[["route"]][["type"]], "exp_affine")
  expect_false(exists(
    ".hypothesis_brma_coherent_draws",
    envir    = asNamespace("RoBMA"),
    inherits = FALSE
  ))

  refs_zero <- .hypothesis_brma_point_refs(
    BayesTools::hypothesis_parse("log_tau_intercept = 0"),
    "log_tau_intercept"
  )
  refs_one <- .hypothesis_brma_point_refs(
    BayesTools::hypothesis_parse("log_tau_intercept = 1"),
    "log_tau_intercept"
  )
  expect_error(
    .hypothesis_brma_check_formula_point_support(refs_zero, target),
    "outside or on the boundary"
  )
  expect_invisible(
    .hypothesis_brma_check_formula_point_support(refs_one, target)
  )

  atomic_samples <- samples
  atomic_atoms <- atoms
  atomic_atoms[["locations"]] <- matrix(
    0, nrow = 1L, ncol = 1L,
    dimnames = list("location", "log_tau_intercept")
  )
  atomic_atoms[["mass"]] <- 0.25
  attr(atomic_samples[["log_tau_intercept"]], "posterior_atoms") <-
    atomic_atoms
  expect_error(
    .hypothesis_brma_exp_affine_certify(
      atomic_samples, "log_tau_intercept", FALSE
    ),
    "posterior is atom-free"
  )
  unproven_samples <- samples
  attr(unproven_samples[["log_tau_intercept"]], "conditional") <- NULL
  attr(
    unproven_samples[["log_tau_intercept"]],
    "resolved_condition_event"
  ) <- NULL
  expect_error(
    .hypothesis_brma_exp_affine_certify(
      unproven_samples, "log_tau_intercept", FALSE
    ),
    "structural evidence for an unconditional posterior"
  )
  atomic_prior_samples <- samples
  attr(atomic_prior_samples, "prior_densities")[[
    "log_tau_intercept"
  ]][["points"]] <- data.frame(x = 0, p = 0.25)
  expect_error(
    .hypothesis_brma_exp_affine_certify(
      atomic_prior_samples, "log_tau_intercept", FALSE
    ),
    "prior is atom-free"
  )
  expect_error(
    .hypothesis_brma_exp_affine_certify(
      samples, "log_tau_intercept", TRUE
    ),
    "conditional product-space"
  )
  expect_error(
    .hypothesis_brma_exp_affine_route_kind(
      BayesTools::hypothesis_parse(c(
        "log_tau_intercept = 0.2",
        "log_tau_intercept > 0.2"
      ))
    ),
    "cannot mix point and region"
  )
})
