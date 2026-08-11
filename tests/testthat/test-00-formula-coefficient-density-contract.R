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
    object   = object,
    selected = selected
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


test_that("factor-level hypotheses resolve exact transformed coordinates", {

  transform <- list(
    schema_version = 1L,
    target_scale   = "original",
    source_names   = "mu_alloc__xXx__ablat[1]",
    target_names   = "mu_alloc__xXx__ablat[1]",
    matrix = matrix(
      0.25,
      nrow = 1L,
      dimnames = list(
        "mu_alloc__xXx__ablat[1]",
        "mu_alloc__xXx__ablat[1]"
      )
    ),
    source_transforms = stats::setNames(
      "identity", "mu_alloc__xXx__ablat[1]"
    ),
    output_transforms = stats::setNames(
      "identity", "mu_alloc__xXx__ablat[1]"
    )
  )
  class(transform) <- c("BayesTools_formula_coefficient_transform", "list")
  selected <- list(
    component = "mods",
    entry = list(
      role              = "formula_coefficient_group",
      formula_parameter = "mu"
    ),
    resolution = list(occurrences = data.frame(
      level          = "random",
      canonical_name = "mu_alloc__xXx__ablat[1]",
      stringsAsFactors = FALSE
    ))
  )
  point_refs <- data.frame(level = "random", value = 0)
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    .package = "BayesTools"
  )

  target <- .hypothesis_brma_formula_coefficient_level_targets(
    object     = list(
      fit = structure(list(), class = "BayesTools_fit")
    ),
    selected   = selected,
    point_refs = point_refs
  )

  expect_named(target, "random")
  expect_identical(
    target[["random"]][["target"]],
    "mu_alloc__xXx__ablat[1]"
  )
  expect_identical(target[["random"]][["route"]][["type"]], "affine")
  expect_identical(
    target[["random"]][["route"]][["weights"]],
    c("mu_alloc__xXx__ablat[1]" = 0.25)
  )
})


test_that("transformed coefficient plots use exact structural weights", {

  transform <- list(
    schema_version = 1L,
    target_scale   = "original",
    source_names   = c("mu_intercept", "mu_x"),
    target_names   = c("mu_intercept", "mu_x"),
    matrix = rbind(
      mu_intercept = c(mu_intercept = 1, mu_x = -2),
      mu_x         = c(mu_intercept = 0, mu_x = 0.5)
    ),
    source_transforms = c(mu_intercept = "identity", mu_x = "identity"),
    output_transforms = c(mu_intercept = "identity", mu_x = "identity")
  )
  class(transform) <- c("BayesTools_formula_coefficient_transform", "list")
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    .package = "BayesTools"
  )
  entry <- list(
    component         = "mods",
    role              = "fixed_coefficient",
    formula_parameter = "mu"
  )

  observed <- .plot_brma_formula_parameter_spec(
    object                    = list(
      fit = structure(list(), class = "BayesTools_fit")
    ),
    parameter                 = "mu_x",
    parameter_entry           = entry,
    standardized_coefficients = FALSE
  )
  standardized <- .plot_brma_formula_parameter_spec(
    object                    = list(
      fit = structure(list(), class = "BayesTools_fit")
    ),
    parameter                 = "mu_x",
    parameter_entry           = entry,
    standardized_coefficients = TRUE
  )

  expect_identical(observed[["type"]], "linear")
  expect_identical(observed[["weights"]], c(mu_x = 0.5))
  expect_null(standardized)
})


test_that("ordinary parameters retain legacy density-scale alignment", {

  entry <- list(
    component         = "mods",
    role              = "fixed_coefficient",
    formula_parameter = NA_character_
  )

  observed <- .plot_brma_formula_parameter_spec(
    object                    = list(),
    parameter                 = "mu",
    parameter_entry           = entry,
    standardized_coefficients = FALSE
  )

  expect_null(observed)
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
  expect_identical(result[["Alternative"]], "log_tau_intercept != 0.2")
  expect_identical(result[["Null"]], "log_tau_intercept = 0.2")
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


test_that("unit log-intercepts retain primitive qCMDE/IWMDE semantics", {

  transform <- list(
    schema_version             = 1L,
    formula_design_version     = 3L,
    parameter_registry_version = 3L,
    parameter                  = "log_tau",
    target_scale               = "original",
    source_names               = "log_tau_intercept",
    target_names               = "log_tau_intercept",
    matrix = matrix(
      1,
      nrow = 1L,
      dimnames = list("log_tau_intercept", "log_tau_intercept")
    ),
    source_transforms = c(log_tau_intercept = "log"),
    output_transforms = c(log_tau_intercept = "exp"),
    dependencies = data.frame(
      target      = "log_tau_intercept",
      source      = "log_tau_intercept",
      coefficient = 1,
      stringsAsFactors = FALSE
    ),
    sources = data.frame(),
    targets = data.frame()
  )
  class(transform) <- c(
    "BayesTools_formula_coefficient_transform",
    "list"
  )
  target <- list(
    formula_parameter = "log_tau",
    target            = "log_tau_intercept",
    target_i          = 1L,
    transform         = transform
  )
  target[["route"]] <- .hypothesis_brma_formula_transform_route(target)

  expect_identical(target[["route"]][["type"]], "identity")
  expect_identical(
    target[["route"]][["weights"]],
    c(log_tau_intercept = 1)
  )
  expect_identical(target[["route"]][["support"]], c(0, Inf))

  density <- structure(list(sentinel = TRUE), class = "prior_linear_density")
  observed_values <- numeric()
  testthat::local_mocked_bindings(
    JAGS_formula_prior_density = function(...) density,
    prior_density_ordinate = function(density, value) {
      observed_values <<- c(observed_values, value)
      list(exact = TRUE)
    },
    .package = "BayesTools"
  )
  target <- .hypothesis_brma_formula_prior_target(
    object      = list(fit = structure(list(), class = "BayesTools_fit")),
    samples     = list(),
    hypothesis  = BayesTools::hypothesis_parse("log_tau_intercept = 1.5"),
    target_info = target
  )

  expect_identical(observed_values, 1.5)
  expect_identical(target[["parameter_spec"]][["type"]], "primitive")
  expect_identical(
    target[["parameter_spec"]][["prior_density"]],
    density
  )
})


test_that("log-intercept support applies on the standardized route", {

  transform <- list(
    schema_version             = 1L,
    formula_design_version     = 3L,
    parameter_registry_version = 3L,
    parameter                  = "log_tau",
    target_scale               = "original",
    source_names               = "log_tau_intercept",
    target_names               = "log_tau_intercept",
    matrix = matrix(
      1,
      nrow = 1L,
      dimnames = list("log_tau_intercept", "log_tau_intercept")
    ),
    source_transforms = c(log_tau_intercept = "log"),
    output_transforms = c(log_tau_intercept = "exp"),
    dependencies = data.frame(),
    sources = data.frame(),
    targets = data.frame()
  )
  class(transform) <- c(
    "BayesTools_formula_coefficient_transform",
    "list"
  )
  selected <- list(
    parameter = "log_tau_intercept",
    component = "scale",
    entry     = list(formula_parameter = "log_tau")
  )
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    .package = "BayesTools"
  )
  target <- .hypothesis_brma_formula_coefficient_target(
    object   = list(fit = structure(list(), class = "BayesTools_fit")),
    selected = selected
  )
  target[["route"]] <- .hypothesis_brma_formula_transform_route(target)

  zero <- .hypothesis_brma_point_refs(
    BayesTools::hypothesis_parse("log_tau_intercept = 0"),
    "log_tau_intercept"
  )
  negative <- .hypothesis_brma_point_refs(
    BayesTools::hypothesis_parse("log_tau_intercept = -1"),
    "log_tau_intercept"
  )
  positive <- .hypothesis_brma_point_refs(
    BayesTools::hypothesis_parse("log_tau_intercept = 1"),
    "log_tau_intercept"
  )
  expect_error(
    .hypothesis_brma_check_formula_point_support(zero, target),
    "outside or on the boundary"
  )
  expect_error(
    .hypothesis_brma_check_formula_point_support(negative, target),
    "outside or on the boundary"
  )
  expect_invisible(
    .hypothesis_brma_check_formula_point_support(positive, target)
  )
})


test_that("implicit exp-affine equality preserves the Bayes factor orientation", {

  parameter <- "log_tau_intercept"
  probabilities <- (seq_len(2001L) - 0.5) / 2001
  posterior <- stats::setNames(
    data.frame(stats::qnorm(probabilities, mean = -1.2, sd = 0.3)),
    parameter
  )
  prior <- stats::setNames(
    data.frame(stats::qnorm(probabilities, mean = -1.6, sd = 0.5)),
    parameter
  )

  evaluate <- function(statement) {

    original <- BayesTools::hypothesis_parse(statement)
    transformed <- .hypothesis_brma_exp_affine_log_hypothesis(original)
    out <- BayesTools::hypothesis_BF(
      posterior      = posterior,
      prior          = prior,
      hypothesis     = transformed,
      parameter      = parameter,
      seed           = 1,
      density_method = "KDE"
    )
    .hypothesis_brma_restore_hypothesis_labels(
      out        = out,
      hypothesis = original
    )
  }

  implicit_equal <- evaluate("log_tau_intercept = 0.2")
  explicit_not_equal <- evaluate(
    "log_tau_intercept != 0.2 vs log_tau_intercept = 0.2"
  )
  explicit_equal <- evaluate(
    "log_tau_intercept = 0.2 vs log_tau_intercept != 0.2"
  )

  expect_identical(
    implicit_equal["Alternative"],
    explicit_not_equal["Alternative"]
  )
  expect_identical(implicit_equal["Null"], explicit_not_equal["Null"])
  expect_identical(
    implicit_equal[["Alternative"]],
    "log_tau_intercept != 0.2"
  )
  expect_identical(implicit_equal[["Null"]], "log_tau_intercept = 0.2")
  expect_identical(
    explicit_equal[["Alternative"]],
    "log_tau_intercept = 0.2"
  )
  expect_identical(
    explicit_equal[["Null"]],
    "log_tau_intercept != 0.2"
  )
  expect_equal(
    attr(implicit_equal, "raw_BF"),
    attr(explicit_not_equal, "raw_BF"),
    tolerance = 1e-12
  )
  expect_equal(
    attr(implicit_equal, "raw_BF") * attr(explicit_equal, "raw_BF"),
    1,
    tolerance = 1e-12
  )
})


test_that("hypothesis results restore public interaction labels", {

  probabilities <- (seq_len(2001L) - 0.5) / 2001
  posterior <- stats::setNames(
    data.frame(stats::qnorm(probabilities, mean = -0.2, sd = 0.8)),
    "mu_interaction"
  )
  prior <- stats::setNames(
    data.frame(stats::qnorm(probabilities)),
    "mu_interaction"
  )
  out <- BayesTools::hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = c("mu_interaction < 0", "mu_interaction < 1"),
    parameter  = "mu_interaction"
  )
  raw_BF <- attr(out, "raw_BF", exact = TRUE)

  display_hypothesis <- BayesTools::hypothesis_parse(c(
    "alloc:ablat[random] < 0",
    "alloc:ablat[random] < alloc:ablat[systematic]"
  ))
  out <- .hypothesis_brma_restore_hypothesis_labels(
    out             = out,
    hypothesis      = display_hypothesis,
    parameter_label = "alloc:ablat"
  )

  expect_identical(out[["Alternative"]], c(
    "alloc:ablat[random] < 0",
    "alloc:ablat[random] < alloc:ablat[systematic]"
  ))
  expect_identical(out[["Null"]], c(
    "alloc:ablat[random] >= 0",
    "alloc:ablat[random] >= alloc:ablat[systematic]"
  ))
  expect_identical(rownames(out), c("alloc:ablat", "alloc:ablat1"))
  expect_false(attr(out, "rownames", exact = TRUE))
  expect_identical(attr(out, "raw_BF", exact = TRUE), raw_BF)
  expect_identical(
    attr(out, "hypothesis_ast", exact = TRUE),
    display_hypothesis
  )
  printed <- capture.output(print(out))
  expect_true(any(grepl("alloc:ablat[random] < 0", printed, fixed = TRUE)))
  expect_false(any(grepl("alloc:ablat1", printed, fixed = TRUE)))
  expect_false(any(grepl("mu_interaction", printed, fixed = TRUE)))
})
