context("BayesTools formula coefficient density contract")

test_that("conditional-normal ordinates retain exact structural classification", {

  priors <- list(
    a = BayesTools::prior("normal", list(0, 1)),
    b = BayesTools::prior("normal", list(0, 1)),
    s = BayesTools::prior("cauchy", list(0, 1), list(0, 5))
  )
  attr(priors$b, "multiply_by") <- "s"
  density <- BayesTools:::.prior_linear_combination_density(
    priors, c(a = 1, b = -1), n_grid = 512
  )
  ordinate <- .iwmde_prior_ordinate_classifications(density, 0)[[1L]]
  expect_identical(ordinate$method, "conditional_normal_mixture")
  expect_true(ordinate$exact)
  expect_true(ordinate$eligible)
  expect_false(ordinate$provenance$integration$exact)
  expect_identical(.iwmde_validate_prior_ordinate(ordinate, 0), ordinate)
  expect_length(.iwmde_ordinate_prior_warnings("mu_intercept", list(ordinate)), 0L)
})


test_that("IWMDE classifies prior ordinates of every route by eligibility", {

  density <- BayesTools:::.prior_linear_combination_density(
    list(
      a = BayesTools::prior("cauchy", list(0, 1)),
      b = BayesTools::prior("t", list(0, 1, 3))
    ),
    c(a = 1, b = 1),
    n_grid = 256
  )
  # A route outside a fixed method vocabulary is classified by the BayesTools
  # exactness rule.
  ordinate <- .iwmde_prior_ordinate_classifications(density, 0.3)[[1L]]
  expect_identical(ordinate$method, "convolution")
  expect_true(ordinate$eligible)
  expect_true(is.na(ordinate$condition))
  expect_identical(.iwmde_validate_prior_ordinate(ordinate, 0.3), ordinate)
  expect_length(.iwmde_ordinate_prior_warnings("mu_intercept", list(ordinate)), 0L)

  # A regular ordinate that the rule refuses warns before estimation.
  testthat::local_mocked_bindings(
    prior_ordinate_status = function(prior_density, values, labels = NULL) {
      out <- eligible_ordinate_status(prior_density, values)
      out[["eligible"]]  <- FALSE
      out[["condition"]] <- "BayesTools_inexact_ordinate"
      out[["reason"]]    <- "Inexact."
      out
    },
    .package = "BayesTools"
  )
  inexact <- .iwmde_prior_ordinate_classifications(density, 0.3)[[1L]]
  expect_identical(inexact$behavior, "regular")
  expect_false(inexact$eligible)
  expect_identical(inexact$condition, "BayesTools_inexact_ordinate")
  expect_match(
    .iwmde_ordinate_prior_warnings("mu_intercept", list(inexact)),
    "could not be classified from deterministic provenance",
    fixed = TRUE
  )
})


test_that("transformed coefficient routes use exact structural weights", {

  transform <- list(
    schema_version             = 2L,
    formula_design_version     = 3L,
    parameter_map_version = 1L,
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
    targets = formula_transform_targets(
      c(mu_intercept = "affine", mu_x = "affine"),
      c(mu_intercept = "identity", mu_x = "identity")
    )
  )
  class(transform) <- c("BayesTools_formula_coefficient_transform", "list")
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    .package = "BayesTools"
  )

  route <- .brma_formula_coefficient_route(
    object   = list(fit = structure(list(sentinel = TRUE), class = "BayesTools_fit")),
    selected = list(
      parameter = "mu_intercept",
      component = "mods",
      entry     = list(formula_parameter = "mu")
    )
  )

  expect_identical(route[["type"]], "affine")
  expect_identical(route[["formula_parameter"]], "mu")
  expect_identical(route[["target"]], "mu_intercept")
  expect_identical(route[["weights"]], c(mu_intercept = 1, mu_x = -2))
  expect_identical(route[["support"]], c(-Inf, Inf))
})


test_that("transformed coefficient plots use exact structural weights", {

  transform <- list(
    schema_version = 2L,
    target_scale   = "original",
    source_names   = c("mu_intercept", "mu_x"),
    target_names   = c("mu_intercept", "mu_x"),
    matrix = rbind(
      mu_intercept = c(mu_intercept = 1, mu_x = -2),
      mu_x         = c(mu_intercept = 0, mu_x = 0.5)
    ),
    source_transforms = c(mu_intercept = "identity", mu_x = "identity"),
    output_transforms = c(mu_intercept = "identity", mu_x = "identity"),
    targets = formula_transform_targets(
      c(mu_intercept = "affine", mu_x = "affine"),
      c(mu_intercept = "identity", mu_x = "identity")
    )
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


test_that("exp-affine KDE requires continuous unconditional structure", {

  condition_event <- structure(list(
    conditional      = character(),
    conditional_rule = "AND",
    families         = list(),
    condition_key    = "<averaged>"
  ), class = "BayesTools_condition_event")
  atoms <- BayesTools::posterior_atom_attribute(
    source = "single_model_structure"
  )
  posterior <- with_draw_metadata(
    c(.15, .25, .35),
    class     = c("mixed_posteriors.simple", "mixed_posteriors"),
    condition = list(
      conditional              = character(),
      condition_key            = "<averaged>",
      resolved_condition_event = condition_event,
      averaged                 = TRUE
    ),
    atoms     = atoms
  )
  prior_density <- structure(list(
    points = data.frame(x = numeric(), p = numeric())
  ), class = c("prior_linear_density", "prior_density"))

  expect_null(.hypothesis_plan_exp_affine_certify(posterior, prior_density))
  # The marginal posterior of a certified target is BayesTools'
  # (test-02-hypothesis.R, "certified exp-affine KDE respects its open
  # support").

  refusal <- function(sample, density = prior_density) {
    out <- .hypothesis_plan_exp_affine_certify(sample, density)
    expect_identical(out[["class"]], c("RoBMA_hypothesis_target", "RoBMA_hypothesis_unavailable"))
    out[["reason"]]
  }
  atomic <- posterior
  BayesTools::posterior_metadata(atomic, "atoms") <- BayesTools::posterior_atom_attribute(
    point_masses = data.frame(x = 0, mass = 0.25)
  )
  expect_match(refusal(atomic), "posterior is atom-free", fixed = TRUE)
  # Unconditioned draws are declared by 'averaged', not by the condition key.
  keyed <- posterior
  BayesTools::posterior_metadata(keyed, "condition") <- list(
    conditional   = character(),
    condition_key = "<averaged>",
    averaged      = FALSE
  )
  expect_match(refusal(keyed), "structural evidence for an unconditional posterior", fixed = TRUE)
  unproven <- posterior
  BayesTools::posterior_metadata(unproven, "condition") <- NULL
  expect_match(refusal(unproven), "structural evidence for an unconditional posterior", fixed = TRUE)
  atomic_prior <- prior_density
  atomic_prior[["points"]] <- data.frame(x = 0, p = 0.25)
  expect_match(refusal(posterior, atomic_prior), "prior is atom-free", fixed = TRUE)
  expect_match(
    refusal(as.numeric(posterior)),
    "require a certified scalar mixed posterior",
    fixed = TRUE
  )
})


test_that("unit log-intercepts are identity maps on their positive support", {

  transform <- list(
    schema_version             = 2L,
    formula_design_version     = 3L,
    parameter_map_version = 1L,
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
    targets = formula_transform_targets(
      c(log_tau_intercept = "identity"),
      c(log_tau_intercept = "exp")
    )
  )
  class(transform) <- c(
    "BayesTools_formula_coefficient_transform",
    "list"
  )
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    .package = "BayesTools"
  )
  route <- .brma_formula_coefficient_route(
    object   = list(fit = structure(list(), class = "BayesTools_fit")),
    selected = list(
      parameter = "log_tau_intercept",
      component = "scale",
      entry     = list(formula_parameter = "log_tau")
    )
  )

  # An identity map keeps the primitive qCMDE/IWMDE target of the
  # coefficient itself (.hypothesis_plan_scalar()).
  expect_identical(route[["type"]], "identity")
  expect_identical(route[["weights"]], c(log_tau_intercept = 1))
  expect_identical(route[["support"]], c(0, Inf))
  description <- "transformed coefficient 'log_tau_intercept'"
  for (value in c(0, -1)) {
    expect_match(
      .hypothesis_plan_support_refusal(value, route[["support"]], description)[["reason"]],
      "outside or on the boundary",
      fixed = TRUE,
      info = value
    )
  }
  expect_null(.hypothesis_plan_support_refusal(1, route[["support"]], description))
})


test_that("variance point statements keep their Bayes factor orientation through the SD", {

  parameter <- "(mu) study: tau2(intercept)"
  sd        <- "(mu) study: tau(intercept)"
  probabilities <- (seq_len(2001L) - 0.5) / 2001
  posterior <- stats::setNames(
    data.frame(stats::qgamma(probabilities, shape = 6, rate = 20)),
    sd
  )
  prior <- stats::setNames(
    data.frame(stats::qgamma(probabilities, shape = 2, rate = 4)),
    sd
  )

  evaluate <- function(statement) {

    original <- BayesTools::hypothesis_rewrite(
      BayesTools::hypothesis_parse(statement),
      c(theta = parameter)
    )
    transformed <- .hypothesis_plan_sd_hypothesis(original, sd)
    expect_equal(
      transformed[["statements"]][[1L]][["left"]][["value"]],
      sqrt(original[["statements"]][[1L]][["left"]][["value"]])
    )
    # Draw-only priors have no structural prior density: BayesTools warns
    # that the point ordinate is estimated from the prior draws.
    warned <- FALSE
    out <- withCallingHandlers(
      BayesTools::hypothesis_BF(
        posterior      = posterior,
        prior          = prior,
        hypothesis     = transformed,
        parameter      = sd,
        seed           = 1,
        density_method = "KDE",
        columns        = "all"
      ),
      BayesTools_inexact_ordinate = function(warning) {
        warned <<- TRUE
        invokeRestart("muffleWarning")
      }
    )
    expect_true(warned, info = statement)
    restored <- .hypothesis_plan_sd_restore(out, original)
    # The densities are those of the variance: the SD densities over 2 * sd.
    expect_equal(
      as.numeric(restored[["prior"]]),
      as.numeric(out[["prior"]]) / (2 * 0.3),
      tolerance = 1e-12
    )
    .hypothesis_brma_restore_hypothesis_labels(
      out        = restored,
      hypothesis = original
    )
  }

  implicit_equal     <- evaluate("theta = 0.09")
  explicit_not_equal <- evaluate("theta != 0.09 vs theta = 0.09")
  explicit_equal     <- evaluate("theta = 0.09 vs theta != 0.09")

  expect_identical(
    implicit_equal["Alternative"],
    explicit_not_equal["Alternative"]
  )
  expect_identical(implicit_equal["Null"], explicit_not_equal["Null"])
  expect_identical(
    implicit_equal[["Alternative"]],
    "(mu) study: tau2(intercept) != 0.09"
  )
  expect_identical(implicit_equal[["Null"]], "(mu) study: tau2(intercept) = 0.09")
  expect_identical(
    explicit_equal[["Alternative"]],
    "(mu) study: tau2(intercept) = 0.09"
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
  expect_identical(rownames(out), c("alloc:ablat (1)", "alloc:ablat (2)"))
  expect_false(attr(out, "rownames", exact = TRUE))
  expect_identical(attr(out, "raw_BF", exact = TRUE), raw_BF)
  expect_identical(
    attr(out, "hypothesis_ast", exact = TRUE),
    display_hypothesis
  )
  printed <- capture.output(print(out))
  expect_true(any(grepl("alloc:ablat[random] < 0", printed, fixed = TRUE)))
  expect_false(any(grepl("alloc:ablat1", printed, fixed = TRUE)))
  expect_false(any(grepl("alloc:ablat (", printed, fixed = TRUE)))
  expect_false(any(grepl("mu_interaction", printed, fixed = TRUE)))
})
