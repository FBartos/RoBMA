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

  testthat::local_mocked_bindings(
    JAGS_formula_prior_density = function(...) density,
    .package = "BayesTools"
  )
  target <- .hypothesis_brma_formula_prior_target(
    object = list(fit = list()), samples = list(),
    hypothesis = BayesTools::hypothesis_parse("mu_intercept = 0"),
    point_values = 0,
    target_info = list(
      formula_parameter = "mu", target = "mu_intercept",
      route = list(type = "affine", weights = c(a = 1, b = -1))
    )
  )
  expect_identical(target$prior_density, density)
  expect_identical(target$parameter_spec$type, "linear")
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


test_that("transformed coefficient hypotheses use exact structural weights", {

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
  density <- structure(list(sentinel = TRUE), class = "prior_linear_density")
  fit <- structure(list(sentinel = TRUE), class = "BayesTools_fit")
  object <- list(fit = fit)
  selected <- list(
    parameter = "mu_intercept",
    component = "mods",
    entry     = list(formula_parameter = "mu")
  )
  context <- structure(list(sentinel = TRUE), class = "prior_density_context")
  samples <- with_draw_metadata(list(), prior_context = context)
  observed_context <- NULL
  observed_values  <- numeric()
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    JAGS_formula_prior_density = function(..., context) {
      observed_context <<- context
      density
    },
    prior_ordinate_status = function(prior_density, values, labels = NULL) {
      observed_values <<- c(observed_values, values)
      eligible_ordinate_status(prior_density, values)
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
    schema_version = 2L,
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
    ),
    targets = formula_transform_targets(
      stats::setNames("affine", "mu_alloc__xXx__ablat[1]"),
      stats::setNames("identity", "mu_alloc__xXx__ablat[1]")
    )
  )
  class(transform) <- c("BayesTools_formula_coefficient_transform", "list")
  # Canonical level names carry the level label; the extraction key maps a
  # direct level cell to its coordinate. The reference level is structural and
  # a mean-difference level combines coordinates.
  level_names <- c(
    random     = "mu_alloc__xXx__ablat[random]",
    alternate  = "mu_alloc__xXx__ablat[alternate]",
    systematic = "mu_alloc__xXx__ablat[systematic]"
  )
  quantities <- data.frame(
    quantity_id      = paste0("q_", names(level_names)),
    canonical_name   = unname(level_names),
    component        = names(level_names),
    status           = c("sampled", "structural", "derived"),
    fixed_value      = c(NA_real_, 0, NA_real_),
    stringsAsFactors = FALSE
  )
  quantities[["extraction_key"]] <- list(
    list(type = "factor_level", dependencies = "mu_alloc__xXx__ablat[1]",
         weights = 1),
    list(type = "factor_level", dependencies = character(), weights = numeric()),
    list(type = "factor_level",
         dependencies = c("mu_alloc__xXx__ablat[1]", "mu_alloc__xXx__ablat[2]"),
         weights = c(0.5, 0.5))
  )
  selected_levels <- function(levels) {

    list(
      parameter = "mu_alloc__xXx__ablat",
      aliases   = list(
        "alloc:ablat"          = "mu_alloc__xXx__ablat",
        "mu_alloc__xXx__ablat" = "mu_alloc__xXx__ablat"
      ),
      component = "mods",
      entry = list(
        role                = "formula_coefficient_group",
        formula_parameter   = "mu",
        member_quantity_ids = list(paste0("q_", names(level_names)))
      ),
      resolution = list(occurrences = data.frame(
        level          = levels,
        canonical_name = unname(level_names[levels]),
        quantity_id    = paste0("q_", levels),
        stringsAsFactors = FALSE
      ))
    )
  }
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    parameter_catalog = function(...) list(quantities = quantities),
    .package = "BayesTools"
  )
  object <- list(fit = structure(list(), class = "BayesTools_fit"))

  target <- .hypothesis_brma_formula_coefficient_level_targets(
    object     = object,
    selected   = selected_levels(c("random", "alternate")),
    point_refs = data.frame(level = c("random", "alternate"), value = 0)
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
  # A level combining several coefficients has no single coefficient prior:
  # the stop names the level and the hypotheses that remain available.
  combined <- tryCatch(
    .hypothesis_brma_formula_coefficient_level_targets(
      object     = object,
      selected   = selected_levels("systematic"),
      point_refs = data.frame(level = "systematic", value = 0)
    ),
    error = conditionMessage
  )
  combined_message <- function(level) {

    paste0(
      "Point hypotheses on factor level 'alloc:ablat[", level, "]' are not ",
      "supported: the level is a linear combination of the fitted contrast ",
      "coefficients (mean-difference, orthonormal, or ordered contrasts), ",
      "not a fitted coefficient itself, and point hypotheses on a single ",
      "level require a level fitted as its own coefficient (treatment or ",
      "independent contrasts). Test the level with a region hypothesis such ",
      "as 'alloc:ablat[", level, "] > 0' or a level contrast such as ",
      "'alloc:ablat[", level, "] = alloc:ablat[random]'."
    )
  }
  expect_identical(combined, combined_message("systematic"))

  # A level whose key is one unit-weight coordinate that a contrast
  # coefficient '{2}' also holds (a unit design row of mean-difference coding
  # in floating point) is not a level cell.
  quantities <- rbind(
    quantities,
    data.frame(
      quantity_id      = c("q_unitrow", "q_coefficient2"),
      canonical_name   = c("mu_alloc__xXx__ablat[unitrow]", "mu_alloc__xXx__ablat{2}"),
      component        = c("unitrow", "{2}"),
      status           = c("derived", "sampled"),
      fixed_value      = NA_real_,
      extraction_key   = I(rep(list(list(
        type = "factor_level", dependencies = "mu_alloc__xXx__ablat[2]",
        weights = 1
      )), 2L)),
      stringsAsFactors = FALSE
    )
  )
  level_names <- c(level_names, unitrow = "mu_alloc__xXx__ablat[unitrow]")
  unit_row <- tryCatch(
    .hypothesis_brma_formula_coefficient_level_targets(
      object     = object,
      selected   = selected_levels("unitrow"),
      point_refs = data.frame(level = "unitrow", value = 0)
    ),
    error = conditionMessage
  )
  expect_identical(unit_row, combined_message("unitrow"))
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
      output_transforms = c(log_tau_intercept = "exp"),
      targets = formula_transform_targets(
        c(log_tau_intercept = "exp_affine"),
        c(log_tau_intercept = "exp")
      )
    )
  )
  class(target[["transform"]]) <- c(
    "BayesTools_formula_coefficient_transform",
    "list"
  )
  density <- structure(list(sentinel = TRUE), class = "prior_linear_density")
  testthat::local_mocked_bindings(
    JAGS_formula_prior_density = function(...) density,
    prior_ordinate_status = eligible_ordinate_status,
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
    output_transforms = c(log_tau_intercept = "exp"),
    targets = formula_transform_targets(
      c(log_tau_intercept = "exp_affine"),
      c(log_tau_intercept = "exp")
    )
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
  samples <- with_draw_metadata(
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
      structure(data.frame(
        Alternative = "log_tau_intercept = -1.6094379124341003",
        Null        = "log_tau_intercept != -1.6094379124341003"
      ), hypothesis_ast = hypothesis)
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
  expect_null(attr(result, "parsed", exact = TRUE))
  expect_equal(
    attr(result, "hypothesis_ast")$statements[[1L]][["left"]][["value"]],
    .2
  )
  expect_identical(
    attr(result, "hypothesis_ast")$statements[[1L]][["left"]][["label"]],
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
  BayesTools::posterior_metadata(
    atomic_samples[["log_tau_intercept"]],
    "atoms"
  ) <- BayesTools::posterior_atom_attribute(
    point_masses = data.frame(x = 0, mass = 0.25)
  )
  expect_error(
    .hypothesis_brma_exp_affine_certify(
      atomic_samples, "log_tau_intercept", FALSE
    ),
    "posterior is atom-free"
  )
  unproven_samples <- samples
  BayesTools::posterior_metadata(
    unproven_samples[["log_tau_intercept"]],
    "condition"
  ) <- NULL
  # Unconditioned draws are declared by 'averaged', not by the condition key.
  keyed_samples <- samples
  BayesTools::posterior_metadata(
    keyed_samples[["log_tau_intercept"]],
    "condition"
  ) <- list(
    conditional   = character(),
    condition_key = "<averaged>",
    averaged      = FALSE
  )
  expect_error(
    .hypothesis_brma_exp_affine_certify(
      keyed_samples, "log_tau_intercept", FALSE
    ),
    "structural evidence for an unconditional posterior"
  )
  expect_error(
    .hypothesis_brma_exp_affine_certify(
      unproven_samples, "log_tau_intercept", FALSE
    ),
    "structural evidence for an unconditional posterior"
  )
  atomic_prior_samples <- samples
  BayesTools::posterior_metadata(atomic_prior_samples, "prior_densities")[[
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
    prior_ordinate_status = function(prior_density, values, labels = NULL) {
      observed_values <<- c(observed_values, values)
      eligible_ordinate_status(prior_density, values)
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
    dependencies = data.frame(),
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
