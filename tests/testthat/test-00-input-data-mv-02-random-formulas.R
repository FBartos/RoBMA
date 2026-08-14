context("brma.mv random-formula input contracts")

test_that("brma.mv supports BayesTools random formulas and owns heterogeneity", {

  dat <- data.frame(
    yi        = c(0.10, 0.20, 0.30, 0.40),
    study     = c("s1", "s1", "s2", "s2"),
    criterion = c("c1", "c2", "c1", "c2")
  )

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    random                    = list(~ 1 | study, ~ 1 | criterion),
    data                      = dat,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  data           <- object[["data"]]
  random_effects <- attr(data[["location"]], "random_effects")
  terms          <- random_effects[["terms"]]

  expect_true(attr(data, "random"))
  expect_s3_class(random_effects, "BayesTools_random_effects")
  expect_null(attr(data[["location"]], "random_terms"))
  expect_match(
    paste(deparse(attr(data[["location"]], "formula")), collapse = " "),
    "Component_1"
  )
  expect_equal(
    vapply(terms, `[[`, character(1), "block_name"),
    c("Component_1", "Component_2")
  )
  expect_s3_class(object[["priors"]][["random"]], "prior_random")
  expect_equal(
    object[["priors"]][["random"]][["allocation"]][[1]][["terms"]],
    c(Component_1 = "Component_1", Component_2 = "Component_2")
  )
  expect_false("mu" %in% names(object[["priors"]][["outcome"]]))
  expect_false("tau" %in% names(object[["priors"]][["outcome"]]))
  expect_true("intercept" %in% names(object[["priors"]][["location"]]))

  syntax <- .create_model_syntax(data, object[["priors"]])
  expect_match(syntax, "yi\\[i\\] ~ dnorm\\(mu\\[i\\],1/\\( sampling_var\\[i\\] \\)\\)")
  expect_false(grepl("pow\\(tau,2\\)", syntax))

  design <- .fitted_formula_design(object, "mu", required = TRUE)
  expect_s3_class(design, "BayesTools_formula_design")
  expect_equal(
    vapply(design[["random_effects"]], `[[`, character(1), "block_name"),
    c("Component_1", "Component_2")
  )
})


test_that("brma.mv names nested random lists deterministically", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30, 0.40),
    district = c("d1", "d1", "d2", "d2"),
    school   = c("s1", "s2", "s1", "s2"),
    study    = c("s1", "s1", "s2", "s2")
  )

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    random                    = list(~ 1 | district / school),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  random_effects <- attr(object[["data"]][["location"]], "random_effects")
  terms          <- random_effects[["terms"]]
  block_names    <- vapply(terms, function(term) term[["block_name"]], character(1))

  expect_equal(
    block_names,
    c("Component_1_school_district", "Component_1_district")
  )
  expect_equal(
    object[["priors"]][["random"]][["allocation"]][[1]][["terms"]],
    c(
      school_district = "Component_1_school_district",
      district        = "Component_1_district"
    )
  )

  named <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    random                    = list(study_component = ~ 1 | study),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  named_terms <- attr(named[["data"]][["location"]], "random_effects")[["terms"]]

  expect_equal(
    vapply(named_terms, function(term) term[["block_name"]], character(1)),
    "study_component"
  )

  hierarchical <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    random                    = list(
      nested = ~ 1 | district / school,
      study  = ~ 1 | study
    ),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  hierarchical_terms      <- attr(
    hierarchical[["data"]][["location"]],
    "random_effects"
  )[["terms"]]
  hierarchical_allocation <- hierarchical[["priors"]][["random"]][["allocation"]]

  expect_equal(
    vapply(hierarchical_terms, function(term) term[["block_name"]], character(1)),
    c("nested_school_district", "nested_district", "study")
  )
  expect_length(hierarchical_allocation, 2L)
  expect_equal(
    hierarchical_allocation[[1]][["terms"]],
    c(nested = "nested", study = "study")
  )
  expect_equal(
    hierarchical_allocation[[2]][["parent"]][["allocation"]],
    "random_total"
  )
  expect_equal(
    hierarchical_allocation[[2]][["parent"]][["component"]],
    "nested"
  )
  expect_equal(
    hierarchical_allocation[[2]][["terms"]],
    c(school_district = "nested_school_district", district = "nested_district")
  )
})


test_that("brma.mv combines fixed moderators and random formulas in location", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30, 0.40),
    study = c("s1", "s1", "s2", "s2"),
    x     = c(0, 1, 2, 3)
  )

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    mods                      = ~ x,
    random                    = ~ 1 | study,
    data                      = dat,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  data <- object[["data"]]

  expect_true(attr(data, "mods"))
  expect_true(attr(data, "random"))
  expect_equal(names(object[["priors"]][["location"]]), c("intercept", "x"))
  expect_equal(names(object[["priors"]][["outcome"]]), character(0))

  design <- .fitted_formula_design(object, "mu", required = TRUE)
  expect_s3_class(design, "BayesTools_formula_design")
  expect_true("x" %in% design[["model_terms"]])
  expect_equal(
    vapply(design[["random_effects"]], `[[`, character(1), "block_name"),
    "study"
  )

  formula_text <- paste(deparse(attr(data[["location"]], "formula")), collapse = " ")
  expect_match(formula_text, "x", fixed = TRUE)
  expect_match(formula_text, "study", fixed = TRUE)

  syntax <- .create_model_syntax(data, object[["priors"]])
  expect_match(syntax, "yi\\[i\\] ~ dnorm\\(mu\\[i\\],1/\\( sampling_var\\[i\\] \\)\\)")
  expect_false(grepl("pow\\(tau,2\\)", syntax))

  identical_name <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    mods                      = ~ study,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  expect_true("study" %in% names(identical_name[["data"]][["location"]]))
  expect_null(.fitted_formula_design(identical_name, "mu", required = FALSE))

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      mods                      = data.frame(study = seq_len(4)),
      random                    = ~ 1 | study,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "different values"
  )
})


test_that("brma.mv supports BayesTools random covariance shapes", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30, 0.40),
    study = c("s1", "s1", "s2", "s2"),
    x     = c(0, 1, 0, 1),
    out   = c("a", "b", "a", "b"),
    time  = c(1, 2, 1, 2)
  )

  shapes <- list(
    id     = list(formula = ~ id(1 | study),             structure = "id"),
    diag   = list(formula = ~ diag(1 + x | study),       structure = "diag"),
    double = list(formula = ~ 1 + x || study,            structure = "diag"),
    us     = list(formula = ~ us(1 + x | study),         structure = "us"),
    un     = list(formula = ~ un(1 + x | study),         structure = "us"),
    cs     = list(formula = ~ cs(out | study),           structure = "cs"),
    hcs    = list(formula = ~ hcs(out | study),          structure = "hcs"),
    ar1    = list(formula = ~ ar1(time | study),         structure = "ar1"),
    ar     = list(formula = ~ ar(time | study),          structure = "ar1"),
    har    = list(formula = ~ har(time | study),         structure = "har"),
    car    = list(formula = ~ car(time | study),         structure = "car"),
    random = list(formula = ~ random(1 | study, name = "Named"), structure = "us"),
    re     = list(formula = ~ re(1 | study, name = "Named"),     structure = "us")
  )

  for (shape in names(shapes)) {
    object <- brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      random                    = shapes[[shape]][["formula"]],
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    )
    terms <- attr(object[["data"]][["location"]], "random_effects")[["terms"]]

    expect_equal(
      unique(vapply(terms, `[[`, character(1), "structure")),
      shapes[[shape]][["structure"]]
    )
    expect_s3_class(object[["priors"]][["random"]], "prior_random")
    expect_false("tau" %in% names(object[["priors"]][["outcome"]]))

    if (shapes[[shape]][["structure"]] %in% c("cs", "hcs", "ar1", "har", "car")) {
      block_name <- names(object[["priors"]][["random"]][["blocks"]])[[1L]]
      covariance <- object[["priors"]][["random"]][["blocks"]][[block_name]][["covariance"]]

      expect_equal(covariance[["rho_scale"]], "rho")
      expect_equal(covariance[["rho"]][["distribution"]], "beta")
      expect_equal(covariance[["rho"]][["parameters"]][["alpha"]], 1)
      expect_equal(covariance[["rho"]][["parameters"]][["beta"]], 1)
    }
  }
})


test_that("brma.mv preserves BayesTools random-effect parameterization policy", {

  dat <- data.frame(
    yi    = seq(-0.20, 0.25, length.out = 10L),
    study = rep(c("s1", "s2"), each = 5L)
  )
  sd_prior <- BayesTools::prior(
    "normal",
    parameters = list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )

  compile <- function(prior_heterogeneity) {

    object <- brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, nrow(dat))),
      random                    = ~ 1 | study,
      data                      = dat,
      measure                   = "GEN",
      prior_heterogeneity       = prior_heterogeneity,
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    )
    design <- .fitted_formula_design(object, "mu", required = TRUE)
    return(design[["random_effects"]][[1L]])
  }

  centered <- compile(BayesTools::prior_random(
    sd = sd_prior,
    parameterization = "centered"
  ))
  automatic <- compile(BayesTools::prior_random(
    sd = sd_prior,
    parameterization = "auto"
  ))
  overridden <- compile(BayesTools::prior_random(
    sd = sd_prior,
    parameterization = "centered",
    study = BayesTools::random_block(parameterization = "noncentered")
  ))

  expect_identical(centered[["parameterization_requested"]], "centered")
  expect_identical(centered[["parameterization_resolved"]], "centered")
  expect_identical(automatic[["parameterization_requested"]], "auto")
  expect_identical(automatic[["parameterization_resolved"]], "centered")
  expect_match(automatic[["parameterization_reason"]], "replication")
  expect_identical(overridden[["parameterization_requested"]], "noncentered")
  expect_identical(overridden[["parameterization_resolved"]], "noncentered")
})


test_that("brma.mv compiles sparse high-dimensional structured effects", {

  n_groups         <- 8L
  levels_per_group <- 5L
  n_outcomes       <- n_groups * levels_per_group
  dat <- data.frame(
    yi      = seq(-0.2, 0.2, length.out = n_outcomes),
    study   = factor(rep(seq_len(n_groups), each = levels_per_group)),
    outcome = factor(seq_len(n_outcomes))
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, nrow(dat))),
    random                    = ~ cs(outcome | study),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  compiled <- BayesTools::JAGS_formula(
    formula       = .create_fit_formula_list(
      data      = object[["data"]],
      parameter = "location"
    ),
    parameter     = "mu",
    data          = .create_fit_formula_data_list(
      data      = object[["data"]],
      parameter = "location"
    ),
    prior_list    = .create_fit_formula_prior_list(
      priors    = object[["priors"]],
      parameter = "location"
    ),
    formula_scale = .data_standardize_continuous_predictors(object[["data"]]),
    prior_random  = object[["priors"]][["random"]],
    random_effects_compile = .object_formula_random_effects_compile(
      object,
      "location"
    )
  )
  term   <- compiled[["formula_design"]][["random_effects"]][[1L]]
  syntax <- compiled[["formula_syntax"]]

  expect_identical(term[["parameterization_resolved"]], "noncentered")
  expect_identical(term[["latent_layout"]][["type"]], "group_local")
  expect_equal(term[["latent_layout"]][["n_local"]], nrow(dat))
  expect_match(syntax, "_xRE_MAPx", fixed = TRUE)
  expect_match(syntax, "_xRE_COLx", fixed = TRUE)
  expect_false(grepl("_xRE_DATAx", syntax, fixed = TRUE))
  expect_false(grepl("_xRE_CORx_L", syntax, fixed = TRUE))
  expect_false(grepl("_xRE_CORx_R", syntax, fixed = TRUE))
})


test_that("brma.mv validates random formula edge cases", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30, 0.40),
    study = c("s1", "s1", "s2", "s2"),
    x     = c(0, 1, 0, 1)
  )
  scale_prior <- BayesTools::prior(
    "normal",
    parameters = list(mean = 0, sd = 2)
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      scale                     = ~ missing_x,
      random                    = ~ 1 | study,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "missing_x"
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      random                    = ~ x | study,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "Plain 'random' terms"
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      scale                     = ~ x,
      random                    = ~(1 | study) + (1 | x),
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "Bare 'random' formulas with multiple random-effect terms"
  )

  expect_silent(
    named_block <- brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      random                    = list(
        study_component = ~ random(1 | study, name = "other_component")
      ),
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
  )
  named_random_effects <- attr(
    named_block[["data"]][["location"]],
    "random_effects"
  )
  expect_equal(
    named_random_effects[["components"]][["study_component"]],
    "other_component"
  )
  expect_equal(
    named_random_effects[["terms"]][[1]][["block_name"]],
    "other_component"
  )

  expect_silent(
    scale_random <- brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      scale                     = ~ x,
      random                    = ~ 1 | study,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    )
  )
  scale_random_design <- .fitted_formula_design(scale_random, "mu", required = TRUE)
  scale_random_binding <- scale_random_design[["random_effects"]][[1]][["sd_binding"]]
  expect_equal(scale_random_binding[["source"]][["name"]], "tau")
  expect_equal(scale_random_binding[["source"]][["shape"]], "row")
  expect_null(scale_random_binding[["source"]][["index"]])
  expect_null(scale_random_binding[["source"]][["expression"]])
  expect_null(scale_random_binding[["source"]][["row_indexed"]])
  expect_match(
    .create_model_syntax(scale_random[["data"]], scale_random[["priors"]]),
    "tau\\[i\\] = exp\\(log_tau\\[i\\]\\)"
  )
  expect_false(.summary_random_components_enabled(scale_random))

  component_scale_random <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    scale                     = list(
      study       = ~ x,
      x_component = ~ x
    ),
    random                    = list(
      study       = ~ 1 | study,
      x_component = ~ 1 | x
    ),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  component_scale_design <- .fitted_formula_design(
    component_scale_random,
    "mu",
    required = TRUE
  )
  component_scale_sources <- stats::setNames(
    vapply(component_scale_design[["random_effects"]], function(term) {
      term[["sd_binding"]][["source"]][["name"]]
    }, character(1)),
    vapply(component_scale_design[["random_effects"]], `[[`, character(1), "block_name")
  )
  expect_equal(component_scale_sources[["study"]], "tau_study")
  expect_equal(component_scale_sources[["x_component"]], "tau_x_component")
  expect_equal(
    .data_scale_formula_parameters(component_scale_random[["data"]]),
    c(study = "log_tau_study", x_component = "log_tau_x_component")
  )
  expect_s3_class(
    .fitted_formula_design(component_scale_random, "log_tau_study", required = TRUE),
    "BayesTools_formula_design"
  )
  component_scale_syntax <- .create_model_syntax(
    component_scale_random[["data"]],
    component_scale_random[["priors"]]
  )
  expect_match(component_scale_syntax, "tau_study\\[i\\] = exp\\(log_tau_study\\[i\\]\\)")
  expect_match(component_scale_syntax, "tau_x_component\\[i\\] = exp\\(log_tau_x_component\\[i\\]\\)")

  default_named_scale_random <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    scale                     = list(
      "Component 1" = ~ x,
      "Component 2" = ~ x
    ),
    random                    = list(~ 1 | study, ~ 1 | x),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  expect_equal(
    names(.data_scale_components(default_named_scale_random[["data"]])),
    c("Component_1", "Component_2")
  )
  expect_equal(
    .data_scale_formula_sources(default_named_scale_random[["data"]]),
    c(Component_1 = "tau_Component_1", Component_2 = "tau_Component_2")
  )

  default_named_scale_prior <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    scale                     = list(
      "Component 1" = ~ x,
      "Component 2" = ~ x
    ),
    random                    = list(~ 1 | study, ~ 1 | x),
    data                      = dat,
    measure                   = "GEN",
    prior_scale               = list(
      "Component 1" = list(x = scale_prior),
      "Component 2" = list(x = scale_prior)
    ),
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  expect_equal(
    names(default_named_scale_prior[["priors"]][["scale"]]),
    c("log_tau_Component_1", "log_tau_Component_2")
  )

  text_named_scale_prior <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    scale                     = list(
      "Study effects" = ~ x,
      "X effects"     = ~ x
    ),
    random                    = list(
      "Study effects" = ~ 1 | study,
      "X effects"     = ~ 1 | x
    ),
    data                      = dat,
    measure                   = "GEN",
    prior_scale               = list(
      "Study effects" = list(x = scale_prior),
      "X effects"     = list(x = scale_prior)
    ),
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  expect_equal(
    names(.data_scale_components(text_named_scale_prior[["data"]])),
    c("Study_effects", "X_effects")
  )
  expect_equal(
    vapply(
      .data_scale_component_specs(text_named_scale_prior[["data"]]),
      `[[`,
      character(1),
      "display_name"
    ),
    c(Study_effects = "Study effects", X_effects = "X effects")
  )

  plain_nested_scale_random <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    scale                     = ~ x,
    random                    = ~ 1 | study / x,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  plain_nested_effects <- attr(
    plain_nested_scale_random[["data"]][["location"]],
    "random_effects"
  )
  expect_equal(
    vapply(plain_nested_effects[["terms"]], `[[`, character(1), "component_label"),
    c("Component_1", "Component_1")
  )
  expect_equal(
    plain_nested_scale_random[["priors"]][["random"]][["allocation"]][[1]][["name"]],
    "total"
  )
  plain_nested_design <- .fitted_formula_design(
    plain_nested_scale_random,
    "mu",
    required = TRUE
  )
  expect_equal(
    vapply(plain_nested_design[["random_effects"]], function(term) {
      term[["sd_binding"]][["source"]][["name"]]
    }, character(1)),
    c("tau", "tau")
  )

  nested_scale_random <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    scale                     = ~ x,
    random                    = list(~ 1 | study / x),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  expect_true(.summary_random_components_enabled(nested_scale_random))
  expect_equal(
    nested_scale_random[["priors"]][["random"]][["allocation"]][[1]][["name"]],
    "total"
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      scale                     = ~ x,
      random                    = list(~ 1 | study, ~ 1 | x),
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "single 'scale' formula is ambiguous"
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      scale                     = list(study = ~ x, missing = ~ x),
      random                    = list(study = ~ 1 | study, x_component = ~ 1 | x),
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "Unknown: missing"
  )

  component_scale_prior <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    scale                     = list(
      study       = ~ x,
      x_component = ~ x
    ),
    random                    = list(
      study       = ~ 1 | study,
      x_component = ~ 1 | x
    ),
    data                      = dat,
    measure                   = "GEN",
    prior_scale               = list(
      study       = list(x = scale_prior),
      x_component = list(x = scale_prior)
    ),
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  expect_equal(
    names(component_scale_prior[["priors"]][["scale"]]),
    c("log_tau_study", "log_tau_x_component")
  )
  expect_true("x" %in% names(component_scale_prior[["priors"]][["scale"]][["log_tau_study"]]))

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      scale                     = list(study = ~ x, x_component = ~ x),
      random                    = list(study = ~ 1 | study, x_component = ~ 1 | x),
      data                      = dat,
      measure                   = "GEN",
      prior_scale               = list(x = scale_prior),
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "prior_scale' names must exactly match"
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      scale                     = list(study = ~ x, x_component = ~ x),
      random                    = list(study = ~ 1 | study, x_component = ~ 1 | x),
      data                      = dat,
      measure                   = "GEN",
      prior_scale               = list(
        study       = scale_prior,
        x_component = scale_prior
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "Each 'prior_scale' component"
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      data                      = dat,
      measure                   = "GEN",
      prior_heterogeneity       = BayesTools::prior_random(
        sd = BayesTools::prior("normal", list(mean = 0, sd = 1), list(0, Inf))
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "prior_random"
  )

  expect_silent(
    object <- brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      random                    = ~ 1 | study,
      data                      = dat,
      measure                   = "GEN",
      prior_heterogeneity       = BayesTools::prior_random(
        sd = BayesTools::prior("normal", list(mean = 0, sd = 1), list(0, Inf))
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    )
  )
  expect_s3_class(object[["priors"]][["random"]], "prior_random")

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      random                    = ~ 1 | study,
      data                      = dat,
      measure                   = "GEN",
      prior_heterogeneity       = BayesTools::prior_random(
        other = BayesTools::random_block(
          sd = BayesTools::prior("normal", list(mean = 0, sd = 1), list(0, Inf))
        )
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "other|study|random"
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      random                    = ~ 1 | study,
      data                      = dat,
      measure                   = "GEN",
      prior_heterogeneity       = BayesTools::prior_random(
        study = BayesTools::random_block(
          sd_source = BayesTools::random_sd_source("tau", shape = "row")
        )
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "Row-shaped SD source 'tau'"
  )

  expect_silent(
    custom_scale_random <- brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      scale                     = ~ x,
      random                    = ~ 1 | study,
      data                      = dat,
      measure                   = "GEN",
      prior_heterogeneity       = BayesTools::prior_random(
        study = BayesTools::random_block(
          sd_source = BayesTools::random_sd_source("tau", shape = "row")
        )
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    )
  )
  custom_scale_binding <- .fitted_formula_design(
    custom_scale_random,
    "mu",
    required = TRUE
  )[["random_effects"]][[1]][["sd_binding"]]
  expect_equal(custom_scale_binding[["source"]][["name"]], "tau")
  expect_equal(custom_scale_binding[["source"]][["shape"]], "row")

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      scale                     = ~ x,
      random                    = ~ 1 | study,
      data                      = dat,
      measure                   = "GEN",
      prior_heterogeneity       = BayesTools::prior_random(
        sd = BayesTools::prior("normal", list(mean = 0, sd = 1), list(0, Inf))
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "row-shaped SD source 'tau'"
  )
})


test_that("brma.mv marginalizes one-to-one random intercepts into likelihood variance", {

  dat <- data.frame(
    yi    = c(0.10, 0.20),
    study = c("s1", "s2")
  )
  V <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2)

  latent <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  latent_design <- .fitted_formula_design(latent, "mu", required = TRUE)
  expect_equal(
    vapply(latent_design[["random_effects"]], `[[`, character(1), "compile_mode"),
    "marginalized"
  )
  latent_syntax <- .create_model_syntax(latent[["data"]], latent[["priors"]])
  expect_match(latent_syntax, "sampling_dependency", fixed = TRUE)
  expect_match(
    latent_syntax,
    "yi\\[i\\] ~ dnorm\\(mu\\[i\\] \\+ sampling_dependency\\[i\\],1/\\( sampling_var\\[i\\] \\+ pow\\(mu__xREx__study_intercept,2\\) \\)\\)"
  )
  latent_posterior <- matrix(
    c(0.50, 1.00),
    ncol     = 1,
    dimnames = list(NULL, "mu__xREx__study_intercept")
  )
  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = latent[["data"]],
      posterior_samples = latent_posterior
    ),
    matrix(c(0.25, 1.00, 0.25, 1.00), nrow = 2)
  )

  formula_args <- .create_jags_formula_args(
    data   = latent[["data"]],
    priors = latent[["priors"]]
  )
  expect_named(formula_args[["formula_list"]], "mu")
  expect_named(formula_args[["formula_random_prior_list"]], "mu")
  expect_named(formula_args[["formula_random_effects_compile_list"]], "mu")

  whitened <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat,
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  whitened_syntax <- .create_model_syntax(whitened[["data"]], whitened[["priors"]])
  expect_match(
    whitened_syntax,
    "whitening_y_1\\[j\\] ~ dnorm\\(whitening_mu_1\\[j\\], 1/\\( whitening_var_1\\[j\\] \\+ pow\\(mu__xREx__study_intercept,2\\) \\)\\)"
  )

  block_mvn <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  block_mvn_syntax <- .create_model_syntax(block_mvn[["data"]], block_mvn[["priors"]])
  expect_match(
    block_mvn_syntax,
    "tau2_observed\\[i\\] = pow\\(mu__xREx__study_intercept,2\\)"
  )
  expect_match(block_mvn_syntax, "dknown_v_mnorm", fixed = TRUE)
  expect_false(grepl("pow\\(tau,2\\)", block_mvn_syntax))
})


test_that("marginalized row SD sources require newdata-shaped source samples", {

  term <- list(
    block_name         = "study",
    sd_parameter_names = NA_character_,
    sd_binding         = list(
      source = list(name = "tau", shape = "row")
    )
  )
  data <- list(
    outcome  = data.frame(yi = seq_len(4)),
    location = data.frame(row = seq_len(4))
  )
  attr(data[["location"]], "marginalized_random_effects") <- list(term)
  posterior_samples <- matrix(
    c(.10, .20, .30, .40),
    nrow     = 1L,
    dimnames = list(NULL, paste0("tau[", 1:4, "]"))
  )

  expect_error(
    .evaluate_marginalized_random_variance(
      data              = data,
      posterior_samples = posterior_samples,
      K                 = 2L
    ),
    "different row count"
  )
  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = data,
      posterior_samples = posterior_samples,
      K                 = 2L,
      source_samples    = list(tau = matrix(c(.50, .60), nrow = 1L))
    ),
    matrix(c(.25, .36), nrow = 1L)
  )
})


test_that("brma.mv can keep one-to-one random intercepts sampled", {

  dat <- data.frame(
    yi    = c(0.10, 0.20),
    study = c("s1", "s2")
  )
  V <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2)

  object <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    marginalize_estimate_level = FALSE,
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  design <- .fitted_formula_design(object, "mu", required = TRUE)
  expect_equal(
    vapply(design[["random_effects"]], `[[`, character(1), "compile_mode"),
    "sampled"
  )

  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])
  expect_match(
    syntax,
    "yi\\[i\\] ~ dnorm\\(mu\\[i\\] \\+ sampling_dependency\\[i\\],1/\\( sampling_var\\[i\\] \\)\\)"
  )
  expect_false(grepl("mu__xREx__study_intercept,2", syntax, fixed = TRUE))
  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = object[["data"]],
      posterior_samples = matrix(numeric(0), nrow = 2, ncol = 0)
    ),
    matrix(0, nrow = 2, ncol = 2)
  )
})


test_that("brma.mv marginalizes nested estimate level with allocation child SD", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30, 0.40),
    study    = c("s1", "s1", "s2", "s2"),
    estimate = paste0("e", 1:4)
  )

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    random                    = ~ 1 | study / estimate,
    data                      = dat,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  design <- .fitted_formula_design(object, "mu", required = TRUE)
  expect_equal(
    vapply(design[["random_effects"]], `[[`, character(1), "block_name"),
    c("estimate_study", "study")
  )
  expect_equal(
    vapply(design[["random_effects"]], `[[`, character(1), "compile_mode"),
    c("marginalized", "sampled")
  )

  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])
  expect_match(
    syntax,
    "sampling_var\\[i\\] \\+ pow\\(mu__xREx__estimate_study_intercept,2\\)"
  )
  expect_false(grepl(
    "pow(mu__xRE_ALLOCx_random_total__total_sd,2)",
    syntax,
    fixed = TRUE
  ))

  posterior <- matrix(
    c(0.20, 0.30),
    ncol     = 1,
    dimnames = list(NULL, "mu__xREx__estimate_study_intercept")
  )
  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = object[["data"]],
      posterior_samples = posterior
    ),
    matrix(rep(c(0.04, 0.09), times = 4), nrow = 2)
  )
  target <- .estimate_log_lik_target_metadata(
    setup     = list(data = object[["data"]], K = nrow(dat)),
    data_hash = .get_outcome_hash(object)
  )
  expect_equal(.data_sampled_random_effect_blocks(object[["data"]]), "study")
  expect_false(.data_has_sampled_estimate_level_random_effects(object[["data"]]))
  expect_equal(target[["random_effects"]], "conditioned")
  expect_equal(target[["estimate_level_random"]], "marginalized")
})


test_that("brma.mv marginalized estimate-level scale uses row SD source", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30, 0.40),
    estimate = paste0("e", 1:4),
    x        = c(0, 1, 0, 1)
  )

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    scale                     = ~ x,
    random                    = ~ 1 | estimate,
    data                      = dat,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])
  expect_match(syntax, "tau\\[i\\] = exp\\(log_tau\\[i\\]\\)")
  expect_match(syntax, "sampling_var\\[i\\] \\+ pow\\(tau\\[i\\],2\\)")

  posterior <- matrix(
    c(
      0.10, 0.20, 0.30, 0.40,
      0.50, 0.60, 0.70, 0.80
    ),
    nrow     = 2,
    byrow    = TRUE,
    dimnames = list(NULL, paste0("tau[", 1:4, "]"))
  )
  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = object[["data"]],
      posterior_samples = posterior
    ),
    posterior^2
  )

  formula_args <- .create_jags_formula_args(
    data   = object[["data"]],
    priors = object[["priors"]]
  )
  expect_true("tau" %in% formula_args[["add_parameters"]])
})


test_that("brma.mv supports partial random scale with sampled random component", {

  dat <- data.frame(
    yi   = c(0.10, 0.20, 0.30, 0.40, 0.15, 0.25),
    vi   = c(0.04, 0.05, 0.06, 0.07, 0.02, 0.03),
    type = factor(
      c("RCT", "RCT", "RCT", "cohort", "cohort", "cohort"),
      levels = c("RCT", "cohort")
    )
  )
  dat$study   <- factor(seq_len(nrow(dat)))
  dat$cohort  <- as.numeric(dat$type == "cohort")
  dat$bias_id <- factor("cohort_bias")

  object <- brma.mv(
    yi                        = yi,
    V                         = vi,
    random                    = list(
      coh_bias    = ~ diag(0 + cohort | bias_id),
      ran_effects = ~ 1 | study
    ),
    scale                     = list(ran_effects = ~ type),
    data                      = dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_equal(
    .data_scale_formula_sources(object[["data"]]),
    c(ran_effects = "tau_ran_effects")
  )
  expect_equal(
    .data_scale_formula_parameters(object[["data"]]),
    c(ran_effects = "log_tau_ran_effects")
  )
  expect_equal(
    .data_sampled_random_effect_blocks(object[["data"]]),
    "coh_bias"
  )

  marginalized_effects <- .data_marginalized_random_effects(object[["data"]])
  expect_length(marginalized_effects, 1L)
  expect_equal(marginalized_effects[[1L]][["block_name"]], "ran_effects")

  design <- .fitted_formula_design(object, "mu", required = TRUE)
  random_terms <- stats::setNames(
    design[["random_effects"]],
    vapply(design[["random_effects"]], `[[`, character(1), "block_name")
  )
  expect_equal(random_terms[["coh_bias"]][["compile_mode"]], "sampled")
  expect_null(random_terms[["coh_bias"]][["sd_binding"]][["source"]])
  expect_equal(random_terms[["ran_effects"]][["compile_mode"]], "marginalized")
  ran_effects_source <- random_terms[["ran_effects"]][["sd_binding"]][["source"]]
  expect_equal(ran_effects_source[["name"]], "tau_ran_effects")
  expect_equal(ran_effects_source[["shape"]], "row")

  random_priors <- object[["priors"]][["random"]][["blocks"]]
  expect_s3_class(random_priors[["coh_bias"]][["sd"]], "prior")
  expect_null(random_priors[["coh_bias"]][["sd_source"]])
  expect_null(random_priors[["ran_effects"]][["sd"]])
  expect_equal(
    random_priors[["ran_effects"]][["sd_source"]][["name"]],
    "tau_ran_effects"
  )

  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])
  expect_match(
    syntax,
    "tau_ran_effects\\[i\\] = exp\\(log_tau_ran_effects\\[i\\]\\)"
  )
  expect_match(
    syntax,
    "sampling_var\\[i\\] \\+ pow\\(tau_ran_effects\\[i\\],2\\)"
  )

  posterior <- matrix(
    c(
      0.10, 0.20, 0.30, 0.40, 0.50, 0.60,
      0.15, 0.25, 0.35, 0.45, 0.55, 0.65
    ),
    nrow     = 2,
    byrow    = TRUE,
    dimnames = list(NULL, paste0("tau_ran_effects[", 1:6, "]"))
  )
  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = object[["data"]],
      posterior_samples = posterior
    ),
    posterior^2
  )

  formula_args <- .create_jags_formula_args(
    data   = object[["data"]],
    priors = object[["priors"]]
  )
  expect_true("tau_ran_effects" %in% formula_args[["add_parameters"]])
})


test_that("brma.mv marginalized random scale applies allocation weights", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30, 0.40),
    study    = c("s1", "s1", "s2", "s2"),
    estimate = paste0("e", 1:4),
    x        = c(0, 1, 0, 1)
  )

  object <- brma.mv(
    yi                         = yi,
    V                          = diag(rep(0.04, 4)),
    scale                      = ~ x,
    random                     = ~ 1 | study / estimate,
    data                       = dat,
    known_v_parameterization   = "block_mvn",
    measure                    = "GEN",
    prior_unit_information_sd  = 1,
    only_priors                = TRUE
  )

  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])
  expect_match(
    syntax,
    "sampling_var\\[i\\] \\+ pow\\(tau\\[i\\] \\* sqrt\\(mu__xRE_ALLOCx_total__weight\\[1\\]\\),2\\)"
  )

  posterior <- matrix(
    c(
      0.10, 0.20, 0.30, 0.40, 0.25, 0.75,
      0.50, 0.60, 0.70, 0.80, 0.60, 0.40
    ),
    nrow     = 2,
    byrow    = TRUE,
    dimnames = list(NULL, c(
      paste0("tau[", 1:4, "]"),
      "mu__xRE_ALLOCx_total__weight[1]",
      "mu__xRE_ALLOCx_total__weight[2]"
    ))
  )
  expected <- posterior[, paste0("tau[", 1:4, "]"), drop = FALSE]^2 *
    posterior[, "mu__xRE_ALLOCx_total__weight[1]"]

  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = object[["data"]],
      posterior_samples = posterior
    ),
    expected
  )
})


test_that("brma.mv rejects ambiguous one-to-one random intercept marginalization", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30, 0.40),
    study    = paste0("s", 1:4),
    estimate = paste0("e", 1:4)
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(rep(0.04, 4)),
      random                    = list(study = ~ 1 | study, estimate = ~ 1 | estimate),
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "Multiple random-effect blocks map one-to-one"
  )
})


test_that("brma.mv random residual methods expose supported known-V targets", {

  dat <- data.frame(
    yi    = c(0.10, 0.20),
    study = c("s1", "s2")
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 2)),
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_silent(
    .check_residual_random_formula_availability(
      object             = object,
      type               = "outcome",
      conditioning_depth = "marginal",
      caller             = "test()"
    )
  )
  expect_silent(
    .check_residual_random_formula_availability(
      object             = object,
      type               = "pearson",
      conditioning_depth = "marginal",
      caller             = "test()"
    )
  )
  expect_silent(
    .check_residual_random_formula_availability(
      object             = object,
      type               = "rstandard",
      conditioning_depth = "marginal",
      caller             = "test()"
    )
  )
  expect_silent(
    .check_residual_random_formula_availability(
      object             = object,
      type               = "pearson",
      conditioning_depth = "estimate",
      caller             = "test()"
    )
  )
  expect_silent(
    .check_residual_random_formula_availability(
      object             = object,
      type               = "rstandard",
      conditioning_depth = "estimate",
      caller             = "test()"
    )
  )
  expect_silent(
    .check_residual_random_formula_availability(
      object             = object,
      type               = "rstudent",
      conditioning_depth = "estimate",
      caller             = "test()"
    )
  )
  expect_error(
    hatvalues(object),
    "posterior|prior|fit|samples|LOO"
  )
  expect_error(
    dfbetas(object),
    "LOO has not been computed|posterior|prior|fit|samples"
  )
  expect_error(
    covratio(object),
    "LOO has not been computed|posterior|prior|fit|samples"
  )
  expect_error(
    dffits(object),
    "LOO has not been computed|posterior|prior|fit|samples"
  )
  expect_error(
    cooks.distance(object),
    "LOO has not been computed|posterior|prior|fit|samples"
  )
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "test()")
  )
  expect_error(
    loo(object),
    "LOO has not been computed"
  )
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "add_waic()")
  )
  expect_silent(
    .check_marglik_available(object, "test()")
  )

  allocation_object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    scale                     = ~ x,
    random                    = ~ 1 | study / estimate,
    data                      = data.frame(
      yi       = c(0.10, 0.20, 0.30, 0.40),
      study    = c("s1", "s1", "s2", "s2"),
      estimate = paste0("e", 1:4),
      x        = c(0, 1, 0, 1)
    ),
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  expect_silent(
    .check_marglik_available(allocation_object, "test()")
  )
})


test_that("brma.mv diagnostic target registry documents implemented semantics", {

  targets <- .brma_mv_diagnostic_target_table()

  expect_s3_class(targets, "data.frame")
  expect_named(
    targets,
    c("method", "status", "target", "known_v_semantics", "known_r_semantics")
  )
  expect_equal(anyDuplicated(targets[["method"]]), 0L)
  expect_true(all(c("implemented", "deferred") %in% targets[["status"]]))

  loo_target <- .brma_mv_target_row("add_loo()/loo()")
  expect_equal(loo_target[["status"]], "implemented")
  expect_equal(loo_target[["target"]], "estimate-unit log-score")
  expect_true(grepl("p(y_i | y_-i", loo_target[["known_v_semantics"]],
                    fixed = TRUE))
  expect_true(grepl("conditioned", loo_target[["known_r_semantics"]],
                    fixed = TRUE))

  marglik_target <- .brma_mv_target_row("add_marglik()/bridge_sampler()")
  expect_equal(marglik_target[["status"]], "implemented")
  expect_equal(marglik_target[["target"]], "full joint fitted likelihood")
  expect_true(grepl("latent random-effect prior",
                    marglik_target[["known_r_semantics"]], fixed = TRUE))

  known_v_object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  known_v_backend <- .brma_mv_known_v_backend_metadata(known_v_object)
  expect_true(known_v_backend[["known_v"]])
  expect_equal(known_v_backend[["known_v_parameterization_requested"]], "auto")
  expect_equal(known_v_backend[["known_v_parameterization"]], "whitened")
  expect_equal(known_v_backend[["known_v_effective_backend"]], "whitened")
  marglik_metadata <- .brma_mv_marglik_target_metadata(known_v_object)
  expect_equal(marglik_metadata[["reported_target"]], "full joint fitted likelihood")
  expect_equal(marglik_metadata[["known_v_parameterization_requested"]], "auto")
  expect_equal(marglik_metadata[["known_v_parameterization"]], "whitened")

  diagonal_object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = c(0.04, 0.09),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  diagonal_backend <- .brma_mv_known_v_backend_metadata(diagonal_object)
  expect_equal(diagonal_backend[["known_v_parameterization"]], "whitened")
  expect_equal(diagonal_backend[["known_v_effective_backend"]], "diagonal")

  dffits_target <- .brma_mv_target_row("dffits()")
  expect_equal(dffits_target[["status"]], "implemented")
  expect_equal(dffits_target[["target"]], "fixed-location fitted-value influence")
  expect_true(grepl("fixed_location_fitted_value", dffits_target[["known_v_semantics"]],
                    fixed = TRUE))
})


test_that("random-effect structure preserves every represented design change", {

  allocation <- list(leaf_index_by_column = 1L)
  term <- list(model_matrix = matrix(c(0, 1e-15, 0), ncol = 1L))
  expect_identical(
    .brma_mv_sd_component_allocation_rows(
      allocation   = allocation,
      term         = term,
      component    = 1L,
      source_shape = "row",
      K            = 3L
    ),
    2L
  )

  intercept_term <- list(
    n_columns    = 1L,
    structure    = "us",
    model_matrix = matrix(1, nrow = 3L),
    group_map    = 1:3,
    n_groups     = 3L
  )
  expect_true(.is_estimate_level_random_intercept(intercept_term, k = 3L))
  intercept_term[["model_matrix"]][2L] <- 1 + .Machine$double.eps
  expect_false(.is_estimate_level_random_intercept(intercept_term, k = 3L))

  object <- list(formula_design = list(mu = list(
    random_effects = list(intercept_term)
  )))
  expect_false(.regplot_mv_random_intercept_only(object))

  multiplier_term <- list(
    row_multiplier      = c(1, 1 + .Machine$double.eps),
    row_multiplier_name = "known_r_multiplier"
  )
  expect_true(
    .marginalized_random_effect_row_multiplier_varying(multiplier_term)
  )
})
