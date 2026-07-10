context("Input handling for brma.mv")

test_that("brma.mv stores and decomposes known V", {

  V <- matrix(
    c(
      0.04, 0.01, 0.00,
      0.01, 0.09, 0.00,
      0.00, 0.00, 0.16
    ),
    nrow = 3,
    byrow = TRUE
  )

  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = V,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  data <- object[["data"]]

  expect_true(attr(data, "known_V"))
  expect_equal(data[["outcome"]][["sei"]], sqrt(diag(V)))

  known_V <- attr(data, "known_V_data")
  expect_equal(known_V[["V"]], V)
  expect_equal(
    diag(known_V[["residual_variance"]], nrow = 3) + tcrossprod(known_V[["B"]]),
    V,
    tolerance = 1e-10
  )
  expect_equal(known_V[["rank"]], ncol(known_V[["B"]]))
})


test_that("brma.mv supports list V input", {

  V1 <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2)
  V2 <- matrix(0.16, nrow = 1)

  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = list(V1, V2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  expect_equal(
    attr(object[["data"]], "known_V_data")[["V"]],
    rbind(cbind(V1, matrix(0, 2, 1)), cbind(matrix(0, 1, 2), V2))
  )
})


test_that("known V list input rejects empty blocks", {

  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = list(matrix(numeric(0), nrow = 0L, ncol = 0L)),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "non-empty square"
  )
})


test_that("known V symmetry tolerance scales with covariance magnitude", {

  V <- matrix(
    c(1e8, 1e4, 1e4 + 1e-4, 2e8),
    nrow = 2L
  )
  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = V,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  expect_equal(
    attr(object[["data"]], "known_V_data")[["V"]],
    (V + t(V)) / 2
  )
})


test_that("V_new list input rejects empty blocks", {

  expect_error(
    .known_v_newdata_prepare(
      V_new = list(matrix(numeric(0), nrow = 0L, ncol = 0L)),
      k     = 0L
    ),
    "V_new.*non-empty square"
  )
})


test_that("V_new symmetry tolerance scales with covariance magnitude", {

  V_new <- matrix(
    c(1e8, 1e4, 1e4 + 1e-4, 2e8),
    nrow = 2L
  )
  out <- .known_v_newdata_prepare(V_new, k = 2L)

  expect_equal(out[["V"]], (V_new + t(V_new)) / 2)
})


test_that("brma.mv supports variance vector V input", {

  vi <- c(0.04, 0.09, 0.16)
  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = vi,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  expect_equal(
    attr(object[["data"]], "known_V_data")[["V"]],
    diag(vi)
  )
  expect_equal(object[["data"]][["outcome"]][["sei"]], sqrt(vi))
})


test_that("brma.mv supports hidden vi and sei diagonal known-V input", {

  vi <- c(0.04, 0.09, 0.16)

  vi_object <- brma.mv(
    yi                        = c(0.10, 0.20, -0.05),
    vi                        = vi,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  sei_object <- brma.mv(
    yi                        = c(0.10, 0.20, -0.05),
    sei                       = sqrt(vi),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  both_object <- brma.mv(
    yi                        = c(0.10, 0.20, -0.05),
    vi                        = vi,
    sei                       = sqrt(vi),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  expect_equal(.data_known_v_data(vi_object[["data"]])[["V"]], diag(vi))
  expect_equal(.data_known_v_data(sei_object[["data"]])[["V"]], diag(vi))
  expect_equal(.data_known_v_data(both_object[["data"]])[["V"]], diag(vi))
  expect_error(
    brma.mv(
      yi        = c(0.10, 0.20, -0.05),
      V         = vi,
      vi        = vi,
      measure   = "GEN",
      only_data = TRUE
    ),
    "Use only one of 'V' and hidden 'vi'/'sei'"
  )
  expect_error(
    brma.mv(
      yi        = c(0.10, 0.20, -0.05),
      vi        = vi,
      sei       = sqrt(vi) + 0.01,
      measure   = "GEN",
      only_data = TRUE
    ),
    "must be consistent"
  )
})


test_that("internal mv data input defaults NULL known-V parameterization to auto", {

  yi <- c(0.10, 0.20)
  V  <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2)

  data <- .check_and_list_data(
    .call                              = quote(brma.mv(yi = yi, V = V)),
    .envir                             = environment(),
    class                              = "mv",
    measure                            = "GEN",
    set_contrast_factor_predictors    = "treatment",
    standardize_continuous_predictors = FALSE,
    random_group_covariance            = NULL,
    known_v_parameterization            = NULL,
    known_v_residual_fraction           = NULL,
    known_v_residual_fraction_specified = FALSE
  )
  known_V <- .data_known_v_data(data)

  expect_equal(known_V[["parameterization_requested"]], "auto")
  expect_equal(known_V[["parameterization"]], "whitened")
})


test_that("brma.mv warns and accepts rank-one all-correlated known V", {

  sei <- c(0.20, 0.30, 0.40)
  V   <- tcrossprod(sei)

  expect_warning(
    object <- brma.mv(
      yi                        = c(0.10, 0.20, 0.30),
      V                         = V,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "positive semidefinite"
  )

  known_V <- attr(object[["data"]], "known_V_data")

  expect_equal(known_V[["V"]], V)
  expect_equal(known_V[["parameterization"]], "whitened")
  expect_equal(crossprod(known_V[["sampling_factor"]]), V, tolerance = 1e-10)
  expect_lte(min(abs(known_V[["diagnostics"]][["min_whitening_variance"]])), 1e-10)

  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20, 0.30),
      V                         = V,
      known_v_parameterization  = "latent",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "cannot use known_v_parameterization = 'latent'"
  )

  dat <- data.frame(
    yi = c(0.10, 0.20, 0.30),
    x  = c(0, 1, 2)
  )
  old_max_block <- getOption("RoBMA.known_v_block_mvn_max_block_size", NULL)
  on.exit({
    options(RoBMA.known_v_block_mvn_max_block_size = old_max_block)
  }, add = TRUE)
  options(RoBMA.known_v_block_mvn_max_block_size = 2L)

  expect_warning(
    scale_object <- brma.mv(
      yi                        = yi,
      V                         = V,
      scale                     = ~ x,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "positive semidefinite"
  )
  expect_equal(
    attr(scale_object[["data"]], "known_V_data")[["parameterization"]],
    "block_mvn"
  )
})


test_that("brma.mv aligns V after subset and missing rows", {

  V <- matrix(
    c(
      0.04, 0.01, 0.00, 0.00,
      0.01, 0.09, 0.02, 0.00,
      0.00, 0.02, 0.16, 0.03,
      0.00, 0.00, 0.03, 0.25
    ),
    nrow = 4,
    byrow = TRUE
  )
  dat <- data.frame(
    yi    = c(0.10, NA, 0.30, 0.40),
    x     = c(0, 1, 1, 0),
    study = c("a", "b", "c", "d"),
    label = paste0("Study ", 1:4)
  )

  expect_warning(
    object <- brma.mv(
      yi                        = yi,
      V                         = V,
      mods                      = ~ x,
      random                    = ~ 1 | study,
      data                      = dat,
      slab                      = label,
      subset                    = c(TRUE, TRUE, TRUE, TRUE),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "removed due to missing values"
  )

  expected_V <- V[c(1, 3, 4), c(1, 3, 4), drop = FALSE]
  expect_equal(attr(object[["data"]], "known_V_data")[["V"]], expected_V)
  expect_equal(object[["data"]][["outcome"]][["yi"]], c(0.10, 0.30, 0.40))
  expect_equal(object[["data"]][["outcome"]][["slab"]],
               c("Study 1", "Study 3", "Study 4"))
  expect_true(attr(object[["data"]], "slab"))
  expect_equal(as.character(object[["data"]][["location"]][["study"]]), c("a", "c", "d"))

  expect_silent(
    subset_object <- brma.mv(
      yi                        = yi,
      V                         = V,
      mods                      = ~ x,
      random                    = ~ 1 | study,
      data                      = dat,
      slab                      = label,
      subset                    = c(TRUE, FALSE, TRUE, FALSE),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
  )
  expect_equal(
    attr(subset_object[["data"]], "known_V_data")[["V"]],
    V[c(1, 3), c(1, 3), drop = FALSE]
  )
  expect_equal(subset_object[["data"]][["outcome"]][["yi"]], c(0.10, 0.30))
  expect_equal(subset_object[["data"]][["outcome"]][["slab"]],
               c("Study 1", "Study 3"))
  expect_equal(as.character(subset_object[["data"]][["location"]][["study"]]),
               c("a", "c"))
})


test_that("brma.mv rejects invalid and unsupported inputs", {

  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = matrix(c(1, 2, 2, 1), nrow = 2),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "positive definite"
  )

  sei <- c(0.20, 0.30, 0.40)
  R   <- matrix(
    c(
      1,  0,  0,
      0,  1, -1,
      0, -1,  1
    ),
    nrow  = 3,
    byrow = TRUE
  )
  V_singular_not_ones <- diag(sei) %*% R %*% diag(sei)
  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20, 0.30),
      V                         = V_singular_not_ones,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "positive definite"
  )

  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = diag(2),
      weights                   = c(1, 1),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "'weights' are not supported"
  )

  dat_unsupported <- data.frame(
    yi    = c(0.10, 0.20),
    study = c("s1", "s2")
  )
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(2),
      cluster                   = study,
      data                      = dat_unsupported,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "'cluster' is not supported"
  )

  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = diag(2),
      prior_bias                = TRUE,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "Selection/publication-bias priors are not supported"
  )

  expect_error(
    brma.mv(
      yi                         = c(0.10, 0.20),
      V                          = diag(2),
      known_v_parameterization   = "block_mvn",
      known_v_residual_fraction  = NA_real_,
      measure                    = "GEN",
      prior_unit_information_sd  = 1,
      only_data                  = TRUE
    ),
    "known_v_residual_fraction"
  )
})


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
    "whitening_y\\[j\\] ~ dnorm\\(whitening_mu\\[j\\], 1/\\( whitening_var\\[j\\] \\+ pow\\(mu__xREx__study_intercept,2\\) \\)\\)"
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
    "tau2_observed\\[i\\] = pow\\(tau_ran_effects\\[i\\],2\\)"
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
    "tau2_observed\\[i\\] = pow\\(tau\\[i\\] \\* sqrt\\(mu__xRE_ALLOCx_total__weight\\[1\\]\\),2\\)"
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
  expect_match(
    .brma_mv_known_v_backend_label(known_v_backend),
    "requested: auto",
    fixed = TRUE
  )
  marglik_metadata <- .brma_mv_marglik_target_metadata(known_v_object)
  expect_equal(marglik_metadata[["reported_target"]], "full joint fitted likelihood")
  expect_equal(marglik_metadata[["known_v_parameterization_requested"]], "auto")
  expect_equal(marglik_metadata[["known_v_parameterization"]], "whitened")

  dffits_target <- .brma_mv_target_row("dffits()")
  expect_equal(dffits_target[["status"]], "implemented")
  expect_equal(dffits_target[["target"]], "fixed-location fitted-value influence")
  expect_true(grepl("fixed_location_fitted_value", dffits_target[["known_v_semantics"]],
                    fixed = TRUE))
})


test_that("brma.mv honors forced known-V backend requests", {

  V <- matrix(c(0.04, 0.015, 0.015, 0.09), nrow = 2)

  for (parameterization in c("latent", "whitened", "block_mvn")) {
    object <- brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = V,
      known_v_parameterization  = parameterization,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    )
    known_V <- attr(object[["data"]], "known_V_data")

    expect_equal(known_V[["parameterization_requested"]], parameterization)
    expect_equal(known_V[["parameterization"]], parameterization)
    expect_equal(.data_known_v_parameterization(object[["data"]]), parameterization)
  }
})


test_that("LOO target comparison rejects missing likelihood target labels", {

  metadata <- list(
    unit               = "estimate",
    conditioning_depth = "estimate",
    data_hash          = "same-data",
    target             = "known_v_estimate"
  )
  loo_a <- structure(list(), class = "loo")
  loo_b <- structure(list(), class = "loo")
  attr(loo_a, "RoBMA_target") <- metadata
  attr(loo_b, "RoBMA_target") <- metadata[names(metadata) != "target"]

  expect_error(
    .check_loo_compare_targets(list(loo_a, loo_b)),
    "without likelihood target labels"
  )
})


test_that("brma.mv adds latent known-V factors to fit syntax", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2),
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  fit_priors <- .create_fit_priors(object[["data"]], object[["priors"]])
  fit_data   <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax     <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_true("sampling_z" %in% names(fit_priors))
  expect_true("sampling_var" %in% names(fit_data))
  expect_true("sampling_B" %in% names(fit_data))
  expect_match(syntax, "sampling_dependency", fixed = TRUE)
  expect_match(syntax, "sampling_var\\[i\\]")
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "test()")
  )
  expect_error(
    .check_log_lik_target_available(object, "cluster", "test()"),
    "unit = 'cluster'"
  )
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "add_waic()")
  )
})


test_that("diagonal brma.mv estimate likelihood matches brma.norm factorization", {

  yi <- c(0.10, 0.20, -0.05)
  vi <- c(0.04, 0.09, 0.16)
  mu <- matrix(
    c(
      0.00, 0.10, -0.10,
      0.05, 0.15, -0.02
    ),
    nrow  = 2,
    byrow = TRUE
  )
  tau <- matrix(
    c(
      0.10, 0.20, 0.05,
      0.05, 0.15, 0.08
    ),
    nrow  = 2,
    byrow = TRUE
  )

  mv_object <- brma.mv(
    yi                        = yi,
    V                         = diag(vi),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  uni_object <- brma.norm(
    yi                        = yi,
    vi                        = vi,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  make_setup <- function(object) {
    list(
      fit               = NULL,
      priors            = object[["priors"]],
      data              = object[["data"]],
      yi                = object[["data"]][["outcome"]][["yi"]],
      sei               = object[["data"]][["outcome"]][["sei"]],
      selection_sei     = object[["data"]][["outcome"]][["sei"]],
      K                 = length(yi),
      S                 = nrow(mu),
      mu                = mu,
      tau_within        = tau,
      tau_between       = matrix(0, nrow = nrow(mu), ncol = length(yi)),
      cluster           = NULL,
      weights           = NULL,
      data_hash         = .get_outcome_hash(object),
      is_weightfunction = FALSE,
      outcome_type      = "norm",
      effect_direction  = "positive",
      posterior_samples = matrix(numeric(0), nrow = nrow(mu), ncol = 0L)
    )
  }

  mv_setup  <- make_setup(mv_object)
  uni_setup <- make_setup(uni_object)

  expect_false(.known_v_estimate_target_uses_backend(mv_object[["data"]]))
  expect_false(.known_v_estimate_target_uses_schur_conditioning(mv_object[["data"]]))
  expect_equal(
    .log_lik_estimate_from_setup(mv_setup),
    .log_lik_estimate_from_setup(uni_setup),
    tolerance = 1e-12
  )

  target <- .estimate_log_lik_target_metadata(
    setup     = mv_setup,
    data_hash = .get_outcome_hash(mv_object)
  )
  expect_equal(target[["target"]], "factorized_estimate")
  expect_false(target[["known_v_estimate_backend"]])
  expect_false(target[["known_v_schur"]])
})


test_that("diagonal brma.mv random estimate target remains factorized", {

  object <- brma.mv(
    yi                         = c(0.10, 0.20, -0.05),
    V                          = diag(c(0.04, 0.09, 0.16)),
    random                     = ~ 1 | study,
    data                       = data.frame(study = c("a", "b", "c")),
    measure                    = "GEN",
    prior_unit_information_sd  = 1,
    marginalize_estimate_level = FALSE,
    only_priors                = TRUE
  )
  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = object[["data"]][["outcome"]][["yi"]],
    sei               = object[["data"]][["outcome"]][["sei"]],
    selection_sei     = object[["data"]][["outcome"]][["sei"]],
    K                 = 3L,
    S                 = 1L,
    mu                = matrix(c(0, 0, 0), nrow = 1),
    tau_within        = matrix(0, nrow = 1, ncol = 3),
    tau_between       = matrix(0, nrow = 1, ncol = 3),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(object),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 1, ncol = 0L)
  )

  target <- .estimate_log_lik_target_metadata(
    setup     = setup,
    data_hash = .get_outcome_hash(object)
  )
  expect_true(target[["known_v_estimate_backend"]])
  expect_false(target[["known_v_schur"]])
  expect_equal(target[["target"]], "factorized_estimate")
  expect_equal(target[["dependency_component_sizes"]], rep(1L, 3L))
  expect_equal(target[["random_effects"]], "conditioned")
})


test_that("brma.mv omits latent factors for diagonal known V", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = diag(c(0.04, 0.09)),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_equal(
    attr(object[["data"]], "known_V_data")[["parameterization"]],
    "whitened"
  )
  expect_equal(
    attr(object[["data"]], "known_V_data")[["parameterization_requested"]],
    "auto"
  )

  fit_priors <- .create_fit_priors(object[["data"]], object[["priors"]])
  fit_data   <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax     <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_false("sampling_z" %in% names(fit_priors))
  expect_false("sampling_rank" %in% names(fit_data))
  expect_false("sampling_B" %in% names(fit_data))
  expect_false(grepl("sampling_dependency", syntax, fixed = TRUE))
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "test()")
  )
})


test_that("brma.mv auto selects exact known-V backend when feasible", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30),
    x     = c(0, 1, 2),
    study = c("s1", "s2", "s3")
  )
  V <- matrix(
    c(
      0.04, 0.02, 0.01,
      0.02, 0.09, 0.03,
      0.01, 0.03, 0.16
    ),
    nrow  = 3,
    byrow = TRUE
  )

  no_scale <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  no_scale_known_V <- attr(no_scale[["data"]], "known_V_data")
  expect_equal(no_scale_known_V[["parameterization_requested"]], "auto")
  expect_equal(no_scale_known_V[["parameterization"]], "whitened")

  random_scalar <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  random_scalar_known_V <- attr(random_scalar[["data"]], "known_V_data")
  expect_equal(random_scalar_known_V[["parameterization_requested"]], "auto")
  expect_equal(random_scalar_known_V[["parameterization"]], "whitened")
  expect_true(.data_has_marginalized_random_effects(random_scalar[["data"]]))
  expect_false(.marginalized_random_effects_row_varying(
    .data_marginalized_random_effects(random_scalar[["data"]])
  ))

  scale_small <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  scale_small_known_V <- attr(scale_small[["data"]], "known_V_data")
  expect_equal(scale_small_known_V[["parameterization_requested"]], "auto")
  expect_equal(scale_small_known_V[["parameterization"]], "block_mvn")

  old_max_block <- getOption("RoBMA.known_v_block_mvn_max_block_size", NULL)
  on.exit({
    options(RoBMA.known_v_block_mvn_max_block_size = old_max_block)
  }, add = TRUE)

  options(RoBMA.known_v_block_mvn_max_block_size = 3L)
  scale_boundary <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  expect_equal(
    attr(scale_boundary[["data"]], "known_V_data")[["parameterization"]],
    "block_mvn"
  )

  options(RoBMA.known_v_block_mvn_max_block_size = Inf)
  scale_inf <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  expect_equal(
    attr(scale_inf[["data"]], "known_V_data")[["parameterization"]],
    "block_mvn"
  )

  random_scale_small <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  random_scale_small_known_V <- attr(random_scale_small[["data"]], "known_V_data")
  expect_equal(random_scale_small_known_V[["parameterization_requested"]], "auto")
  expect_equal(random_scale_small_known_V[["parameterization"]], "block_mvn")
  expect_true(.data_has_marginalized_random_effects(random_scale_small[["data"]]))
  expect_true(.marginalized_random_effects_row_varying(
    .data_marginalized_random_effects(random_scale_small[["data"]])
  ))

  options(RoBMA.known_v_block_mvn_max_block_size = 2L)

  scale_large <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  scale_large_known_V <- attr(scale_large[["data"]], "known_V_data")
  expect_equal(scale_large_known_V[["parameterization_requested"]], "auto")
  expect_equal(scale_large_known_V[["parameterization"]], "latent")

  random_scale_large <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  random_scale_large_known_V <- attr(random_scale_large[["data"]], "known_V_data")
  expect_equal(random_scale_large_known_V[["parameterization_requested"]], "auto")
  expect_equal(random_scale_large_known_V[["parameterization"]], "latent")

  explicit_block <- brma.mv(
    yi                        = yi,
    V                         = V,
    scale                     = ~ x,
    data                      = dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  explicit_known_V <- attr(explicit_block[["data"]], "known_V_data")
  expect_equal(explicit_known_V[["parameterization_requested"]], "block_mvn")
  expect_equal(explicit_known_V[["parameterization"]], "block_mvn")

  options(RoBMA.known_v_block_mvn_max_block_size = NA_real_)
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = V,
      scale                     = ~ x,
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "single positive"
  )
})


test_that("brma.mv prepares whitened known-V backend", {

  V <- matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2)

  expect_warning(
    object <- brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = V,
      known_v_parameterization  = "whitened",
      known_v_residual_fraction = 0.50,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "disregarded"
  )

  known_V <- attr(object[["data"]], "known_V_data")
  expect_equal(known_V[["parameterization"]], "whitened")
  expect_true(known_V[["correlated"]])
  expect_equal(known_V[["rank"]], 0L)
  expect_equal(known_V[["residual_fraction_requested"]], 0.50)

  rotated_V <- known_V[["whitening_matrix"]] %*% V %*%
    t(known_V[["whitening_matrix"]])
  expect_equal(
    rotated_V,
    diag(known_V[["whitening_variance"]], nrow = 2),
    tolerance = 1e-10
  )
})


test_that("brma.mv auto known-V update preserves residual fraction", {

  old_options <- options(RoBMA.known_v_block_mvn_max_block_size = 1L)
  on.exit(options(old_options))

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"), levels = c("e1", "e2", "e3"))
  )
  V <- matrix(
    c(
      0.040, 0.018, 0.012,
      0.018, 0.090, 0.024,
      0.012, 0.024, 0.160
    ),
    nrow  = 3,
    byrow = TRUE
  )
  R <- diag(c(1, 4, 9))
  dimnames(R) <- list(levels(dat[["estimate"]]), levels(dat[["estimate"]]))

  expect_warning(
    object <- brma.mv(
      yi                          = yi,
      V                           = V,
      random                      = ~ 1 | estimate,
      R                           = R,
      Rscale                      = "none",
      data                        = dat,
      known_v_parameterization    = "auto",
      known_v_residual_fraction   = 0.35,
      measure                     = "GEN",
      prior_unit_information_sd   = 1,
      only_priors                 = TRUE
    ),
    NA
  )

  known_V <- attr(object[["data"]], "known_V_data")
  expect_equal(known_V[["parameterization_requested"]], "auto")
  expect_equal(known_V[["parameterization"]], "latent")
  expect_equal(known_V[["residual_fraction_requested"]], 0.35)

  data <- brma.mv(
    yi                        = yi,
    V                         = diag(c(0.04, 0.09, 0.16)),
    data                      = data.frame(yi = c(0.10, 0.20, 0.30)),
    known_v_parameterization  = "auto",
    known_v_residual_fraction = 0.40,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )[["data"]]
  term <- list(
    block_name           = "study",
    sd_parameter_names   = "tau",
    row_multiplier       = c(1, 4, 9),
    row_multiplier_name  = "known_r_multiplier"
  )

  updated_known_V <- attr(.known_v_auto_update_for_marginalized_random(
    data  = data,
    terms = list(term)
  ), "known_V_data")

  expect_equal(updated_known_V[["parameterization_requested"]], "auto")
  expect_equal(updated_known_V[["parameterization"]], "block_mvn")
  expect_equal(updated_known_V[["residual_fraction_requested"]], 0.40)
})


test_that("brma.mv marginalized known R contributes known-V row variance", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"), levels = c("e1", "e2", "e3"))
  )
  R <- diag(c(4, 9, 16))
  dimnames(R) <- list(levels(dat[["estimate"]]), levels(dat[["estimate"]]))

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(dat)),
    random                    = ~ 1 | estimate,
    R                         = R,
    Rscale                    = "none",
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  known_V  <- attr(object[["data"]], "known_V_data")
  term     <- .data_marginalized_random_effects(object[["data"]])[[1L]]
  sd_name  <- term[["sd_parameter_names"]]
  row_name <- term[["row_multiplier_name"]]
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax   <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_equal(known_V[["parameterization_requested"]], "auto")
  expect_equal(known_V[["parameterization"]], "block_mvn")
  expect_equal(unname(term[["row_multiplier"]]), c(4, 9, 16))
  expect_equal(unname(fit_data[[row_name]]), c(4, 9, 16))
  expect_match(
    .data_marginalized_random_variance_expression(object[["data"]]),
    paste0("pow\\(", sd_name, ",2\\) \\* ", row_name, "\\[i\\]")
  )
  expect_match(syntax, paste0("tau2_observed\\[i\\] = pow\\(", sd_name))
  expect_match(syntax, paste0(row_name, "\\[i\\]"))

  posterior <- matrix(
    0.20,
    nrow     = 1L,
    dimnames = list(NULL, sd_name)
  )
  expected_extra_variance <- matrix(
    0.20^2 * c(4, 9, 16),
    nrow = 1L
  )
  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = object[["data"]],
      posterior_samples = posterior
    ),
    expected_extra_variance
  )

  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = dat[["yi"]],
    sei               = object[["data"]][["outcome"]][["sei"]],
    selection_sei     = object[["data"]][["outcome"]][["sei"]],
    K                 = 3L,
    S                 = 1L,
    mu                = matrix(0, nrow = 1L, ncol = 3L),
    tau_within        = matrix(0, nrow = 1L, ncol = 3L),
    tau_between       = matrix(0, nrow = 1L, ncol = 3L),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(object),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = posterior
  )
  expected_log_lik <- matrix(
    stats::dnorm(
      dat[["yi"]],
      mean = 0,
      sd   = sqrt(0.01 + as.numeric(expected_extra_variance)),
      log  = TRUE
    ),
    nrow = 1L
  )
  parameters <- list(mu = 0)
  parameters[[sd_name]] <- 0.20

  expect_equal(
    .known_v_extra_variance_from_setup(setup),
    expected_extra_variance
  )
  expect_equal(
    .log_lik_estimate_from_setup(setup),
    expected_log_lik,
    tolerance = 1e-12
  )
  expect_equal(
    .marglik_known_v_extra_variance(
      parameters         = parameters,
      model_data         = object[["data"]],
      bridge_context     = NULL,
      tau_within_samples = matrix(0, nrow = 1L, ncol = 3L),
      is_random          = TRUE,
      K                  = 3L
    ),
    expected_extra_variance
  )
  expect_equal(
    .log_posterior(
      parameters                 = parameters,
      data                       = fit_data,
      is_mods                    = FALSE,
      is_scale                   = FALSE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "block_mvn",
      model_data                 = object[["data"]]
    ),
    sum(expected_log_lik),
    tolerance = 1e-12
  )

  target <- .estimate_log_lik_target_metadata(
    setup     = setup,
    data_hash = .get_outcome_hash(object)
  )
  expect_true(target[["known_r"]])
  expect_equal(target[["known_r_blocks"]], "estimate")
  expect_true(grepl("diagonal tau^2 row multipliers",
                    target[["known_r_semantics"]], fixed = TRUE))
  expect_equal(target[["random_effects"]], "none")
  expect_equal(target[["estimate_level_random"]], "marginalized")
})


test_that("brma.mv known R row multipliers enter known-V backends consistently", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"), levels = c("e1", "e2", "e3")),
    x        = c(0, 1, 2)
  )
  R <- diag(c(4, 9, 16))
  dimnames(R) <- list(levels(dat[["estimate"]]), levels(dat[["estimate"]]))
  V <- matrix(
    c(
      0.04, 0.01, 0.005,
      0.01, 0.09, 0.002,
      0.005, 0.002, 0.16
    ),
    nrow  = 3,
    byrow = TRUE
  )

  latent <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | estimate,
    R                         = R,
    Rscale                    = "none",
    data                      = dat,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  latent_term   <- .data_marginalized_random_effects(latent[["data"]])[[1L]]
  latent_syntax <- .create_model_syntax(latent[["data"]], latent[["priors"]])
  expect_match(latent_syntax, "sampling_var\\[i\\]")
  expect_match(latent_syntax, latent_term[["row_multiplier_name"]],
               fixed = TRUE)

  whitened <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | estimate,
    R                         = R,
    Rscale                    = "cor",
    data                      = dat,
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  whitened_term   <- .data_marginalized_random_effects(whitened[["data"]])[[1L]]
  whitened_syntax <- .create_model_syntax(whitened[["data"]], whitened[["priors"]])
  expect_false(.marginalized_random_effects_row_varying(list(whitened_term)))
  expect_equal(unname(whitened_term[["row_multiplier"]]), c(1, 1, 1))
  expect_match(whitened_syntax, "whitening_var\\[j\\]")
  expect_match(whitened_syntax, whitened_term[["row_multiplier_name"]],
               fixed = TRUE)

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = V,
      random                    = ~ 1 | estimate,
      R                         = R,
      Rscale                    = "none",
      data                      = dat,
      known_v_parameterization  = "whitened",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "row-varying marginalized estimate-level random effects"
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = diag(0.01, nrow(dat)),
      scale                     = ~ x,
      random                    = ~ 1 | estimate,
      R                         = R,
      Rscale                    = "none",
      data                      = dat,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "row-indexed external SD sources"
  )
})


test_that("brma.mv adds whitened known-V likelihood to fit syntax", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  fit_priors <- .create_fit_priors(object[["data"]], object[["priors"]])
  fit_data   <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax     <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_false("sampling_z" %in% names(fit_priors))
  expect_true(all(c("whitening_y", "whitening_var", "whitening_matrix") %in% names(fit_data)))
  expect_false("sampling_var" %in% names(fit_data))
  expect_false("sampling_B" %in% names(fit_data))
  expect_match(syntax, "mu_observed", fixed = TRUE)
  expect_match(syntax, "whitening_mu", fixed = TRUE)
  expect_match(syntax, "whitening_y", fixed = TRUE)
  expect_false(grepl("sampling_dependency", syntax, fixed = TRUE))
})


test_that("brma.mv prepares block-MVN known-V backend", {

  V <- matrix(
    c(
      0.04, 0.03, 0.00,
      0.03, 0.09, 0.00,
      0.00, 0.00, 0.16
    ),
    nrow = 3,
    byrow = TRUE
  )

  object <- brma.mv(
    yi                       = c(0.10, 0.20, 0.30),
    V                        = V,
    known_v_parameterization = "block_mvn",
    measure                  = "GEN",
    prior_unit_information_sd = 1,
    only_data                = TRUE
  )

  known_V <- attr(object[["data"]], "known_V_data")
  expect_equal(known_V[["parameterization"]], "block_mvn")
  expect_true(known_V[["correlated"]])
  expect_equal(known_V[["rank"]], 0L)
  expect_equal(known_V[["residual_variance"]], diag(V))
  expect_equal(known_V[["sampling_factor"]], chol(V))
  expect_null(known_V[["residual_fraction_requested"]])

  expect_length(known_V[["block_mvn_blocks"]], 2L)
  expect_equal(known_V[["block_mvn_blocks"]][[1]][["index"]], 1:2)
  expect_equal(
    known_V[["block_mvn_blocks"]][[1]][["v_lower"]],
    V[1:2, 1:2][lower.tri(V[1:2, 1:2], diag = TRUE)]
  )
  expect_equal(known_V[["block_mvn_blocks"]][[2]][["index"]], 3L)
})


test_that("brma.mv adds block-MVN known-V likelihood to fit syntax", {

  object <- brma.mv(
    yi                       = c(0.10, 0.20, 0.30),
    V                        = matrix(
      c(
        0.04, 0.03, 0.00,
        0.03, 0.09, 0.00,
        0.00, 0.00, 0.16
      ),
      nrow = 3,
      byrow = TRUE
    ),
    known_v_parameterization = "block_mvn",
    measure                  = "GEN",
    prior_unit_information_sd = 1,
    only_priors              = TRUE
  )

  fit_priors <- .create_fit_priors(object[["data"]], object[["priors"]])
  fit_data   <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax     <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_false("sampling_z" %in% names(fit_priors))
  expect_false(any(c("sampling_var", "sampling_B") %in% names(fit_data)))
  expect_false(any(c("whitening_y", "whitening_var", "whitening_matrix") %in% names(fit_data)))
  expect_true(all(c("yi", "known_v_var", "known_v_y_1", "known_v_lower_1") %in% names(fit_data)))
  expect_equal(fit_data[["known_v_y_1"]], c(0.10, 0.20))
  expect_equal(length(fit_data[["known_v_lower_1"]]), 3L)

  expect_match(syntax, "mu_observed", fixed = TRUE)
  expect_match(syntax, "tau2_observed", fixed = TRUE)
  expect_match(syntax, "dknown_v_mnorm", fixed = TRUE)
  expect_match(syntax, "known_v_y_1\\[1:2\\] ~ dknown_v_mnorm")
  expect_match(syntax, "yi\\[3\\] ~ dnorm\\(mu_observed\\[3\\]")
  expect_false(grepl("sampling_dependency", syntax, fixed = TRUE))
})


test_that("brma.mv omits scalar known-V data for all-correlated block-MVN fits", {

  object <- brma.mv(
    yi                       = c(0.10, 0.20),
    V                        = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
    known_v_parameterization = "block_mvn",
    measure                  = "GEN",
    prior_unit_information_sd = 1,
    only_priors              = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax   <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_false("yi" %in% names(fit_data))
  expect_false("known_v_var" %in% names(fit_data))
  expect_true(all(c("known_v_y_1", "known_v_lower_1") %in% names(fit_data)))
  expect_false(grepl("yi\\[", syntax))
  expect_match(syntax, "known_v_y_1\\[1:2\\] ~ dknown_v_mnorm")
})


test_that("brma.mv block-MVN packing handles mixed-sign covariance blocks", {

  R <- matrix(
    c(
      1.00, -0.65,  0.45,
     -0.65,  1.00, -0.40,
      0.45, -0.40,  1.00
    ),
    nrow = 3,
    byrow = TRUE
  )
  V_block <- diag(c(0.18, 0.24, 0.30)) %*% R %*% diag(c(0.18, 0.24, 0.30))
  V <- rbind(
    cbind(V_block, matrix(0, 3, 1)),
    cbind(matrix(0, 1, 3), matrix(0.025, 1, 1))
  )

  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30, 0.40),
    V                         = V,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  blocks <- attr(object[["data"]], "known_V_data")[["block_mvn_blocks"]]

  expect_length(blocks, 2L)
  expect_equal(blocks[[1]][["index"]], 1:3)
  expect_equal(blocks[[1]][["v_lower"]], V_block[lower.tri(V_block, diag = TRUE)])
  expect_equal(blocks[[2]][["index"]], 4L)
})


test_that("brma.mv block-MVN accepts near-singular positive known-V blocks", {

  rho <- 0.99995
  R   <- matrix(rho, nrow = 3, ncol = 3)
  diag(R) <- 1
  V <- diag(c(0.10, 0.20, 0.40)) %*% R %*% diag(c(0.10, 0.20, 0.40))

  object <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = V,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )

  known_V <- attr(object[["data"]], "known_V_data")

  expect_equal(known_V[["block_mvn_blocks"]][[1]][["index"]], 1:3)
  expect_true(all(is.finite(known_V[["block_mvn_blocks"]][[1]][["v_lower"]])))
  expect_gt(known_V[["diagnostics"]][["min_eigenvalue"]], 0)

  rho_tight <- 1 - 1e-10
  V_tight <- 0.04 * matrix(c(
    1, rho_tight,
    rho_tight, 1
  ), nrow = 2, byrow = TRUE)

  tight_object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = V_tight,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  tight_known_V <- attr(tight_object[["data"]], "known_V_data")

  expect_equal(tight_known_V[["block_mvn_blocks"]][[1]][["index"]], 1:2)
  expect_gt(tight_known_V[["diagnostics"]][["min_eigenvalue"]], 0)
})


test_that("brma.mv allows block-MVN known-V scale models", {

  dat <- data.frame(
    yi = c(0.10, 0.20, 0.30),
    x  = c(0, 1, 2)
  )

  expect_silent(
    object <- brma.mv(
      yi                       = yi,
      V                        = matrix(
        c(
          0.04, 0.03, 0.00,
          0.03, 0.09, 0.00,
          0.00, 0.00, 0.16
        ),
        nrow = 3,
        byrow = TRUE
      ),
      scale                    = ~ x,
      data                     = dat,
      known_v_parameterization = "block_mvn",
      measure                  = "GEN",
      prior_unit_information_sd = 1,
      only_priors              = TRUE
    )
  )

  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])
  expect_match(syntax, "tau\\[i\\] = exp\\(log_tau\\[i\\]\\)")
  expect_match(syntax, "tau2_observed\\[i\\] = pow\\(tau\\[i\\],2\\)")
  expect_match(syntax, "dknown_v_mnorm", fixed = TRUE)

  random_object <- brma.mv(
    yi                       = yi,
    V                        = matrix(
      c(
        0.04, 0.03, 0.00,
        0.03, 0.09, 0.00,
        0.00, 0.00, 0.16
      ),
      nrow = 3,
      byrow = TRUE
    ),
    scale                    = ~ x,
    random                   = ~ 1 | study,
    data                     = data.frame(dat, study = c("a", "a", "b")),
    known_v_parameterization = "block_mvn",
    measure                  = "GEN",
    prior_unit_information_sd = 1,
    only_priors              = TRUE
  )
  random_syntax <- .create_model_syntax(random_object[["data"]], random_object[["priors"]])
  random_design <- .fitted_formula_design(random_object, "mu", required = TRUE)
  expect_match(random_syntax, "tau\\[i\\] = exp\\(log_tau\\[i\\]\\)")
  expect_match(random_syntax, "tau2_observed\\[i\\] = 0")
  expect_match(random_syntax, "dknown_v_mnorm", fixed = TRUE)
  expect_equal(
    random_design[["random_effects"]][[1]][["sd_binding"]][["source"]][["shape"]],
    "row"
  )
})


test_that("brma.mv rejects unsupported whitened known-V combinations", {

  dat <- data.frame(
    yi = c(0.10, 0.20),
    x  = c(0, 1)
  )

  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
      scale                     = ~ x,
      data                      = dat,
      known_v_parameterization  = "whitened",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_data                 = TRUE
    ),
    "scale regression"
  )
  expect_error(
    .known_v_prepare(
      V                                   = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
      keep_rows                           = c(TRUE, TRUE),
      known_v_parameterization            = "whitened",
      known_v_residual_fraction           = NULL,
      known_v_residual_fraction_specified = FALSE,
      known_v_is_scale                    = TRUE
    ),
    "scale regression"
  )
})


test_that("known-V preparation can suppress duplicate singular warnings", {

  V <- matrix(0.04, nrow = 2L, ncol = 2L)

  expect_warning(
    .known_v_prepare(
      V                                   = V,
      keep_rows                           = c(TRUE, TRUE),
      known_v_parameterization            = "whitened",
      known_v_residual_fraction           = NULL,
      known_v_residual_fraction_specified = FALSE,
      warn_singular                       = TRUE
    ),
    "singular"
  )
  expect_silent(
    .known_v_prepare(
      V                                   = V,
      keep_rows                           = c(TRUE, TRUE),
      known_v_parameterization            = "whitened",
      known_v_residual_fraction           = NULL,
      known_v_residual_fraction_specified = FALSE,
      warn_singular                       = FALSE
    )
  )
})


test_that("brma.mv known-V models pass bridge marginal likelihood availability", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = diag(c(0.04, 0.09)),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_silent(
    .check_marglik_available(object, "test()")
  )
})


test_that("bridge SD-source spec expands row sources and attaches bounds", {

  posterior <- matrix(
    1:6,
    nrow = 2L,
    dimnames = list(NULL, c("tau[1]", "tau[2]", "sigma"))
  )

  spec <- .marglik_bridge_sd_source_spec(
    add_parameters = c("tau", "sigma"),
    fit            = posterior,
    K              = 2L
  )

  expect_equal(spec[["parameters"]], c("tau[1]", "tau[2]", "sigma"))
  expect_equal(names(spec[["bounds"]][["lb"]]), spec[["parameters"]])
  expect_equal(names(spec[["bounds"]][["ub"]]), spec[["parameters"]])
  expect_equal(unname(spec[["bounds"]][["lb"]]), c(0, 0, 0))
  expect_true(all(is.infinite(spec[["bounds"]][["ub"]])))

  empty <- .marglik_bridge_sd_source_spec(
    add_parameters = character(),
    fit            = posterior,
    K              = 2L
  )
  expect_equal(empty[["parameters"]], character())
  expect_null(empty[["bounds"]])
})


test_that("brma.mv known-V bridge log posterior matches exact normal targets", {

  V <- matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2)
  parameters <- list(mu = 0, tau = 0.10)

  block_mvn <- brma.mv(
    yi                       = c(0.10, 0.20),
    V                        = V,
    known_v_parameterization = "block_mvn",
    measure                  = "GEN",
    prior_unit_information_sd = 1,
    only_priors              = TRUE
  )
  block_fit_data <- .create_fit_data(block_mvn[["data"]], block_mvn[["priors"]])
  block_expected <- .marglik_mvn_log_density(
    y          = c(0.10, 0.20),
    mean       = c(0, 0),
    covariance = V + diag(parameters[["tau"]]^2, nrow = 2)
  )
  expect_equal(
    .log_posterior(
      parameters                 = parameters,
      data                       = block_fit_data,
      is_mods                    = FALSE,
      is_scale                   = FALSE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "block_mvn",
      model_data                 = block_mvn[["data"]]
    ),
    block_expected,
    tolerance = 1e-10
  )
  expect_error(
    .log_posterior(
      parameters                 = parameters,
      data                       = block_fit_data,
      is_mods                    = FALSE,
      is_scale                   = FALSE,
      is_multilevel              = FALSE,
      is_weights                 = TRUE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "block_mvn",
      model_data                 = block_mvn[["data"]]
    ),
    "Known-V bridge likelihoods"
  )
  block_setup <- list(
    fit               = NULL,
    priors            = block_mvn[["priors"]],
    data              = block_mvn[["data"]],
    yi                = block_mvn[["data"]][["outcome"]][["yi"]],
    sei               = block_mvn[["data"]][["outcome"]][["sei"]],
    selection_sei     = block_mvn[["data"]][["outcome"]][["sei"]],
    K                 = 2L,
    S                 = 1L,
    mu                = matrix(c(0, 0), nrow = 1),
    tau_within        = matrix(c(parameters[["tau"]], parameters[["tau"]]),
                               nrow = 1),
    tau_between       = matrix(0, nrow = 1, ncol = 2),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(block_mvn),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 1, ncol = 0L)
  )
  expect_equal(
    .log_lik_estimate_sum_from_setup(block_setup),
    block_expected,
    tolerance = 1e-8
  )

  whitened <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = V,
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  whitened_fit_data <- .create_fit_data(whitened[["data"]], whitened[["priors"]])
  expect_equal(
    .log_posterior(
      parameters                 = parameters,
      data                       = whitened_fit_data,
      is_mods                    = FALSE,
      is_scale                   = FALSE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "whitened",
      model_data                 = whitened[["data"]]
    ),
    block_expected,
    tolerance = 1e-10
  )

  rank_one_sei        <- c(0.20, 0.30, 0.40)
  rank_one_V          <- tcrossprod(rank_one_sei)
  rank_one_parameters <- list(mu = 0, tau = 0.10)
  rank_one_expected   <- .marglik_mvn_log_density(
    y          = c(0.10, 0.20, -0.05),
    mean       = rep(0, 3),
    covariance = rank_one_V + diag(rank_one_parameters[["tau"]]^2, nrow = 3)
  )

  expect_warning(
    rank_one_block <- brma.mv(
      yi                        = c(0.10, 0.20, -0.05),
      V                         = rank_one_V,
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )
  expect_equal(
    .log_posterior(
      parameters                 = rank_one_parameters,
      data                       = .create_fit_data(
        rank_one_block[["data"]],
        rank_one_block[["priors"]]
      ),
      is_mods                    = FALSE,
      is_scale                   = FALSE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "block_mvn",
      model_data                 = rank_one_block[["data"]]
    ),
    rank_one_expected,
    tolerance = 1e-10
  )

  expect_warning(
    rank_one_whitened <- brma.mv(
      yi                        = c(0.10, 0.20, -0.05),
      V                         = rank_one_V,
      known_v_parameterization  = "whitened",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )
  expect_equal(
    .log_posterior(
      parameters                 = rank_one_parameters,
      data                       = .create_fit_data(
        rank_one_whitened[["data"]],
        rank_one_whitened[["priors"]]
      ),
      is_mods                    = FALSE,
      is_scale                   = FALSE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "whitened",
      model_data                 = rank_one_whitened[["data"]]
    ),
    rank_one_expected,
    tolerance = 1e-10
  )

  dat_scale <- data.frame(
    yi = c(0.10, 0.20, -0.05),
    x  = c(-1, 0, 1)
  )
  V_scale <- matrix(
    c(
      0.05, 0.015, 0.010,
      0.015, 0.08, 0.020,
      0.010, 0.020, 0.07
    ),
    nrow  = 3,
    byrow = TRUE
  )
  scale_tau <- c(0.04, 0.10, 0.18)
  scale_parameters <- list(
    mu      = 0,
    log_tau = log(scale_tau)
  )
  scale_block <- brma.mv(
    yi                        = yi,
    V                         = V_scale,
    scale                     = ~ x,
    data                      = dat_scale,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  scale_expected <- .marglik_mvn_log_density(
    y          = dat_scale[["yi"]],
    mean       = rep(0, 3),
    covariance = V_scale + diag(scale_tau^2, nrow = 3)
  )
  expect_equal(
    .log_posterior(
      parameters                 = scale_parameters,
      data                       = .create_fit_data(scale_block[["data"]], scale_block[["priors"]]),
      is_mods                    = FALSE,
      is_scale                   = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "block_mvn",
      model_data                 = scale_block[["data"]]
    ),
    scale_expected,
    tolerance = 1e-10
  )

  dat_random <- data.frame(
    yi    = c(0.10, 0.20),
    study = c("s1", "s2")
  )
  random_sd <- 0.12
  random_parameters <- list(mu = c(0, 0))
  random_parameters[["mu__xREx__study_intercept"]] <- random_sd
  random_expected <- .marglik_mvn_log_density(
    y          = dat_random[["yi"]],
    mean       = c(0, 0),
    covariance = V + diag(random_sd^2, nrow = 2)
  )

  random_block <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat_random,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  expect_equal(
    .log_posterior(
      parameters                 = random_parameters,
      data                       = .create_fit_data(random_block[["data"]], random_block[["priors"]]),
      is_mods                    = FALSE,
      is_scale                   = FALSE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "block_mvn",
      model_data                 = random_block[["data"]]
    ),
    random_expected,
    tolerance = 1e-10
  )
  expect_true(.marglik_needs_bridge_context(random_block[["data"]]))

  random_whitened <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat_random,
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  row_parameters <- list(mu = c(0, 0))
  attr(row_parameters, "posterior_samples") <- matrix(
    random_sd,
    nrow     = 1L,
    dimnames = list(NULL, "mu__xREx__study_intercept")
  )
  expect_equal(
    unname(.parameters_as_sample_matrix(row_parameters)[1L, "mu__xREx__study_intercept"]),
    random_sd
  )
  expect_equal(
    .log_posterior(
      parameters                 = row_parameters,
      data                       = .create_fit_data(random_whitened[["data"]], random_whitened[["priors"]]),
      is_mods                    = FALSE,
      is_scale                   = FALSE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "whitened",
      model_data                 = random_whitened[["data"]]
    ),
    random_expected,
    tolerance = 1e-10
  )

  allocation_dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30, 0.40),
    study    = c("s1", "s1", "s2", "s2"),
    estimate = paste0("e", 1:4),
    x        = c(0, 1, 0, 1)
  )
  allocation_V <- diag(rep(0.04, 4))
  allocation_tau <- c(0.10, 0.20, 0.30, 0.40)
  allocation_block <- brma.mv(
    yi                        = yi,
    V                         = allocation_V,
    scale                     = ~ x,
    random                    = ~ 1 | study / estimate,
    data                      = allocation_dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  allocation_term <- .data_marginalized_random_effects(
    allocation_block[["data"]]
  )[[1L]]
  allocation_factor <- .marginalized_random_effect_allocation_factors(
    allocation_term
  )[[1L]]
  allocation_weight_name <- allocation_factor[["weight_name"]]
  allocation_bridge_context <- structure(
    list(nodes = c(
      stats::setNames(allocation_tau, paste0("tau[", seq_along(allocation_tau), "]")),
      stats::setNames(
        c(0.25, 0.75),
        paste0(allocation_weight_name, "[", 1:2, "]")
      )
    )),
    class = c("BayesTools_bridge_context", "list")
  )
  allocation_parameters <- list(mu = 0)
  allocation_posterior <- .marglik_bridge_posterior_samples(
    parameters     = allocation_parameters,
    bridge_context = allocation_bridge_context
  )
  expect_equal(
    unname(allocation_posterior[1L, paste0(allocation_weight_name, "[1]")]),
    0.25
  )
  allocation_expected <- .marglik_mvn_log_density(
    y          = allocation_dat[["yi"]],
    mean       = rep(0, 4),
    covariance = allocation_V + diag(allocation_tau^2 * 0.25, nrow = 4)
  )
  expect_true(.marglik_needs_bridge_context(allocation_block[["data"]]))
  expect_equal(
    .log_posterior(
      parameters                 = allocation_parameters,
      data                       = .create_fit_data(allocation_block[["data"]], allocation_block[["priors"]]),
      is_mods                    = FALSE,
      is_scale                   = TRUE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "block_mvn",
      model_data                 = allocation_block[["data"]],
      bridge_context             = allocation_bridge_context
    ),
    allocation_expected,
    tolerance = 1e-10
  )

  sd_component_dat <- transform(
    allocation_dat,
    criterion = c("c1", "c2", "c1", "c2")
  )
  sd_component_block <- brma.mv(
    yi                        = yi,
    V                         = allocation_V,
    random                    = list(~ 1 | study / estimate, ~ 1 | criterion),
    data                      = sd_component_dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  sd_component_terms <- .data_marginalized_random_effects(
    sd_component_block[["data"]]
  )
  sd_component_term       <- sd_component_terms[[1L]]
  sd_component_binding    <- sd_component_term[["sd_binding"]]
  sd_component_allocation <- sd_component_binding[["allocations"]][[1L]]
  sd_component_allocation[["target"]]               <- "sd_component"
  sd_component_allocation[["parent_factors"]]       <- sd_component_allocation[["factors"]]
  sd_component_allocation[["leaf_index_by_column"]] <- 1L
  sd_component_allocation[["weight_name"]]          <- "mu__xRE_ALLOCx_leaf__weight"
  sd_component_allocation[["scale"]]                <- "mean_variance"
  sd_component_allocation[["n_targets"]]            <- 3L
  sd_component_binding[["allocations"]][[1L]]       <- sd_component_allocation
  sd_component_term[["sd_parameter_names"]]         <- NA_character_
  sd_component_term[["sd_binding"]]                 <- sd_component_binding
  sd_component_terms[[1L]]      <- sd_component_term
  sd_component_location         <- sd_component_block[["data"]][["location"]]
  attr(sd_component_location, "marginalized_random_effects") <- sd_component_terms
  sd_component_block[["data"]][["location"]] <- sd_component_location

  parent_factors <- sd_component_allocation[["parent_factors"]]
  parent_nodes   <- stats::setNames(
    c(0.25, 0.40),
    vapply(parent_factors, function(factor) {
      paste0(factor[["weight_name"]], "[", factor[["index"]], "]")
    }, character(1))
  )
  sd_component_context <- structure(
    list(nodes = c(
      stats::setNames(
        0.50,
        sd_component_allocation[["source"]][["name"]]
      ),
      parent_nodes,
      stats::setNames(
        0.20,
        paste0(sd_component_allocation[["weight_name"]], "[1]")
      )
    )),
    class = c("BayesTools_bridge_context", "list")
  )
  sd_component_parameters <- list(mu = 0)
  sd_component_posterior  <- .marglik_bridge_posterior_samples(
    parameters     = sd_component_parameters,
    bridge_context = sd_component_context
  )
  sd_component_variance   <- 0.50^2 * 0.25 * 0.40 * (3 * 0.20)
  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = sd_component_block[["data"]],
      posterior_samples = sd_component_posterior
    ),
    matrix(sd_component_variance, nrow = 1L, ncol = 4L)
  )
  sd_component_expected   <- .marglik_mvn_log_density(
    y          = sd_component_dat[["yi"]],
    mean       = rep(0, 4),
    covariance = allocation_V + diag(sd_component_variance, nrow = 4)
  )
  expect_equal(
    .log_posterior(
      parameters                 = sd_component_parameters,
      data                       = .create_fit_data(sd_component_block[["data"]], sd_component_block[["priors"]]),
      is_mods                    = FALSE,
      is_scale                   = FALSE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      known_v_parameterization   = "block_mvn",
      model_data                 = sd_component_block[["data"]],
      bridge_context             = sd_component_context
    ),
    sd_component_expected,
    tolerance = 1e-10
  )
})


test_that("brma.mv preserves tiny nonzero known-V covariance", {

  V <- matrix(c(0.04, 1e-12, 1e-12, 0.09), nrow = 2)

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = V,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  known_V <- attr(object[["data"]], "known_V_data")

  expect_length(known_V[["block_indices"]], 1L)
  expect_gt(known_V[["rank"]], 0L)
  expect_equal(
    diag(known_V[["residual_variance"]], nrow = 2) + tcrossprod(known_V[["B"]]),
    V,
    tolerance = 1e-10
  )
})


test_that("brma.mv known-V estimate log-likelihood uses Schur conditional target", {

  V <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2)
  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = V,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  mu <- matrix(
    c(
      0.00, 0.10,
      0.05, 0.15
    ),
    nrow  = 2,
    byrow = TRUE
  )
  tau <- matrix(
    c(
      0.10, 0.20,
      0.05, 0.15
    ),
    nrow  = 2,
    byrow = TRUE
  )
  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = object[["data"]][["outcome"]][["yi"]],
    sei               = object[["data"]][["outcome"]][["sei"]],
    selection_sei     = object[["data"]][["outcome"]][["sei"]],
    K                 = 2L,
    S                 = 2L,
    mu                = mu,
    tau_within        = tau,
    tau_between       = matrix(0, nrow = 2, ncol = 2),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(object),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 2, ncol = 0L)
  )

  log_lik <- .log_lik_estimate_from_setup(setup)
  expected <- matrix(NA_real_, nrow = 2, ncol = 2)
  for (s in seq_len(2)) {
    covariance <- V + diag(tau[s, ]^2, nrow = 2)
    for (i in seq_len(2)) {
      j <- setdiff(seq_len(2), i)
      conditional_mean <- mu[s, i] +
        covariance[i, j] / covariance[j, j] * (setup[["yi"]][j] - mu[s, j])
      conditional_variance <- covariance[i, i] -
        covariance[i, j]^2 / covariance[j, j]
      expected[s, i] <- stats::dnorm(
        setup[["yi"]][i],
        mean = conditional_mean,
        sd   = sqrt(conditional_variance),
        log  = TRUE
      )
    }
  }

  expect_equal(log_lik, expected, tolerance = 1e-12)
  expected_joint <- vapply(seq_len(2), function(s) {
    .marglik_mvn_log_density(
      y          = setup[["yi"]],
      mean       = mu[s, ],
      covariance = V + diag(tau[s, ]^2, nrow = 2)
    )
  }, numeric(1))
  expect_equal(
    .log_lik_estimate_sum_from_setup(setup),
    expected_joint,
    tolerance = 1e-12
  )

  setup_negative <- setup
  setup_negative[["effect_direction"]] <- "negative"
  log_lik_negative <- .log_lik_estimate_from_setup(setup_negative)
  expected_negative <- matrix(NA_real_, nrow = 2, ncol = 2)
  for (s in seq_len(2)) {
    covariance <- V + diag(tau[s, ]^2, nrow = 2)
    for (i in seq_len(2)) {
      j <- setdiff(seq_len(2), i)
      conditional_mean <- -mu[s, i] +
        covariance[i, j] / covariance[j, j] * (-setup[["yi"]][j] + mu[s, j])
      conditional_variance <- covariance[i, i] -
        covariance[i, j]^2 / covariance[j, j]
      expected_negative[s, i] <- stats::dnorm(
        -setup[["yi"]][i],
        mean = conditional_mean,
        sd   = sqrt(conditional_variance),
        log  = TRUE
      )
    }
  }
  expect_equal(log_lik_negative, expected_negative, tolerance = 1e-12)

  target <- .estimate_log_lik_target_metadata(
    setup     = setup,
    data_hash = .get_outcome_hash(object)
  )
  expect_true(target[["known_v_estimate_backend"]])
  expect_true(target[["known_v_schur"]])
  expect_equal(target[["target"]], "known_v_estimate")
})


test_that("known-V dependency blocks are canonicalized from V metadata", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20, -0.05),
    V                         = matrix(
      c(
        0.04, 0.01, 0.00,
        0.01, 0.09, 0.00,
        0.00, 0.00, 0.16
      ),
      nrow = 3,
      byrow = TRUE
    ),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  data_missing_blocks <- object[["data"]]
  known_V             <- .data_known_v_data(data_missing_blocks)
  known_V[["block_indices"]] <- NULL
  attr(data_missing_blocks, "known_V_data") <- known_V

  blocks <- .known_v_dependency_blocks(data_missing_blocks, K = 3L)
  expect_equal(blocks, list(1:2, 3L))

  data_bad_blocks <- object[["data"]]
  known_V         <- .data_known_v_data(data_bad_blocks)
  known_V[["block_indices"]] <- list(c(1L, 1L), 3L)
  attr(data_bad_blocks, "known_V_data") <- known_V
  expect_error(
    .known_v_dependency_blocks(data_bad_blocks, K = 3L),
    "partition"
  )
})


test_that("singular PSD known-V Cholesky targets fail with targeted messages", {

  sei <- c(0.20, 0.30, 0.40)
  V   <- tcrossprod(sei)

  expect_warning(
    object <- brma.mv(
      yi                        = c(0.10, 0.20, -0.05),
      V                         = V,
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )

  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = object[["data"]][["outcome"]][["yi"]],
    sei               = object[["data"]][["outcome"]][["sei"]],
    selection_sei     = object[["data"]][["outcome"]][["sei"]],
    K                 = 3L,
    S                 = 1L,
    mu                = matrix(c(0, 0, 0), nrow = 1),
    tau_within        = matrix(0, nrow = 1, ncol = 3),
    tau_between       = matrix(0, nrow = 1, ncol = 3),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(object),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 1, ncol = 0L)
  )

  expect_error(
    .log_lik_estimate_from_setup(setup),
    "positive semidefinite"
  )
  expect_error(
    .marglik_mvn_log_density(
      y          = c(0.10, 0.20, -0.05),
      mean       = rep(0, 3),
      covariance = V
    ),
    "positive semidefinite"
  )
})


test_that("brma.mv Schur estimate target matches partitioned MVN conditionals", {

  V <- matrix(c(
    0.10, 0.03, 0.02,
    0.03, 0.12, 0.04,
    0.02, 0.04, 0.11
  ), nrow = 3)
  yi <- c(0.10, -0.05, 0.20)
  mu <- matrix(c(0.00, -0.02, 0.10), nrow = 1)
  tau <- matrix(c(0.05, 0.10, 0.08), nrow = 1)

  object <- brma.mv(
    yi                        = yi,
    V                         = V,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = yi,
    sei               = object[["data"]][["outcome"]][["sei"]],
    selection_sei     = object[["data"]][["outcome"]][["sei"]],
    K                 = 3L,
    S                 = 1L,
    mu                = mu,
    tau_within        = tau,
    tau_between       = matrix(0, nrow = 1, ncol = 3),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(object),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 1, ncol = 0L)
  )

  covariance <- V + diag(as.numeric(tau)^2, nrow = 3)
  expected   <- matrix(NA_real_, nrow = 1, ncol = 3)
  for (i in seq_len(3)) {
    j <- setdiff(seq_len(3), i)
    conditional_mean <- mu[1, i] +
      covariance[i, j, drop = FALSE] %*%
      solve(covariance[j, j, drop = FALSE], yi[j] - mu[1, j])
    conditional_variance <- covariance[i, i] -
      covariance[i, j, drop = FALSE] %*%
      solve(covariance[j, j, drop = FALSE],
            covariance[j, i, drop = FALSE])
    expected[1, i] <- stats::dnorm(
      yi[i],
      mean = conditional_mean,
      sd   = sqrt(conditional_variance),
      log  = TRUE
    )
  }

  expect_true(.known_v_estimate_target_uses_schur_conditioning(object[["data"]]))
  expect_equal(.log_lik_estimate_from_setup(setup), expected,
               tolerance = 1e-12)
})


test_that("brma.mv known-V target hash includes covariance structure", {

  object_1 <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.00, 0.00, 0.09), nrow = 2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  object_2 <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_false(identical(
    .get_outcome_hash(object_1),
    .get_outcome_hash(object_2)
  ))
})


test_that("brma.mv allows estimate log-likelihood diagnostics for correlated V", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_silent(
    .check_log_lik_target_available(object, "estimate", "test()")
  )
  expect_error(
    .check_log_lik_target_available(object, "cluster", "test()"),
    "unit = 'cluster'"
  )
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "add_waic()")
  )
})


test_that("brma.mv allows estimate log-likelihood diagnostics for correlated whitened V", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_silent(
    .check_log_lik_target_available(object, "estimate", "test()")
  )
  expect_error(
    .check_log_lik_target_available(object, "cluster", "test()"),
    "unit = 'cluster'"
  )
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "add_waic()")
  )
})


test_that("brma.mv regplot PI/SI no longer use the old feature guard", {

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    mods                      = ~ x,
    data                      = data.frame(
      yi = c(0.10, 0.20, 0.30, 0.40),
      x  = c(0, 1, 0, 1)
    ),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_error(
    regplot(object, pi = TRUE, as_data = TRUE),
    "Posterior samples are required"
  )
  expect_error(
    regplot(object, si = TRUE, as_data = TRUE),
    "Posterior samples are required"
  )
})
