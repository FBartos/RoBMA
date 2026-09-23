test_that("random semantic point hypotheses reject structural atoms", {

  source_prior <- BayesTools::prior_mixture(
    prior_list = list(
      BayesTools::prior("spike", parameters = list(location = 0.5)),
      BayesTools::prior("beta", parameters = list(alpha = 1, beta = 1))
    ),
    is_null = c(TRUE, FALSE)
  )
  prior_values <- c(rep(0.5, 50), seq(0.51, 1, length.out = 50))
  posterior_values <- c(rep(0.5, 25), seq(0.51, 1, length.out = 75))
  expect_equal(mean(prior_values == 0.5), 0.50)
  expect_equal(mean(posterior_values == 0.5), 0.25)

  selected <- function(values) {
    list(
      entry = list(term = "variance proportion"),
      spec = list(
        quantity     = "var_prop",
        source_parameter = "rho",
        label            = "total: tau2_prop(study)"
      ),
      samples      = matrix(values, ncol = 1L),
      prior        = BayesTools::prior("beta", list(alpha = 1, beta = 1)),
      source_prior = source_prior
    )
  }
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(
        object, parameter, standardized_coefficients, prior = FALSE, ...) {
      if (prior) selected(prior_values) else selected(posterior_values)
    },
    .package = "RoBMA"
  )

  point <- BayesTools::hypothesis_parse("theta = 0.5")
  for (density_method in c("KDE", "normal")) {
    expect_error(
      .hypothesis_brma_random(
        object                    = list(),
        parameter                 = "theta",
        hypothesis                = point,
        standardized_coefficients = FALSE,
        conditional               = FALSE,
        logBF                     = FALSE,
        BF01                      = FALSE,
        seed                      = 1,
        density_method            = density_method,
        n_samples                 = 100,
        columns                   = "default"
      ),
      "induced prior/posterior contains a point mass",
      fixed = TRUE
    )
  }

  region <- .hypothesis_brma_random(
    object                    = list(),
    parameter                 = "theta",
    hypothesis                = BayesTools::hypothesis_parse(
      "theta <= 0.75 vs theta > 0.75"
    ),
    standardized_coefficients = FALSE,
    conditional               = FALSE,
    logBF                     = FALSE,
    BF01                      = FALSE,
    seed                      = 1,
    density_method            = "KDE",
    n_samples                 = 100,
    columns                   = "default"
  )
  expect_s3_class(region, "BayesTools_hypothesis_BF")
  expect_true(is.finite(attr(region, "raw_BF")))
})


test_that("hypothesis discovery suppresses atomic random point routes", {

  source_prior <- BayesTools::prior_mixture(
    prior_list = list(
      BayesTools::prior("spike", parameters = list(location = 0.5)),
      BayesTools::prior("beta", parameters = list(alpha = 1, beta = 1))
    ),
    is_null = c(TRUE, FALSE)
  )
  fit <- list()
  attr(fit, "prior_list") <- list(rho = source_prior)
  object <- structure(list(fit = fit), class = "brma")
  specs <- data.frame(
    parameter         = "theta",
    label             = "total: tau2_prop(study)",
    quantity      = "var_prop",
    formula_parameter = "tau",
    block             = NA_character_,
    grouping          = "study",
    structure         = NA_character_,
    allocation        = "total",
    random_component  = "study",
    source_type       = "identity",
    source_parameter  = "rho",
    source_prior_name = "rho",
    source_transform  = "identity",
    source_scale      = 1,
    stringsAsFactors  = FALSE
  )
  specs[["display_transform"]] <- I(list(list(type = "identity")))
  testthat::local_mocked_bindings(
    .brma_parameter_catalog = function(object) {
      data.frame(
        alias      = "total: tau2_prop(study)",
        parameter  = "theta",
        component  = "random",
        term       = "variance proportion",
        stringsAsFactors = FALSE
      )
    },
    .brma_parameter_catalog_metadata = function(object) list(entries = NULL),
    .brma_random_parameter_bundle = function(object, ...) {
      list(
        samples = matrix(seq(0.1, 0.9, length.out = 20), ncol = 1L),
        specs   = specs,
        priors  = list(theta = BayesTools::prior(
          "beta", list(alpha = 1, beta = 1)
        ))
      )
    },
    .package = "RoBMA"
  )

  out <- hypothesis_quantities(object)
  expect_false(out[["point_test"]])
  expect_identical(out[["point_test_methods"]], "")
  expect_true(out[["direction_test"]])
  expect_identical(out[["direction_test_methods"]], "KDE, normal")
  expect_match(out[["reason"]], "induced prior/posterior contains a point mass")
})


test_that("realized allocation gates define aggregate and proportion atoms", {

  gate_prior <- function(indicator, probability) {
    out <- BayesTools::prior(
      "spike",
      parameters = list(location = probability)
    )
    attr(out, "random_allocation_indicator") <- indicator
    out
  }
  fit <- list()
  attr(fit, "prior_list") <- list(
    gate_study = gate_prior("gate_study", 0.5),
    gate_drug  = gate_prior("gate_drug", 0.5)
  )
  object <- list(fit = fit)
  allocation <- list(
    scale          = "total_variance",
    n_targets      = 2L,
    inclusion      = list(
      list(index = 1L, indicator_name = "gate_study"),
      list(index = 2L, indicator_name = "gate_drug")
    ),
    parent_factors = list()
  )
  raw_samples <- cbind(
    gate_study = c(0, 1, 0, 1),
    gate_drug  = c(0, 0, 1, 1)
  )

  total_selected <- list(
    spec = list(quantity = "sd_total", allocation_index = NA_integer_),
    allocation_definition = allocation
  )
  total_metadata <-
    .brma_random_parameter_allocation_gate_metadata(total_selected)
  total_state <- .brma_random_parameter_allocation_gate_state(
    total_metadata,
    raw_samples
  )
  total_prior <- .brma_random_parameter_allocation_gate_prior(
    object,
    total_metadata
  )
  expect_identical(total_state[["point_zero"]], c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(total_state[["continuous"]], c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(total_prior[["points"]], data.frame(x = 0, p = 0.25))
  expect_equal(total_prior[["continuous_mass"]], 0.75)

  proportion_selected <- list(
    spec = list(quantity = "var_prop", allocation_index = 1L),
    allocation_definition = allocation
  )
  proportion_metadata <-
    .brma_random_parameter_allocation_gate_metadata(proportion_selected)
  proportion_state <- .brma_random_parameter_allocation_gate_state(
    proportion_metadata,
    raw_samples
  )
  proportion_prior <- .brma_random_parameter_allocation_gate_prior(
    object,
    proportion_metadata
  )
  expect_identical(proportion_state[["defined"]], c(FALSE, TRUE, TRUE, TRUE))
  expect_identical(proportion_state[["point_zero"]], c(FALSE, FALSE, TRUE, FALSE))
  expect_identical(proportion_state[["point_one"]], c(FALSE, TRUE, FALSE, FALSE))
  expect_identical(proportion_state[["continuous"]], c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(
    proportion_prior[["points"]],
    data.frame(x = c(0, 1), p = c(1 / 3, 1 / 3))
  )
  expect_equal(proportion_prior[["continuous_mass"]], 1 / 3)

  reason <- .brma_random_parameter_point_test_reason(
    spec = list(
      quantity         = "var_prop",
      source_parameter = NA_character_,
      label            = "tau2_prop(study)"
    ),
    allocation_gate_prior = proportion_prior
  )
  expect_match(reason, "contains a point mass", fixed = TRUE)

  inherited_allocation <- allocation
  inherited_allocation[["inclusion"]] <- list()
  inherited_allocation[["parent_factors"]] <- list(
    list(inclusion_name = "gate_study")
  )
  inherited_selected <- proportion_selected
  inherited_selected[["allocation_definition"]] <- inherited_allocation
  inherited_metadata <-
    .brma_random_parameter_allocation_gate_metadata(inherited_selected)
  inherited_state <- .brma_random_parameter_allocation_gate_state(
    inherited_metadata,
    raw_samples
  )
  inherited_prior <- .brma_random_parameter_allocation_gate_prior(
    object,
    inherited_metadata
  )
  expect_identical(inherited_state[["defined"]], c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(
    inherited_state[["continuous"]],
    inherited_state[["defined"]]
  )
  expect_equal(nrow(inherited_prior[["points"]]), 0L)
  expect_equal(inherited_prior[["continuous_mass"]], 1)
})


test_that("allocation prior density preserves structural boundary masses", {

  samples <- c(0, 1, seq(0.01, 0.99, length.out = 2000L))
  continuous <- c(FALSE, FALSE, rep(TRUE, 2000L))
  density <- .brma_random_parameter_prior_density(
    samples      = samples,
    support      = c(0, 1),
    continuous   = continuous,
    point_masses = data.frame(x = c(0, 1), p = c(1 / 3, 1 / 3))
  )

  expect_s3_class(density, "prior_linear_density")
  expect_equal(density[["points"]][["p"]], c(1 / 3, 1 / 3))
  expect_equal(density[["density"]][["mass"]], 1 / 3)
})
