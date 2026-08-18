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
        label            = "total: var_prop(study)"
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
    label             = "total: var_prop(study)",
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
        alias      = "total: var_prop(study)",
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
