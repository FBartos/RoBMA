make_random_metadata_prior <- function(
    summary_type, label, parameter = "mu", block = NA_character_,
    grouping = NA_character_, allocation = NA_character_,
    component = NA_character_) {

  prior <- BayesTools::prior(
    distribution = "normal",
    parameters   = list(mean = 0, sd = 1)
  )
  attr(prior, "random_summary")         <- summary_type
  attr(prior, "random_summary_label")   <- label
  attr(prior, "parameter")              <- parameter
  attr(prior, "random_block")           <- block
  attr(prior, "random_grouping_factor") <- grouping
  attr(prior, "random_allocation")      <- allocation
  attr(prior, "random_component")       <- component
  prior
}


test_that("random-parameter sources follow canonical metadata", {

  raw_total <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
  attr(raw_total, "random_sd_total")  <- TRUE
  attr(raw_total, "parameter")        <- "mu"
  attr(raw_total, "random_allocation") <- "total"

  raw_prior <- function() {
    BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
  }
  prior_list <- list(
    total_source = raw_total,
    weight       = raw_prior(),
    sd_a         = raw_prior(),
    sd_b         = raw_prior(),
    rho          = raw_prior()
  )

  summaries <- list(
    sd_total = make_random_metadata_prior(
      "sd_total", "sd_total", allocation = "total"
    ),
    variance = make_random_metadata_prior(
      "var_frac", "var_frac", allocation = "total", component = "a"
    ),
    sd_a = make_random_metadata_prior(
      "sd", "sd(a)", block = "study", grouping = "study", component = "a"
    ),
    sd_b = make_random_metadata_prior(
      "sd", "sd(b)", block = "study", grouping = "study", component = "b"
    ),
    rho = make_random_metadata_prior(
      "rho", "rho(study)", block = "study", grouping = "study"
    ),
    correlation = make_random_metadata_prior(
      "cor", "cor(a, b)", block = "study", grouping = "study"
    )
  )
  attr(summaries[["variance"]], "random_allocation_metadata") <- list(
    weight_name = "weight"
  )

  formula_design <- list(mu = list(
    parameter      = "mu",
    random_effects = list(list(
      block_name         = "study",
      group_label        = "study",
      sd_parameter_names = c("sd_a", "sd_b"),
      sd_leaves          = list(leaf_terms = c(sd_a = "a", sd_b = "b")),
      correlation        = list(
        type        = "rho",
        rho_scale   = "rho",
        rho_name    = "rho",
        sample_name = "rho"
      )
    ))
  ))
  fit <- structure(
    list(),
    prior_list     = prior_list,
    formula_design = formula_design
  )

  names(summaries) <- paste0("(mu) ", names(summaries))
  specs   <- .brma_random_parameter_specs(summaries)
  sources <- .brma_random_parameter_source_names(
    fit            = fit,
    specs          = specs,
    summary_priors = summaries
  )

  expect_equal(
    sources,
    c("total_source", "weight", "sd_a", "sd_b", "rho", NA_character_)
  )
})


test_that("random-parameter source mapping rejects ambiguous metadata", {

  prior <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
  summary_prior <- make_random_metadata_prior(
    "sd", "sd(a)", block = "study", grouping = "study", component = "a"
  )
  summaries        <- list("(mu) sd(a)" = summary_prior)
  formula_design   <- list(mu = list(
    parameter      = "mu",
    random_effects = list(
      list(
        block_name         = "study",
        group_label        = "study",
        sd_parameter_names = c("sd_a", "sd_b")
      ),
      list(
        block_name         = "study",
        group_label        = "study",
        sd_parameter_names = "sd_c"
      )
    )
  ))
  fit <- structure(
    list(),
    prior_list     = list(sd_a = prior, sd_b = prior, sd_c = prior),
    formula_design = formula_design
  )

  expect_identical(
    .brma_random_parameter_source_names(
      fit            = fit,
      specs          = .brma_random_parameter_specs(summaries),
      summary_priors = summaries
    ),
    NA_character_
  )
})


test_that("random-parameter metadata preserves derived and semantic support", {

  prior <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
  rho_summary <- make_random_metadata_prior(
    "rho", "rho(study)", block = "study", grouping = "study"
  )
  summaries      <- list("(mu) rho(study)" = rho_summary)
  formula_design <- list(mu = list(
    parameter      = "mu",
    random_effects = list(list(
      block_name  = "study",
      group_label = "study",
      correlation = list(
        type         = "rho",
        rho_scale    = "rho",
        rho_name     = "rho",
        sample_name  = "rho",
        sample_fixed = 0
      )
    ))
  ))
  fit <- structure(
    list(),
    prior_list     = list(rho = prior),
    formula_design = formula_design
  )

  expect_identical(
    .brma_random_parameter_source_names(
      fit            = fit,
      specs          = .brma_random_parameter_specs(summaries),
      summary_priors = summaries
    ),
    NA_character_
  )

  simplex <- BayesTools::prior(
    "dirichlet",
    parameters = list(alpha = c(1, 1))
  )
  expect_equal(
    .brma_random_parameter_support(list(summary_type = "var_frac"),
                                   source_prior = simplex),
    c(0, 1)
  )

  multiplier <- make_random_metadata_prior(
    "sd_multiplier", "sd_mult(allocation: study)"
  )
  attr(multiplier, "random_allocation_metadata") <- list(
    scale     = "mean_variance",
    n_targets = 4L
  )
  expect_equal(
    .brma_random_parameter_support(
      list(summary_type = "sd_multiplier"),
      prior = multiplier
    ),
    c(0, 2)
  )
  attr(multiplier, "random_allocation_metadata") <- list(
    scale     = "total_variance",
    n_targets = 4L
  )
  expect_equal(
    .brma_random_parameter_support(
      list(summary_type = "sd_multiplier"),
      prior = multiplier
    ),
    c(0, 1)
  )
})
