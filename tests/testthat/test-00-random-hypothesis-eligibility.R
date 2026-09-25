.random_discovery_constant_draws_object <- function(prior) {

  result <- BayesTools::JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = data.frame(study = factor(c("a", "a", "b", "b"))),
    prior_list = list(intercept = BayesTools::prior("normal", list(0, 1))),
    prior_random = BayesTools::prior_random(sd = prior)
  )
  sd_name <- result[["formula_design"]][["random_effects"]][[1L]][["sd_parameter_names"]]
  samples <- cbind(mu_intercept = seq(-.2, .2, length.out = 20L), sd = .2)
  colnames(samples)[2L] <- sd_name
  fit <- coda::mcmc.list(coda::mcmc(samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = result[["formula_design"]])
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)

  structure(list(
    fit  = fit,
    data = structure(list(), random = TRUE)
  ), class = c("brma.mv", "brma"))
}


test_that("random hypothesis eligibility uses structural rather than sampled constancy", {

  for (fixed in c(FALSE, TRUE)) {
    prior <- if (fixed) {
      BayesTools::prior("spike", list(location = .2))
    } else {
      BayesTools::prior("gamma", list(shape = 2, rate = 2))
    }
    object <- .random_discovery_constant_draws_object(prior)
    quantities <- hypothesis_quantities(object)
    random <- quantities[quantities[["component"]] == "random", , drop = FALSE]

    expect_gt(nrow(random), 0L)
    expect_identical(unique(random[["direction_test"]]), !fixed)
    expect_identical(unique(random[["point_test"]]), !fixed)
    expect_identical(any(grepl("fixed by the fitted model", random[["reason"]])), fixed)
  }
})


test_that("sampled random quantities take region prior odds from the exact prior", {

  # The posterior draws of the study SD are constant (0.2), yet the quantity
  # is sampled: its region hypotheses use the prior odds of its canonical
  # prior density, Gamma(2, 2), not prior draws.
  object <- .random_discovery_constant_draws_object(
    BayesTools::prior("gamma", list(shape = 2, rate = 2))
  )
  quantities <- hypothesis_quantities(object)
  parameter  <- quantities[["parameter"]][quantities[["component"]] == "random"][[1L]]
  out <- suppressWarnings(.hypothesis_brma_random(
    object = object, parameter = parameter,
    hypothesis = BayesTools::hypothesis_rewrite(
      BayesTools::hypothesis_parse("theta > .5 vs theta <= .5"),
      c(theta = parameter)
    ),
    standardized_coefficients = FALSE, conditional = FALSE,
    logBF = FALSE, BF01 = FALSE, seed = 1, density_method = "KDE",
    n_samples = 20L, columns = "all"
  ))
  expect_equal(
    as.numeric(out[["prior"]]),
    stats::pgamma(.5, 2, 2, lower.tail = FALSE) / stats::pgamma(.5, 2, 2),
    tolerance = 1e-10
  )
})
