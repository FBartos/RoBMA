.pr82_random_mixture_metadata <- function(external = FALSE, fixed_component = FALSE) {

  mixture <- BayesTools::prior_mixture(list(
    if (fixed_component) BayesTools::prior("point", list(.2)) else
      BayesTools::prior("normal", list(0, .25), truncation = list(0, Inf)),
    BayesTools::prior("normal", list(0, .75), truncation = list(0, Inf))
  ), is_null = c(FALSE, FALSE))
  random_prior <- if (external) {
    BayesTools::prior_random(study = BayesTools::random_block(
      sd_source = BayesTools::random_sd_source("tau")
    ))
  } else BayesTools::prior_random(sd = mixture)
  design <- BayesTools::JAGS_formula(
    ~ 1 + random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu", data = data.frame(study = factor(c("a", "a", "b", "b"))),
    prior_list = list(intercept = BayesTools::prior_mixture(list(
      BayesTools::prior("point", list(0)), BayesTools::prior("normal", list(0, 1))
    ), is_null = c(TRUE, FALSE))), prior_random = random_prior
  )
  source <- if (external) "tau" else design[["formula_design"]][["random_effects"]][[1L]][["sd_parameter_names"]]
  source <- unique(as.character(source))
  priors <- design[["prior_list"]]
  if (external) priors[[source]] <- mixture
  samples <- cbind(mu_intercept = c(.2, .3), mu_intercept_indicator = 2,
                   source = c(if (fixed_component) .2 else .4, .5),
                   source_indicator = c(1, 2))
  colnames(samples)[3:4] <- c(source, paste0(source, "_indicator"))
  fit <- coda::mcmc.list(coda::mcmc(samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "formula_design") <- list(mu = design[["formula_design"]])
  attr(fit, "prior_list") <- priors
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)
  list(context = list(
    data = structure(list(), random = TRUE), formula_fit = fit,
    indicator_names = c("mu_intercept_indicator", paste0(source, "_indicator")),
    flat_prior_list = priors, prior_cache = new.env(parent = emptyenv())
  ), samples = samples, source = source)
}

test_that("compiled random and external SD mixtures keep their original branch indices", {

  for (external in c(FALSE, TRUE)) for (fixed_component in c(FALSE, TRUE)) {
    fixture <- .pr82_random_mixture_metadata(external, fixed_component)
    preserved <- .iwmde_random_prior_indicators(fixture[["context"]])
    indicator <- paste0(fixture[["source"]], "_indicator")
    expect_identical(preserved, indicator)
    samples <- .iwmde_localize_active_branch_samples(
      fixture[["context"]], fixture[["samples"]],
      list(priors = list(location = list(intercept = BayesTools::prior("normal", list(0, 1)))),
           preserved_indicators = preserved)
    )
    expect_identical(samples[, indicator], c(1, 2))
    expect_identical(samples[, "mu_intercept_indicator"], c(1, 1))
    expect_identical(samples[, fixture[["source"]]], fixture[["samples"]][, fixture[["source"]]])
  }
})
