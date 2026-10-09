# Attach a synthetic posterior to a constructor-only object: compile its
# formulas, then finalize the fit through the BayesTools fit contract.
.replay_object <- function(object, samples) {

  args <- .create_jags_formula_args(object[["data"]], object[["priors"]])
  compiled <- lapply(names(args[["formula_list"]]), function(parameter) {
    BayesTools::JAGS_formula(
      formula = args[["formula_list"]][[parameter]], parameter = parameter,
      data = args[["formula_data_list"]][[parameter]],
      prior_list = args[["formula_prior_list"]][[parameter]],
      formula_scale = args[["formula_scale_list"]][[parameter]],
      prior_random = args[["formula_random_prior_list"]][[parameter]],
      random_effects_compile = args[["formula_random_effects_compile_list"]][[parameter]]
    )
  })
  names(compiled) <- names(args[["formula_list"]])
  design <- lapply(compiled, `[[`, "formula_design")
  priors <- list()
  for (item in compiled) priors <- c(priors, item[["prior_list"]])

  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "list")
  )
  attr(fit, "formula_design") <- design
  attr(fit, "formula_scale")  <- lapply(compiled, `[[`, "formula_scale")
  attr(fit, "prior_list")     <- c(priors, .create_fit_priors(object[["data"]],
                                                                object[["priors"]]))
  object[["formula_design"]] <- design
  object[["fit"]]            <- as_bayestools_fit(fit)
  object
}

test_that("known-V estimate residuals remove the declared PET offset", {

  dat <- data.frame(yi = c(.1, .5, .3, .8), sei = c(.1, .3, .2, .4), id = factor(1:4))
  object <- bPET.mv(
    yi = yi, V = sei^2, data = dat, random = ~ 1 | id, measure = "GEN",
    prior_unit_information_sd = 1, prior_bias = BayesTools::prior_PET("normal", list(0, 1)),
    only_priors = TRUE, silent = TRUE
  )
  samples <- rbind(
    c(mu_intercept = .2, mu__xREx__id_intercept = .5, PET = 1.0),
    c(mu_intercept = .1, mu__xREx__id_intercept = .3, PET = -.5)
  )
  object <- .replay_object(object, samples)

  # The posterior mean of the latent effect given y = theta + e. With the bias
  # kept (bias_adjusted = FALSE) the PET term belongs to the location; when it
  # is removed, it is the declared offset of the observation
  # y = mu + u + PET * sei + e.
  blup <- function(bias_adjusted) {
    colMeans(t(vapply(1:2, function(draw) {
      bias   <- samples[draw, "PET"] * dat$sei
      mu     <- samples[draw, "mu_intercept"] + if (bias_adjusted) 0 else bias
      offset <- if (bias_adjusted) bias else 0
      tau2   <- samples[draw, "mu__xREx__id_intercept"]^2
      mu + tau2 / (tau2 + dat$sei^2) * (dat$yi - offset - mu)
    }, numeric(4L))))
  }
  for (bias_adjusted in c(FALSE, TRUE)) {
    residuals_estimate <- residuals(
      object, type = "outcome", conditioning_depth = "estimate",
      bias_adjusted = bias_adjusted
    )
    expect_equal(unname(as.numeric(residuals_estimate)),
                 dat$yi - blup(bias_adjusted), tolerance = 1e-12)
    predicted <- colMeans(as.matrix(predict(
      object, type = "blup", bias_adjusted = bias_adjusted, quiet = TRUE
    )))
    expect_equal(unname(as.numeric(residuals_estimate)), unname(dat$yi - predicted),
                 tolerance = 1e-12)
  }
})
