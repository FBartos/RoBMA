context("Random-effect quantities of standardized coefficients")

# The values of the predictor of the stand-in models below: with
# 'formula_scale', its coefficients are fitted on the standardized predictor
# (the sample SD of these values), so that the original-scale SD of a random
# slope is its fitted-scale SD times 1 / sd(x_values).
.scaled_slope_x_values <- function() rep(c(2, 4, 6, 9), each = 3L)

# A fitted-object stand-in with the random slope (0 + x || g) whose SD has a
# Uniform(1, 2) prior on the fitted scale; its draws are samples of that prior
# (the posterior is the prior, the Bayes factors are near 1). Without
# 'formula_scale' the predictor is not standardized and both scales agree.
.scaled_slope_object <- function(formula_scale = TRUE, n = 2000L, seed = 42L) {

  data <- data.frame(
    x = .scaled_slope_x_values(),
    g = factor(rep(1:4, each = 3L))
  )
  compiled <- BayesTools::JAGS_formula(
    ~ 1 + x + (0 + x || g), "mu", data,
    list(
      intercept = BayesTools::prior("normal", list(0, 1)),
      x         = BayesTools::prior("normal", list(0, 1))
    ),
    prior_random = BayesTools::prior_random(g = BayesTools::random_block(
      sd = BayesTools::prior("uniform", list(1, 2))
    )),
    formula_scale = formula_scale
  )
  sd_name <- unique(
    compiled[["formula_design"]][["random_effects"]][[1L]][["sd_parameter_names"]]
  )
  set.seed(seed)
  draws <- cbind(mu_intercept = stats::rnorm(n), mu_x = stats::rnorm(n))
  draws <- cbind(draws, stats::setNames(
    data.frame(stats::runif(n, 1, 2)), sd_name
  ))
  fit <- coda::mcmc.list(coda::mcmc(as.matrix(draws)))
  attr(fit, "prior_list")     <- compiled[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = compiled[["formula_design"]])
  attr(fit, "formula_scale")  <- list(mu = compiled[["formula_scale"]])

  structure(list(
    fit  = as_bayestools_fit(fit),
    data = structure(list(), random = TRUE)
  ), class = c("brma.mv", "brma"))
}

# The log Bayes factor of 'q = value' for the draws of a quantity with the prior
# object that declares its law and bounds: a computation of the draws that
# reads no support of the catalog.
.scaled_slope_reference <- function(draws, prior, value) {

  attr(BayesTools::hypothesis_BF(
    as.numeric(draws), prior, paste0("q = ", format(value, digits = 17)),
    parameter = "q", seed = 1
  ), "raw_log_BF")
}


test_that("the support and the Bayes factors of standardized random slopes are those of the fitted scale", {

  object <- .scaled_slope_object()
  original_to_fitted <- 1 / stats::sd(.scaled_slope_x_values())

  for (standardized in c(TRUE, FALSE)) {
    # The bounds of the Uniform(1, 2) prior of the SD on the scale of the draws.
    k    <- if (standardized) 1 else original_to_fitted
    info <- paste("standardized_coefficients =", standardized)

    for (quantity in c("sd", "var")) {
      selector <- paste0("(mu) ", quantity, "(x)")
      power    <- if (identical(quantity, "sd")) 1 else 2
      bounds   <- k^power * c(1, 2^power)
      selected <- .brma_random_parameter_select(
        object, selector, standardized_coefficients = standardized
      )
      draws <- .brma_random_parameter_mixed_posterior(
        object, selector, standardized_coefficients = standardized,
        selected = selected
      )[[1L]]
      expect_equal(
        .brma_random_parameter_support(selected), bounds,
        info = paste(info, selector)
      )
      expect_equal(
        BayesTools::posterior_metadata(draws, "support")[["bounds"]], bounds,
        info = paste(info, selector)
      )
      expect_true(all(draws >= bounds[1L] & draws <= bounds[2L]),
                  info = paste(info, selector))
    }

    # The Savage-Dickey ratio reflects the posterior density at the bounds of
    # the prior; a variance is tested through its SD.
    sd_draws <- .brma_random_parameter_mixed_posterior(
      object, "(mu) sd(x)", standardized_coefficients = standardized
    )[[1L]]
    for (value in c(1.05, 1.5, 1.95)) {
      reference <- .scaled_slope_reference(
        sd_draws, BayesTools::prior("uniform", list(k, 2 * k)), k * value
      )
      for (quantity in c("sd", "var")) {
        point <- (k * value)^(if (identical(quantity, "sd")) 1 else 2)
        got <- hypothesis(
          object, paste0("`(mu) ", quantity, "(x)` = ", format(point, digits = 17)),
          standardized_coefficients = standardized, density_method = "KDE", seed = 1
        )
        expect_equal(
          attr(got, "raw_log_BF"), reference, tolerance = 1e-8,
          info = paste(info, quantity, value)
        )
      }
    }

    # Discovery tests the SD and the variance with a point hypothesis.
    quantities <- hypothesis_quantities(object, transform_scaled = !standardized)
    random     <- quantities[quantities[["parameter"]] %in%
                               c("(mu) tau(x)", "(mu) tau2(x)"), , drop = FALSE]
    expect_gt(nrow(random), 0L)
    expect_true(all(random[["point_test"]]), info = info)
  }
})


test_that("the scales of the random-effect quantities of unscaled predictors agree", {

  object <- .scaled_slope_object(formula_scale = FALSE)

  selected <- lapply(c(TRUE, FALSE), function(standardized) {
    .brma_random_parameter_select(
      object, "(mu) sd(x)", standardized_coefficients = standardized
    )
  })
  for (selection in selected) {
    expect_equal(.brma_random_parameter_support(selection), c(1, 2))
  }
  draws <- lapply(c(TRUE, FALSE), function(standardized) {
    .brma_random_parameter_mixed_posterior(
      object, "(mu) sd(x)", standardized_coefficients = standardized
    )[[1L]]
  })
  expect_identical(as.numeric(draws[[1L]]), as.numeric(draws[[2L]]))
  expect_identical(
    BayesTools::posterior_metadata(draws[[1L]], "support"),
    BayesTools::posterior_metadata(draws[[2L]], "support")
  )
  for (standardized in c(TRUE, FALSE)) {
    got <- hypothesis(
      object, "`(mu) sd(x)` = 1.5",
      standardized_coefficients = standardized, density_method = "KDE", seed = 1
    )
    expect_equal(
      attr(got, "raw_log_BF"),
      .scaled_slope_reference(draws[[1L]], BayesTools::prior("uniform", list(1, 2)), 1.5),
      tolerance = 1e-8
    )
  }
})


# A stand-in whose intercept SD is fixed (a spike) on the fitted scale and whose
# slope SD is sampled: the block has a scaled predictor, so the catalog of the
# fit (the original scale) takes both SDs for sampled, while the fitted scale
# fixes the intercept SD.
.fixed_intercept_sd_object <- function(n = 1000L, seed = 2L) {

  data <- data.frame(
    x = .scaled_slope_x_values(),
    g = factor(rep(1:4, each = 3L))
  )
  compiled <- BayesTools::JAGS_formula(
    ~ 1 + x + (1 + x || g), "mu", data,
    list(
      intercept = BayesTools::prior("normal", list(0, 1)),
      x         = BayesTools::prior("normal", list(0, 1))
    ),
    prior_random = BayesTools::prior_random(g = BayesTools::random_block(
      sd = BayesTools::prior("spike", list(1.5))
    )),
    formula_scale = TRUE
  )
  sd_names <- unique(
    compiled[["formula_design"]][["random_effects"]][[1L]][["sd_parameter_names"]]
  )
  compiled[["prior_list"]][[sd_names[2L]]] <- BayesTools::prior("uniform", list(1, 2))
  draws <- BayesTools:::.generate_prior_sample_matrix(
    compiled[["prior_list"]], n, seed = seed
  )
  fit <- coda::mcmc.list(coda::mcmc(draws))
  attr(fit, "prior_list")     <- compiled[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = compiled[["formula_design"]])
  attr(fit, "formula_scale")  <- list(mu = compiled[["formula_scale"]])

  structure(list(
    fit  = as_bayestools_fit(fit),
    data = structure(list(), random = TRUE)
  ), class = c("brma.mv", "brma"))
}


test_that("standardized random-effect quantities are classified on the fitted scale", {

  object <- .fixed_intercept_sd_object()

  # The intercept SD is fixed on the fitted scale and sampled on the original
  # scale (the intercept of a block with a scaled slope).
  standardized <- .brma_random_parameter_select(
    object, "(mu) sd(intercept)", standardized_coefficients = TRUE
  )
  original <- .brma_random_parameter_select(
    object, "(mu) sd(intercept)", standardized_coefficients = FALSE
  )
  expect_identical(standardized[["spec"]][["status"]], "structural")
  expect_identical(original[["spec"]][["status"]], "sampled")
  expect_identical(standardized[["entry"]][["selection"]][["quantity_id"]],
                   original[["entry"]][["selection"]][["quantity_id"]])
  expect_identical(
    standardized[["entry"]][["selection"]],
    BayesTools::parameter_catalog_resolve(
      BayesTools::parameter_catalog(object[["fit"]], transform_scaled = FALSE),
      "(mu) sd(intercept)"
    )
  )
  expect_identical(
    original[["entry"]][["selection"]],
    BayesTools::parameter_catalog_resolve(
      BayesTools::parameter_catalog(object[["fit"]]), "(mu) sd(intercept)"
    )
  )

  # A fixed quantity is refused as one, not as a quantity without a prior mass.
  expect_error(
    hypothesis(
      object, "`(mu) sd(intercept)` > 1", standardized_coefficients = TRUE,
      density_method = "KDE", seed = 1
    ),
    "fixed by the fitted model", class = "BayesTools_hypothesis_fixed"
  )

  # The slope SD is sampled on both scales, and the draws of a quantity do not
  # depend on the catalog it is read from.
  for (standardized in c(TRUE, FALSE)) {
    selected <- .brma_random_parameter_select(
      object, "(mu) sd(x)", standardized_coefficients = standardized
    )
    expect_identical(selected[["spec"]][["status"]], "sampled")
  }
})


# The nested allocation of a gamma SD over two blocks and its split into the
# components of a correlated block with a scaled slope (the correlation is
# monitored through its Cholesky factor only); the draws are samples of the
# priors.
.cholesky_only_allocation_object <- function(n = 1000L, seed = 42L) {

  data <- data.frame(
    study = factor(rep(1:4, each = 3L)),
    drug  = factor(rep(1:3, 4L)),
    x     = rep(c(2, 4, 6), 4L)
  )
  compiled <- BayesTools::JAGS_formula(
    ~ 1 + random(1 + x | study, name = "study", covariance = "us") +
      random(1 | drug, name = "drug", covariance = "diag"),
    "mu", data,
    list(intercept = BayesTools::prior("normal", list(0, 1))),
    formula_scale = TRUE,
    prior_random  = BayesTools::prior_random(
      root  = BayesTools::random_variance_allocation(
        name = "root", terms = c(study = "study", drug = "drug"),
        sd = BayesTools::prior("gamma", list(2, 2))
      ),
      split = BayesTools::random_variance_allocation(
        name = "split", terms = "study", target = "sd_component",
        parent = BayesTools::allocation_ref("root", "study")
      ),
      study = BayesTools::random_block(cor = BayesTools::prior_lkj(2))
    )
  )
  draws <- BayesTools:::.generate_prior_sample_matrix(
    compiled[["prior_list"]], n, seed = seed
  )
  term      <- compiled[["formula_design"]][["random_effects"]][[1L]]
  primitive <- term[["correlation"]][["primitive_names"]]
  draws     <- cbind(draws, matrix(
    stats::rbeta(nrow(draws), 2, 2), ncol = 1L,
    dimnames = list(NULL, primitive)
  ))
  r         <- 2 * draws[, primitive] - 1
  cholesky  <- cbind(1, r, 0, sqrt(1 - r^2))
  colnames(cholesky) <- as.vector(
    BayesTools:::.bt_random_effect_cholesky_names(term, 2L)
  )
  draws <- cbind(draws, cholesky[, !colnames(cholesky) %in% colnames(draws), drop = FALSE])
  # A valid Cholesky-monitored fit need not monitor the primitive variables.
  draws <- draws[, !colnames(draws) %in% primitive, drop = FALSE]

  chains  <- coda::mcmc.list(coda::mcmc(draws))
  attr(chains, "prior_list")     <- compiled[["prior_list"]]
  attr(chains, "formula_design") <- list(mu = compiled[["formula_design"]])
  attr(chains, "formula_scale")  <- list(mu = compiled[["formula_scale"]])
  carrier <- as_bayestools_fit(chains)
  fit <- structure(
    list(mcmc = chains, sample = nrow(draws), summary.pars = list(mutate = NULL),
         monitor = colnames(draws)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  for (name in c("prior_list", "formula_design", "formula_scale", "parameter_map")) {
    attr(fit, name) <- attr(carrier, name)
  }
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)

  structure(list(
    fit  = fit,
    data = structure(list(), random = TRUE)
  ), class = c("brma.mv", "brma"))
}


test_that("standardized random-effect hypotheses on prior draws extract the fitted-scale quantity", {

  object <- .cholesky_only_allocation_object()
  statement <- "`(mu) study: sd(x)` > .5"

  # The region mass of the quantity comes from its prior draws, extracted from
  # the prior samples of the fitted scale.
  got <- hypothesis(
    object, statement, standardized_coefficients = TRUE,
    density_method = "KDE", seed = 1, n_samples = 1000
  )
  # The same fitted-scale quantity of the same fit through BayesTools.
  expected <- BayesTools::hypothesis_BF(
    object[["fit"]], "(mu) study: sd(x) > .5", seed = 1, n_samples = 1000
  )
  expect_equal(attr(got, "raw_log_BF"), attr(expected, "raw_log_BF"), tolerance = 1e-8)

  selected <- .brma_random_parameter_select(
    object, "(mu) study: sd(x)", standardized_coefficients = TRUE,
    prior = TRUE, n_prior_samples = 1000, seed = 1
  )
  expect_identical(nrow(selected[["samples"]]), 1000L)
  expect_true(all(is.finite(selected[["samples"]])))
})
