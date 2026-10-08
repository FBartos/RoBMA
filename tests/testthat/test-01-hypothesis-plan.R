context("Cached hypothesis-plan fits")

source(testthat::test_path("common-functions.R"))

# Keep the complete builders here: this fitting source owns their cache hash.
test_that("hypothesis plan fits are prepared for cached post-fit tests", {

  skip_on_cran()
  skip_refit_if_cached("hypothesis-plan")
  keys <- hypothesis_plan_fit_names()
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") ||
    is_false_env("ROBMA_TEST_SKIP_REFIT")
  needs <- setNames(vapply(keys, function(key) {
    refit || !is_cached_fit_valid(key)
  }, logical(1)), names(keys))

  fit <- function(data, mods, prior_mods = NULL, contrast = "treatment",
                  scale = NULL) {
    arguments <- list(
      yi = quote(yi), sei = quote(sei), mods = mods, data = data,
      measure = "SMD", set_contrast_factor_predictors = contrast,
      chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
      silent = TRUE
    )
    if (!is.null(prior_mods)) {
      arguments[["prior_mods"]] <- prior_mods
    }
    if (!is.null(scale)) {
      arguments[["scale"]] <- scale
    }
    suppressWarnings(do.call(brma, arguments))
  }
  set.seed(1)
  k <- 48L
  factors <- data.frame(
    g1  = factor(rep(c(5, 10, 20), length.out = k), levels = c(5, 10, 20)),
    g2  = factor(rep(1:4, length.out = k)[sample.int(k)], levels = 1:4),
    sei = stats::runif(k, 0.1, 0.3)
  )
  factors[["yi"]] <- stats::rnorm(
    k,
    c(0, 0.2, 0.4)[as.integer(factors[["g1"]])] +
      c(0, 0.1, 0.3, 0.5)[as.integer(factors[["g2"]])],
    factors[["sei"]]
  )
  fits <- lapply(
    c(treatment = "treatment", meandif = "meandif",
      orthonormal = "orthonormal", independent = "independent"),
    function(contrast) {
      if (!needs[[contrast]]) {
        return(load_fit(keys[[contrast]]))
      }
      fit(
        data     = factors,
        mods     = if (identical(contrast, "independent")) ~ g1 + g2 - 1 else ~ g1 + g2,
        contrast = contrast
      )
    }
  )
  for (n_levels in 2:4) {
    set.seed(3)
    labels <- c("a", "b", "c", "d")[seq_len(n_levels)]
    ordered <- data.frame(
      g   = factor(rep(labels, length.out = k), levels = labels),
      sei = stats::runif(k, 0.1, 0.3)
    )
    ordered[["yi"]] <- stats::rnorm(
      k, seq(0, 0.3, length.out = n_levels)[as.integer(ordered[["g"]])],
      ordered[["sei"]]
    )
    name <- paste0("ordered_", n_levels)
    if (!needs[[name]]) {
      fits[[name]] <- load_fit(keys[[name]])
      next
    }
    fits[[name]] <- fit(
      data       = ordered,
      mods       = ~ g,
      prior_mods = list(g = BayesTools::prior_ordered(
        BayesTools::prior("normal", list(0, 1))
      ))
    )
  }
  set.seed(5)
  scaled <- data.frame(
    x   = stats::rnorm(k, 1, 2),
    g   = factor(rep(c("a", "b"), length.out = k)),
    sei = stats::runif(k, 0.1, 0.3)
  )
  scaled[["yi"]] <- stats::rnorm(k, 0.1 * scaled[["x"]], scaled[["sei"]])
  if (needs[["scaled_interaction"]]) {
    fits[["scaled_interaction"]] <- fit(data = scaled, mods = ~ x * g)
  } else {
    fits[["scaled_interaction"]] <- load_fit(keys[["scaled_interaction"]])
  }
  # A factor term in the location and the scale formula: its alias 'g1'
  # names a quantity of each component.
  if (needs[["location_scale"]]) {
    fits[["location_scale"]] <- fit(
      data = factors, mods = ~ g1, scale = ~ g1, contrast = "meandif"
    )
  } else {
    fits[["location_scale"]] <- load_fit(keys[["location_scale"]])
  }
  if (needs[["model_averaged"]]) {
    fits[["model_averaged"]] <- suppressWarnings(BMA.norm(
      yi = yi, sei = sei, mods = ~ g1, data = factors, measure = "SMD",
      chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
      silent = TRUE
    ))
  } else {
    fits[["model_averaged"]] <- load_fit(keys[["model_averaged"]])
  }
  # A RoBMA ensemble with publication-bias components; its bias rows are not
  # hypothesis targets.
  if (needs[["robma_mixture"]]) {
    fits[["robma_mixture"]] <- suppressWarnings(RoBMA(
      yi = yi, sei = sei, mods = ~ g1, data = factors, measure = "SMD",
      chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
      silent = TRUE
    ))
  } else {
    fits[["robma_mixture"]] <- load_fit(keys[["robma_mixture"]])
  }

  if (needs[["two_scale"]]) {
    data <- data.frame(
      yi     = c(0.08, 0.13, 0.18, 0.20, 0.01, 0.05),
      study  = rep(c("s1", "s2", "s3"), each = 2L),
      effect = rep(c("a", "b"), 3L),
      x      = c(0, 1, 0, 1, 0, 1),
      g      = factor(c("u", "v", "v", "u", "u", "v"))
    )
    V <- kronecker(diag(3L), matrix(c(0.04, 0.018, 0.018, 0.05), nrow = 2L))
    fits[["two_scale"]] <- suppressWarnings(brma.mv(
      yi = yi, V = V, data = data, measure = "GEN", mods = ~ g,
      random = list(study = ~ 1 | study, effect = ~ 1 | study:effect),
      scale  = list(study = ~ x, effect = ~ x),
      prior_unit_information_sd = 1,
      chains = 1, sample = 500, burnin = 200, adapt = 200, seed = 1,
      silent = TRUE,
      convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
    ))
  } else {
    fits[["two_scale"]] <- load_fit(keys[["two_scale"]])
  }

  if (needs[["meandif_mt"]]) {
    set.seed(1)
    k <- 48L
    data <- data.frame(
      g1  = factor(rep(c(5, 10, 20), length.out = k), levels = c(5, 10, 20)),
      sei = stats::runif(k, 0.1, 0.3)
    )
    data[["yi"]] <- stats::rnorm(
      k, c(0, 0.2, 0.4)[as.integer(data[["g1"]])], data[["sei"]]
    )
    fits[["meandif_mt"]] <- suppressWarnings(brma(
      yi = yi, sei = sei, mods = ~ g1, data = data, measure = "SMD",
      set_contrast_factor_predictors = "meandif",
      prior_mods = list(g1 = BayesTools::prior_factor(
        "mt", list(location = 0, scale = 0.5, df = 3), contrast = "meandif"
      )),
      chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
      silent = TRUE
    ))
  } else {
    fits[["meandif_mt"]] <- load_fit(keys[["meandif_mt"]])
  }

  if (needs[["treatment"]]) {
    save_fit("hypothesis_plan_treatment", fits[["treatment"]])
  }
  if (needs[["meandif"]]) {
    save_fit("hypothesis_plan_meandif", fits[["meandif"]])
  }
  if (needs[["orthonormal"]]) {
    save_fit("hypothesis_plan_orthonormal", fits[["orthonormal"]])
  }
  if (needs[["independent"]]) {
    save_fit("hypothesis_plan_independent", fits[["independent"]])
  }
  if (needs[["ordered_2"]]) {
    save_fit("hypothesis_plan_ordered_2", fits[["ordered_2"]])
  }
  if (needs[["ordered_3"]]) {
    save_fit("hypothesis_plan_ordered_3", fits[["ordered_3"]])
  }
  if (needs[["ordered_4"]]) {
    save_fit("hypothesis_plan_ordered_4", fits[["ordered_4"]])
  }
  if (needs[["scaled_interaction"]]) {
    save_fit("hypothesis_plan_scaled_interaction", fits[["scaled_interaction"]])
  }
  if (needs[["location_scale"]]) {
    save_fit("hypothesis_plan_location_scale", fits[["location_scale"]])
  }
  if (needs[["model_averaged"]]) {
    save_fit("hypothesis_plan_model_averaged", fits[["model_averaged"]])
  }
  if (needs[["robma_mixture"]]) {
    save_fit("hypothesis_plan_robma_mixture", fits[["robma_mixture"]])
  }
  if (needs[["two_scale"]]) {
    save_fit("hypothesis_plan_two_scale", fits[["two_scale"]])
  }
  if (needs[["meandif_mt"]]) {
    save_fit("hypothesis_plan_meandif_mt", fits[["meandif_mt"]])
  }

  for (key in keys) {
    expect_identical(validate_cached_fit(key, deep = TRUE), character())
  }
})
