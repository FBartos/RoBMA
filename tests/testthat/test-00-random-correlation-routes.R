context("Random-slope correlation routes")

# A scaled us(1 + x | study) block: the original-scale intercept-slope
# correlation is derived from the LKJ correlation and the allocation-derived
# SDs. Small single-chain fit; the expectations are identities on its draws.
.random_correlation_cache <- new.env(parent = emptyenv())

# 'heterogeneity' is NULL for the default allocation of the heterogeneity
# over the block's SDs, or a prior_heterogeneity() specification.
.random_correlation_fit <- function(heterogeneity = NULL) {

  key <- if (is.null(heterogeneity)) "fit" else "fit_direct"
  if (!is.null(.random_correlation_cache[[key]])) {
    return(.random_correlation_cache[[key]])
  }

  set.seed(2)
  k   <- 30L
  dat <- data.frame(
    study = rep(sprintf("s%02d", seq_len(10L)), each = 3L),
    x     = stats::rnorm(k)
  )
  dat[["yi"]] <- 0.2 + 0.1 * dat[["x"]] +
    stats::rnorm(10L, 0, 0.1)[as.integer(factor(dat[["study"]]))] +
    stats::rnorm(k, 0, 0.15)
  args <- list(
    yi = quote(yi), V = diag(rep(0.0225, k)), random = ~ us(1 + x | study),
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    chains = 1, sample = 300, burnin = 100, adapt = 100, seed = 1,
    silent = TRUE
  )
  args[["prior_heterogeneity"]] <- heterogeneity
  fit <- suppressWarnings(do.call(brma.mv, args))
  .random_correlation_cache[[key]] <- fit

  return(fit)
}

# The same scaled block under BMA.mv, whose default heterogeneity null gates
# the block: a gate-only allocation, split over the block's SDs by an
# sd_component child, multiplies every SD of the block. The data have no
# study effects, so some draws switch the gate off.
.random_correlation_gated_fit <- function() {

  if (!is.null(.random_correlation_cache[["gated"]])) {
    return(.random_correlation_cache[["gated"]])
  }

  set.seed(2)
  k   <- 30L
  dat <- data.frame(
    study = rep(sprintf("s%02d", seq_len(10L)), each = 3L),
    x     = stats::rnorm(k)
  )
  dat[["yi"]] <- 0.2 + 0.1 * dat[["x"]] + stats::rnorm(k, 0, 0.15)
  fit <- suppressWarnings(BMA.mv(
    yi = yi, V = diag(rep(0.0225, k)), random = ~ us(1 + x | study),
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    chains = 1, sample = 300, burnin = 100, adapt = 100, seed = 1,
    silent = TRUE
  ))
  .random_correlation_cache[["gated"]] <- fit

  return(fit)
}

# Draws 3 and 4 allocate the whole block variance to the intercept: the SD
# of the scaled slope is zero and the original-scale correlation undefined.
.random_correlation_zero_sd_fit <- function() {

  fit   <- .random_correlation_fit()
  chain <- fit[["fit"]][["mcmc"]][[1L]]
  chain[3:4, "mu__xRE_ALLOCx_heterogeneity__weight[1]"] <- 1
  chain[3:4, "mu__xRE_ALLOCx_heterogeneity__weight[2]"] <- 0
  chain[3:4, "mu__xREx__study_intercept"] <-
    chain[3:4, "mu__xRE_ALLOCx_heterogeneity__allocation_sd"]
  chain[3:4, "mu__xREx__study_x"] <- 0
  fit[["fit"]][["mcmc"]][[1L]] <- chain

  return(fit)
}

.random_correlation_selection <- function(fit, name) {

  BayesTools::parameter_catalog_resolve(
    BayesTools::parameter_catalog(fit[["fit"]]),
    alias     = name,
    namespace = "mu"
  )
}

.random_correlation_odds <- function(p) p / (1 - p)


test_that("correlation hypotheses use the correlation draws", {

  skip_on_cran()
  fit       <- .random_correlation_fit()
  selection <- .random_correlation_selection(fit, "(mu) cor(intercept,x)")
  posterior <- as.numeric(as.matrix(BayesTools::parameter_draws(fit[["fit"]], selection)))

  # Prior correlation draws: parameter_draws() on the prior draws of the
  # same seed, with the allocation-derived SDs (the fitted-scale SD
  # quantities evaluated from the allocation sources alone) supplied as the
  # declared sources. The prior draws carry the deterministic SD monitors
  # the correlation depends on, equal to these SDs.
  sd_names <- c("mu__xREx__study_intercept", "mu__xREx__study_x")
  raw <- BayesTools::transform_prior_samples(
    fit[["fit"]], n_samples = 10000, seed = 1, formula_scale = list()
  )
  expect_true(all(sd_names %in% colnames(raw)))
  sources   <- raw[, setdiff(colnames(raw), sd_names), drop = FALSE]
  prior_fit <- BayesTools::JAGS_with_draws(fit[["fit"]], coda::mcmc.list(coda::mcmc(sources)))
  attr(prior_fit, "formula_scale") <- list()
  sds <- vapply(c(intercept = "(mu) sd(intercept)", x = "(mu) sd(x)"), function(name) {
    as.numeric(as.matrix(BayesTools::parameter_draws(
      prior_fit,
      .random_correlation_selection(fit, name)
    )))
  }, numeric(nrow(raw)))
  colnames(sds) <- sd_names
  expect_identical(unname(raw[, sd_names]), unname(sds))
  prior <- as.numeric(as.matrix(BayesTools::parameter_draws(
    fit[["fit"]], selection, model_samples = cbind(sources, sds)
  )))

  regions <- list(
    "rho(intercept,x) > 0" = function(x) x > 0,
    "rho(intercept,x) > -0.5 & rho(intercept,x) < 0.5" = function(x) x > -0.5 & x < 0.5
  )
  for (hypothesis in names(regions)) {
    region <- regions[[hypothesis]]
    result <- suppressWarnings(hypothesis(fit, hypothesis, columns = "all", seed = 1))
    # Region tests report odds.
    expect_equal(
      result[["posterior"]], .random_correlation_odds(mean(region(posterior))),
      tolerance = 1e-12, info = hypothesis
    )
    expect_equal(
      result[["prior"]], .random_correlation_odds(mean(region(prior))),
      tolerance = 1e-12, info = hypothesis
    )
  }

  # On the fitted scale the correlation is the LKJ(1) correlation, uniform on
  # (-1, 1): both regions have prior probability 1/2 (Monte Carlo SE 0.005
  # with 10,000 prior draws; tolerance 4 SE).
  for (hypothesis in names(regions)) {
    result <- suppressWarnings(hypothesis(
      fit, hypothesis, columns = "all", seed = 1,
      standardized_coefficients = TRUE
    ))
    prior_probability <- result[["prior"]] / (1 + result[["prior"]])
    expect_lt(abs(prior_probability - 0.5), 0.02)
  }
})


test_that("correlation hypotheses use the defined draws of a declared quantity", {

  skip_on_cran()
  fit <- .random_correlation_zero_sd_fit()

  # The random-parameter bundle keeps the undefined-draw declaration.
  selected <- .brma_random_parameter_select(fit, "rho(intercept,x)")
  expect_identical(which(is.na(selected[["samples"]][, 1L])), 3:4)
  expect_identical(
    BayesTools::posterior_metadata(selected[["samples"]], "undefined_draws"),
    stats::setNames("correlation", colnames(selected[["samples"]]))
  )

  selection <- .random_correlation_selection(fit, "(mu) cor(intercept,x)")
  posterior <- as.numeric(as.matrix(BayesTools::parameter_draws(fit[["fit"]], selection)))
  defined   <- posterior[!is.na(posterior)]
  expect_length(defined, 298L)

  result <- suppressWarnings(hypothesis(
    fit, "rho(intercept,x) > 0", columns = "all", seed = 1
  ))
  expect_equal(
    result[["posterior"]], .random_correlation_odds(mean(defined > 0)),
    tolerance = 1e-12
  )
  expect_true(paste0(
    "rho(intercept,x): computed from 298 of 300 posterior draws where the ",
    "correlation is defined, i.e. both SDs are positive."
  ) %in% attr(result, "footnotes"))

  # Plots of the correlation use the same defined draws.
  plotted <- .brma_random_parameter_mixed_posterior(fit, "rho(intercept,x)")
  expect_identical(as.numeric(plotted[[1L]]), defined)
})


test_that("summary footnotes follow the renamed random-effect rows", {

  skip_on_cran()
  fit       <- .random_correlation_zero_sd_fit()
  estimates <- summary(fit)[["estimates_random"]]
  footnotes <- attr(estimates, "footnotes")

  expect_true("rho(intercept,x)" %in% rownames(estimates))
  expect_false("cor(intercept,x)" %in% rownames(estimates))
  expect_identical(
    footnotes[["rho(intercept,x)"]],
    paste0(
      "rho(intercept,x): summarized over 298 of 300 draws where the ",
      "correlation is defined, i.e. both SDs are positive."
    )
  )
  expect_false(any(grepl("cor(intercept,x)", footnotes, fixed = TRUE)))
})


test_that("correlation plots draw the exact LKJ prior on the fitted scale", {

  skip_on_cran()
  fit <- .random_correlation_fit()

  # The 2 x 2 block has an LKJ(1) prior: the fitted-scale correlation r has
  # (r + 1) / 2 ~ Beta(1, 1), the density 1/2 on (-1, 1).
  posterior_only <- plot(fit, parameter = "rho(intercept,x)",
                         standardized_coefficients = TRUE, plot_type = "ggplot")
  with_prior <- plot(fit, parameter = "rho(intercept,x)", prior = TRUE,
                     standardized_coefficients = TRUE, plot_type = "ggplot")
  posterior_layers <- ggplot2::ggplot_build(posterior_only)[["data"]]
  layers <- ggplot2::ggplot_build(with_prior)[["data"]]
  expect_length(layers, length(posterior_layers) + 1L)
  prior_layer <- layers[[1L]]
  interior    <- prior_layer[["x"]] > -1 & prior_layer[["x"]] < 1
  expect_gt(sum(interior), 100L)
  expect_equal(
    prior_layer[["y"]][interior],
    stats::dbeta((prior_layer[["x"]][interior] + 1) / 2, 1, 1) / 2,
    tolerance = 1e-10
  )
})


test_that("correlation plots without an exact prior density draw the posterior alone", {

  skip_on_cran()
  # Half-normal priors on the block's SDs: the original-scale correlation of
  # the scaled block mixes the LKJ correlation with the SDs and has no exact
  # prior density.
  fit <- .random_correlation_fit(BayesTools::prior_random(
    sd = BayesTools::prior("normal", list(0, .5), list(0, Inf))
  ))
  posterior_only <- plot(fit, parameter = "rho(intercept,x)", plot_type = "ggplot")
  expect_warning(
    with_prior <- plot(fit, parameter = "rho(intercept,x)", prior = TRUE,
                       plot_type = "ggplot"),
    class = "BayesTools_prior_curve_unavailable"
  )
  expect_identical(
    ggplot2::ggplot_build(with_prior)[["data"]],
    ggplot2::ggplot_build(posterior_only)[["data"]]
  )
})


test_that("qCMDE/IWMDE stops for the correlation name the requested operation", {

  skip_on_cran()
  fit <- .random_correlation_fit()

  # The scaled block's original-scale correlation has no scalar source
  # coordinate for a qCMDE/IWMDE ordinate; the refusal names the operation.
  expect_identical(
    .brma_random_parameter_density_target(
      fit, "rho(intercept,x)", operation = "point hypotheses"
    )[["reason"]],
    paste0(
      "qCMDE/IWMDE point hypotheses are not available for random-effect ",
      "quantity 'rho(intercept,x)' because it has no supported scalar ",
      "random-component coordinate. Use density_method = 'KDE'."
    )
  )
  # Point hypotheses on it are refused for every method by its plan: the
  # correlation is atom-free (see the next test), but the original-scale
  # correlation of a scaled block mixes the LKJ correlation with the block's
  # SDs and has no prior density, so the Savage-Dickey ratio is undefined.
  for (method in c("KDE", "qCMDE", "IWMDE")) {
    expect_error(
      suppressWarnings(hypothesis(
        fit, "rho(intercept,x) = 0", density_method = method, seed = 1
      )),
      class = "RoBMA_hypothesis_target",
      info  = method
    )
  }
  expect_error(
    plot(fit, parameter = "rho(intercept,x)", density_method = "qCMDE"),
    paste0(
      "qCMDE/IWMDE plots are not available for random-effect quantity ",
      "'rho(intercept,x)'"
    ),
    fixed = TRUE
  )
})


test_that("original-scale correlations of allocated blocks are atom-free and plot", {

  skip_on_cran()
  # BayesTools declares the original-scale correlation of an allocated us()
  # block atom-free when no SD of the block has a point mass other than a
  # gate that scales all of them: the ungated default allocation of brma.mv
  # and the gated default allocation of BMA.mv. A draw with the gate off has
  # all SDs 0 and an undefined correlation; the defined draws are
  # continuous. The correlation has no prior density on the original scale.
  fits <- list(
    ungated = .random_correlation_fit(),
    gated   = .random_correlation_gated_fit()
  )
  for (name in names(fits)) {
    fit     <- fits[[name]]
    samples <- .brma_random_parameter_mixed_posterior(
      fit, "rho(intercept,x)"
    )[[1L]]
    expect_true(BayesTools::posterior_atoms_free(samples), info = name)
    expect_null(
      BayesTools::posterior_metadata(samples, "prior_density"),
      info = name
    )

    # The defined draws are the draws with the gate on.
    draws <- as.matrix(coda::as.mcmc(fit[["fit"]]))
    gate  <- grep("include_component_1_indicator$", colnames(draws), value = TRUE)
    expect_length(gate, if (identical(name, "gated")) 1L else 0L)
    n_defined <- if (length(gate) == 0L) nrow(draws) else sum(draws[, gate] == 1)
    expect_length(samples, n_defined)
    expect_false(anyNA(samples), info = name)
    if (identical(name, "gated")) {
      expect_lt(n_defined, nrow(draws))
    }

    # Plots draw the posterior; without a prior density the prior curve is
    # left out with a classed warning.
    posterior_only <- plot(fit, parameter = "rho(intercept,x)", plot_type = "ggplot")
    expect_warning(
      with_prior <- plot(fit, parameter = "rho(intercept,x)", prior = TRUE,
                         plot_type = "ggplot"),
      class = "BayesTools_prior_curve_unavailable"
    )
    expect_identical(
      ggplot2::ggplot_build(with_prior)[["data"]],
      ggplot2::ggplot_build(posterior_only)[["data"]],
      info = name
    )

    # Region hypotheses use the defined draws; point hypotheses are refused
    # for every method, as hypothesis_quantities() renders.
    region <- suppressWarnings(hypothesis(
      fit, "rho(intercept,x) > 0", columns = "all", seed = 1
    ))
    expect_equal(
      region[["posterior"]],
      .random_correlation_odds(mean(as.numeric(samples) > 0)),
      tolerance = 1e-12, info = name
    )
    for (method in c("KDE", "qCMDE", "IWMDE")) {
      expect_error(
        suppressWarnings(hypothesis(
          fit, "rho(intercept,x) = 0", density_method = method, seed = 1
        )),
        class = "RoBMA_hypothesis_target",
        info  = paste(name, method)
      )
    }
    quantities <- hypothesis_quantities(fit)
    row <- quantities[quantities[["alias"]] == "rho(intercept,x)", , drop = FALSE]
    expect_identical(nrow(row), 1L)
    expect_false(row[["point_test"]], info = name)
    expect_true(row[["direction_test"]], info = name)

    # On the fitted scale the correlation is the LKJ(1) correlation of the
    # 2 x 2 block, uniform on (-1, 1): its exact prior ordinate is 1/2.
    point <- suppressWarnings(hypothesis(
      fit, "rho(intercept,x) = 0", density_method = "KDE", columns = "all",
      seed = 1, standardized_coefficients = TRUE
    ))
    expect_equal(point[["prior"]], 0.5, tolerance = 1e-10, info = name)
  }
})
