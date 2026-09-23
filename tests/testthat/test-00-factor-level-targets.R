context("Factor-level hypothesis targets and plots")

# Level labels are selectors, never JAGS coordinate positions. The labels
# 5/10/20 differ from every coordinate name, and the labels 1:4 coincide with
# the coordinate names of other levels (treatment level 3 is coordinate 2).
# Small single-chain fits suffice: every expectation below is an identity on
# the same posterior draws.
.factor_level_target_cache <- new.env(parent = emptyenv())

.factor_level_target_fits <- function() {

  if (!is.null(.factor_level_target_cache[["fits"]])) {
    return(.factor_level_target_cache[["fits"]])
  }

  set.seed(1)
  k    <- 48L
  data <- data.frame(
    g1  = factor(rep(c(5, 10, 20), length.out = k), levels = c(5, 10, 20)),
    g2  = factor(rep(1:4, length.out = k)[sample.int(k)], levels = 1:4),
    sei = stats::runif(k, 0.1, 0.3)
  )
  data[["yi"]] <- stats::rnorm(
    k,
    c(0, 0.2, 0.4)[as.integer(data[["g1"]])] +
      c(0, 0.1, 0.3, 0.5)[as.integer(data[["g2"]])],
    data[["sei"]]
  )
  fits <- lapply(c(treatment = "treatment", meandif = "meandif"), function(contrast) {
    suppressWarnings(brma(
      yi = yi, sei = sei, mods = ~ g1 + g2, data = data, measure = "SMD",
      set_contrast_factor_predictors = contrast,
      chains = 1, sample = 1000, burnin = 200, adapt = 100,
      seed = 1, silent = TRUE
    ))
  })
  .factor_level_target_cache[["fits"]] <- fits

  return(fits)
}


# The coefficient transformation target that hypothesis() pairs with a level,
# resolved through the same steps as hypothesis.brma().
.factor_level_prior_target <- function(fit, hypothesis) {

  metadata <- .brma_parameter_catalog_metadata(fit)
  parsed   <- BayesTools::hypothesis_parse(
    hypothesis,
    catalog        = metadata[["catalog"]],
    simplify_names = TRUE
  )
  selected <- .hypothesis_brma_select_parameter(
    object     = fit,
    hypothesis = parsed,
    component  = "auto",
    metadata   = metadata
  )
  rewritten <- .hypothesis_brma_rewrite(
    hypothesis = parsed,
    aliases    = selected[["aliases"]],
    parameter  = selected[["parameter"]]
  )
  point_refs <- .hypothesis_brma_point_refs(
    hypothesis     = rewritten,
    parameter      = selected[["parameter"]],
    require_direct = TRUE
  )
  targets <- .hypothesis_brma_formula_coefficient_level_targets(
    object     = fit,
    selected   = selected,
    point_refs = point_refs
  )

  return(vapply(targets, `[[`, character(1), "target"))
}


test_that("factor-level point hypotheses use the level's own fitted coordinate", {

  skip_on_cran()
  fit   <- .factor_level_target_fits()[["treatment"]]
  mcmc  <- as.matrix(fit[["fit"]][["mcmc"]])
  cases <- data.frame(
    hypothesis = c("g1[10] = 0", "g1[20] = 0", "g2[2] = 0", "g2[3] = 0", "g2[4] = 0"),
    term       = c("g1", "g1", "g2", "g2", "g2"),
    level      = c("10", "20", "2", "3", "4"),
    coordinate = c("mu_g1[1]", "mu_g1[2]", "mu_g2[1]", "mu_g2[2]", "mu_g2[3]"),
    stringsAsFactors = FALSE
  )
  compared <- character()

  for (i in seq_len(nrow(cases))) {
    # Savage-Dickey ratio of the treatment coordinate: its normal prior
    # density over the unbounded Gaussian KDE of its draws at the null. The
    # comparison needs a null inside the KDE grid (levels 20 and 4 lie far
    # from zero; their ordinates are extrapolated).
    prior     <- fit[["priors"]][["mods"]][[cases[["term"]][i]]]
    posterior <- stats::density(mcmc[, cases[["coordinate"]][i]])
    height    <- stats::approx(posterior[["x"]], posterior[["y"]], xout = 0)[["y"]]
    if (!is.na(height)) {
      direct_BF <- stats::dnorm(
        0,
        prior[["parameters"]][["mean"]],
        prior[["parameters"]][["sd"]]
      ) / height
      result <- suppressWarnings(hypothesis(
        fit, cases[["hypothesis"]][i], density_method = "KDE"
      ))
      expect_equal(
        attr(result, "raw_BF"), direct_BF,
        tolerance = 1e-10, info = cases[["hypothesis"]][i]
      )
      compared <- c(compared, cases[["hypothesis"]][i])
    }
    expect_identical(
      .factor_level_prior_target(fit, cases[["hypothesis"]][i]),
      stats::setNames(cases[["coordinate"]][i], cases[["level"]][i])
    )
  }
  expect_identical(compared, c("g1[10] = 0", "g2[2] = 0", "g2[3] = 0"))
})


test_that("mean-difference level point hypotheses stop without a single coordinate", {

  skip_on_cran()
  fit <- .factor_level_target_fits()[["meandif"]]

  for (hypothesis in c("g1[10] = 0", "g2[2] = 0")) {
    expect_error(
      suppressWarnings(hypothesis(fit, hypothesis, density_method = "KDE")),
      "is absent from the fitted coefficient transformation",
      fixed = TRUE,
      info = hypothesis
    )
  }
})


test_that("factor-cell plots draw the selected level for every contrast", {

  skip_on_cran()
  fits  <- .factor_level_target_fits()
  cases <- data.frame(
    parameter = c("g1[10]", "g1[20]", "g2[2]", "g2[3]", "g2[4]"),
    term      = c("g1", "g1", "g2", "g2", "g2"),
    level     = c(2L, 3L, 2L, 3L, 4L),
    stringsAsFactors = FALSE
  )

  for (contrast in names(fits)) {
    fit  <- fits[[contrast]]
    mcmc <- as.matrix(fit[["fit"]][["mcmc"]])
    for (i in seq_len(nrow(cases))) {
      term        <- cases[["term"]][i]
      n_levels    <- if (identical(term, "g1")) 3L else 4L
      coordinates <- paste0("mu_", term, "[", seq_len(n_levels - 1L), "]")
      # Level values from the fitted coordinates and the contrast matrix.
      contrasts <- if (identical(contrast, "treatment")) {
        stats::contr.treatment(n_levels)
      } else {
        BayesTools::contr.meandif(n_levels)
      }
      expected <- as.numeric(
        mcmc[, coordinates, drop = FALSE] %*% contrasts[cases[["level"]][i], ]
      )

      entry <- .brma_parameter_select_entry(
        fit, cases[["parameter"]][i], allow_factor_cells = TRUE
      )
      cell <- .plot_brma_factor_cell_samples(
        object                    = fit,
        entry                     = entry,
        standardized_coefficients = FALSE,
        conditional               = FALSE,
        precomputed               = FALSE
      )
      expect_equal(
        as.numeric(cell[["samples"]][[1L]]), expected,
        tolerance = 1e-10, info = paste(contrast, cases[["parameter"]][i])
      )
    }
  }

  plotted <- suppressWarnings(plot(
    fits[["meandif"]], parameter = "g1[10]", plot_type = "ggplot",
    density_method = "KDE"
  ))
  expect_s3_class(plotted, "ggplot")
})
