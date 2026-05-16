# ============================================================================ #
# test-02-iwmde-api.R
# ============================================================================ #

context("IWMDE public API")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


test_that("IWMDE diagnostics compute for representative cached fits", {

  fit_names <- c(
    "bcg_meta-analysis",
    "bcg_meta-regression",
    "bangertdrowns2004_location-scale",
    "konstantopoulos2011_3lvl",
    "dat.lehmann2018-PET",
    "dat.lehmann2018-PEESE",
    "dat.lehmann2018-3PSM",
    "bcg_glmm",
    "nielweise2008_glmm",
    "bcg_BMA.glmm"
  )
  .skip_if_missing_raw_fits(fit_names)

  cases <- list(
    "bcg_meta-analysis"                 = c("mu", "tau"),
    "bcg_meta-regression"               = c("mu_intercept", "mu_ablat", "tau"),
    "bangertdrowns2004_location-scale"  = c("mu_intercept", "log_tau_intercept"),
    "konstantopoulos2011_3lvl"          = c("mu", "tau", "rho", "gamma[1]"),
    "dat.lehmann2018-PET"               = c("mu", "tau", "PET"),
    "dat.lehmann2018-PEESE"             = c("mu", "tau", "PEESE"),
    "dat.lehmann2018-3PSM"              = c("mu", "tau"),
    "bcg_glmm"                          = c("mu", "tau", "pi[1]", "theta[1]"),
    "nielweise2008_glmm"                = c("mu", "tau", "phi[1]", "theta[1]"),
    "bcg_BMA.glmm"                      = c("mu", "tau", "pi[1]", "theta[1]")
  )

  for (name in names(cases)) {
    out <- plot_iwmde_diagnostics(
      object      = load_fit(name, validate = FALSE),
      parameters  = cases[[name]],
      n_points    = 20,
      max_samples = 20,
      plot        = FALSE,
      as_data     = TRUE,
      seed        = 1
    )

    .expect_iwmde_ok(out, cases[[name]])
  }
})

test_that("plain Chen IWMDE diagnostics compute without q-grid normalization", {

  .skip_if_missing_raw_fits("bcg_meta-analysis")

  out <- plot_iwmde_diagnostics(
    object      = load_fit("bcg_meta-analysis", validate = FALSE),
    parameters  = c("mu", "tau"),
    n_points    = 20,
    max_samples = 50,
    method      = "iwmde",
    plot        = FALSE,
    as_data     = TRUE,
    seed        = 1
  )

  expect_named(out, c("mu", "tau"))
  for (parameter in names(out)) {
    diagnostic <- out[[parameter]]
    expect_equal(diagnostic[["status"]], "ok")
    expect_equal(diagnostic[["diagnostics"]][["estimator"]], "iwmde")
    expect_true(startsWith(diagnostic[["diagnostics"]][["weight_method"]], "chen_"))
    expect_equal(diagnostic[["diagnostics"]][["normalization_points"]], 0L)
    expect_true(is.na(diagnostic[["diagnostics"]][["normalization_integral"]]))
    expect_true(all(is.finite(diagnostic[["iwmde"]][["y"]])))
    expect_true(all(diagnostic[["iwmde"]][["y"]] >= 0))
    expect_gt(diagnostic[["diagnostics"]][["integral"]], 0)
  }
})

test_that("plain Chen IWMDE uses moment-matched recommended weights", {

  .skip_if_missing_raw_fits(c("bcg_meta-analysis", "konstantopoulos2011_3lvl"))

  meta_out <- plot_iwmde_diagnostics(
    object      = load_fit("bcg_meta-analysis", validate = FALSE),
    parameters  = c("mu", "tau"),
    n_points    = 20,
    max_samples = 50,
    method      = "iwmde",
    plot        = FALSE,
    as_data     = TRUE,
    seed        = 1
  )
  rho_out <- plot_iwmde_diagnostics(
    object      = load_fit("konstantopoulos2011_3lvl", validate = FALSE),
    parameters  = "rho",
    n_points    = 20,
    max_samples = 30,
    method      = "iwmde",
    plot        = FALSE,
    as_data     = TRUE,
    seed        = 1
  )

  expect_equal(meta_out[["mu"]][["diagnostics"]][["weight_method"]],
               "chen_conditional_normal")
  expect_equal(meta_out[["tau"]][["diagnostics"]][["weight_method"]],
               "chen_exponential")
  expect_equal(rho_out[["rho"]][["diagnostics"]][["weight_method"]],
               "chen_power")
})

test_that("IWMDE diagnostics compute for estimated marginal means", {

  .skip_if_missing_raw_fits("bcg_meta-regression2")

  fit <- load_fit("bcg_meta-regression2", validate = FALSE)
  mm  <- marginal_means(fit, n_samples = 1000)
  out <- plot_iwmde_marginal_means_diagnostics(
    object                = fit,
    parameter             = "alloc",
    marginal_means_object = mm,
    n_points              = 20,
    max_samples           = 20,
    plot                  = FALSE,
    as_data               = TRUE,
    seed                  = 1
  )

  .expect_iwmde_ok(out, names(out))
  expect_named(out, c("alloc: alternate", "alloc: random", "alloc: systematic"))
  expect_false(is.null(mm[["iwmde_signature"]]))
})

test_that("IWMDE reuses identical marginal-mean linear targets", {

  .skip_if_missing_raw_fits("bcg_meta-regression4")

  fit <- load_fit("bcg_meta-regression4", validate = FALSE)
  mm  <- marginal_means(fit, n_samples = 1000)
  out <- plot_iwmde_marginal_means_diagnostics(
    object                = fit,
    parameter             = c("intercept", "alloc", "year_before1969", "alloc:year_before1969"),
    levels                = c("intercept", "alternate", "FALSE", "alternate, FALSE"),
    marginal_means_object = mm,
    n_points              = 20,
    max_samples           = 20,
    plot                  = FALSE,
    as_data               = TRUE,
    seed                  = 1
  )

  duplicate_targets <- c(
    "intercept: intercept",
    "alloc: alternate",
    "year_before1969: FALSE",
    "alloc:year_before1969: alternate, FALSE"
  )
  .expect_iwmde_ok(out[duplicate_targets], duplicate_targets)

  baseline <- out[[duplicate_targets[1L]]]
  for (target in duplicate_targets[-1L]) {
    expect_equal(out[[target]][["target_key"]], baseline[["target_key"]])
    expect_equal(out[[target]][["active_rows"]], baseline[["active_rows"]])
    expect_equal(out[[target]][["iwmde"]][["x"]], baseline[["iwmde"]][["x"]])
    expect_equal(out[[target]][["iwmde"]][["y"]], baseline[["iwmde"]][["y"]])
    expect_equal(out[[target]][["parameter"]], target)
  }
})

test_that("IWMDE marginal means reject incompatible precomputed objects", {

  .skip_if_missing_raw_fits(c("bcg_meta-regression2", "bcg_meta-regression2b"))

  fit    <- load_fit("bcg_meta-regression2", validate = FALSE)
  fit_b  <- load_fit("bcg_meta-regression2b", validate = FALSE)
  mm     <- marginal_means(fit, n_samples = 1000)
  mm_old <- mm
  mm_old[["iwmde_signature"]] <- NULL

  expect_error(
    plot_iwmde_marginal_means_diagnostics(
      object                = fit_b,
      marginal_means_object = mm,
      plot                  = FALSE,
      as_data               = TRUE
    ),
    "was not computed from"
  )
  expect_error(
    plot_iwmde_marginal_means_diagnostics(
      object                = fit,
      marginal_means_object = mm_old,
      plot                  = FALSE,
      as_data               = TRUE
    ),
    "does not contain an IWMDE fit signature"
  )
  expect_error(
    plot_iwmde_marginal_means_diagnostics(
      object                = fit,
      marginal_means_object = marginal_means(fit, n_samples = 1000, transform = "EXP"),
      plot                  = FALSE,
      as_data               = TRUE
    ),
    "linear predictor scale"
  )
})

test_that("IWMDE diagnostics handle RoBMA mixture branches", {

  .skip_if_missing_raw_fits("dat.lehmann2018_RoBMA")

  out <- plot_iwmde_diagnostics(
    object      = load_fit("dat.lehmann2018_RoBMA", validate = FALSE),
    parameters  = c("mu", "tau", "PET", "PEESE", "omega[1]", "eta[1]"),
    n_points    = 20,
    max_samples = 20,
    plot        = FALSE,
    as_data     = TRUE,
    seed        = 1
  )

  .expect_iwmde_ok(out[c("mu", "tau", "PET", "PEESE")], c("mu", "tau", "PET", "PEESE"))
  expect_equal(out[["omega[1]"]][["status"]], "unsupported")
  expect_equal(out[["eta[1]"]][["status"]], "unsupported")
  expect_gt(nrow(out[["PET"]][["point_masses"]]), 0)
  expect_gt(nrow(out[["PEESE"]][["point_masses"]]), 0)
  expect_equal(
    sum(out[["PET"]][["histogram"]][["density"]] * diff(out[["PET"]][["histogram"]][["breaks"]])),
    out[["PET"]][["active_mass"]],
    tolerance = 1e-8
  )
})

test_that("IWMDE diagnostics draw base plot matrix", {

  .skip_if_missing_raw_fits("dat.lehmann2018_RoBMA")

  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = file, width = 8, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot_iwmde_diagnostics(
    object      = load_fit("dat.lehmann2018_RoBMA", validate = FALSE),
    parameters  = c("mu", "tau", "PET", "PEESE"),
    n_points    = 20,
    max_samples = 20,
    seed        = 1
  ))
})
