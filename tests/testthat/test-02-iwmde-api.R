# ============================================================================ #
# test-02-iwmde-api.R
# ============================================================================ #

context("IWMDE public API")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


test_that("density methods are case-insensitive", {

  expect_equal(.density_method_normalize("kde"), "KDE")
  expect_equal(.density_method_normalize("QcMdE"), "qCMDE")
  expect_equal(.density_method_normalize("iwmde"), "IWMDE")
  expect_equal(.density_method_normalize("NORMAL", allow_normal = TRUE), "normal")
  expect_equal(.density_method_normalize_precomputed("qcmde"), "qCMDE")
  expect_equal(.density_method_normalize_precomputed("Iwmde"), "IWMDE")
})


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
      as_data     = TRUE
    )

    .expect_iwmde_ok(out, cases[[name]])
  }
})

test_that("IWMDE diagnostics can include requested density support values", {

  .skip_if_missing_raw_fits("bcg_meta-analysis")

  fit     <- load_fit("bcg_meta-analysis", validate = FALSE)
  context <- .iwmde_context(fit)
  base <- .iwmde_parameter_diagnostic(
    context              = context,
    parameter            = "mu",
    n_points             = 20,
    max_samples          = 20,
    normalization_points = 20,
    normalization_prob   = .999,
    diagnostic_cache     = .iwmde_diagnostic_cache()
  )
  include_value <- base[["xlim"]][2] + .05 * diff(base[["xlim"]])
  out <- .iwmde_parameter_diagnostic(
    context              = context,
    parameter            = "mu",
    n_points             = 20,
    max_samples          = 20,
    normalization_points = 20,
    normalization_prob   = .999,
    include_values       = include_value,
    diagnostic_cache     = .iwmde_diagnostic_cache()
  )

  expect_equal(out[["status"]], "ok")
  expect_true(include_value >= min(out[["iwmde"]][["x"]]))
  expect_true(include_value <= max(out[["iwmde"]][["x"]]))
  expect_true(any(abs(out[["iwmde"]][["x"]] - include_value) <=
                    sqrt(.Machine$double.eps) * max(1, abs(include_value))))
  expect_true(out[["diagnostics"]][["bf_included"]])
  expect_equal(out[["diagnostics"]][["bf_value"]], include_value)
  expect_true(is.finite(out[["diagnostics"]][["bf_ordinate"]]))
  expect_true(is.finite(out[["diagnostics"]][["bf_error_percent"]]))
})

test_that("IWMDE BF diagnostics retain boundary evaluation point", {

  evaluation_value <- 1e-6
  density <- list(
    x                = 0,
    evaluation_x     = evaluation_value,
    y                = 2,
    mcse             = .1,
    relative_mcse    = .05,
    finite_terms     = 40L,
    ess              = 25,
    max_weight_share = .2,
    max_log_ratio    = 1
  )

  bf_diagnostics <- .iwmde_density_bf_diagnostics(density, 0)
  diagnostic <- list(
    status      = "ok",
    diagnostics = c(
      bf_diagnostics,
      list(
        estimator              = "q_grid_cmde",
        weight_method          = "conditional_grid",
        active_mass            = 1,
        normalization_integral = 1,
        normalization_range    = c(evaluation_value / 2, evaluation_value * 2)
      )
    )
  )
  ordinate <- .iwmde_posterior_ordinate_attribute(diagnostic, "qCMDE")

  expect_true(bf_diagnostics[["bf_included"]])
  expect_equal(bf_diagnostics[["bf_value"]], 0)
  expect_equal(bf_diagnostics[["bf_evaluation_value"]], evaluation_value)
  expect_equal(ordinate[["value"]], 0)
  expect_equal(ordinate[["evaluation_value"]], evaluation_value)
  expect_true(.iwmde_posterior_ordinate_supports_bf(ordinate))
})

test_that("IWMDE diagnostics warn on unused dots", {

  .skip_if_missing_raw_fits("bcg_meta-analysis")

  expect_warning(
    plot_iwmde_diagnostics(
      object      = load_fit("bcg_meta-analysis", validate = FALSE),
      parameters  = "mu",
      n_points    = 20,
      max_samples = 20,
      method      = "iwmde",
      plot        = FALSE
    ),
    "Unused argument.*method"
  )
})

test_that("old IWMDE variant labels are rejected", {

  .skip_if_missing_raw_fits(c("bcg_meta-analysis", "bcg_meta-regression2"))

  expect_error(
    plot_iwmde_diagnostics(
      object         = load_fit("bcg_meta-analysis", validate = FALSE),
      parameters     = "mu",
      n_points       = 20,
      max_samples    = 20,
      density_method = "IWMDE-Chen",
      plot           = FALSE
    ),
    "must be one of"
  )
  expect_error(
    plot(
      load_fit("bcg_meta-analysis", validate = FALSE),
      parameter      = "mu",
      density_method = "IWMDE-CMDE"
    ),
    "must be one of"
  )
  expect_error(
    marginal_means(
      load_fit("bcg_meta-regression2", validate = FALSE),
      n_samples      = 1000,
      density_method = "IWMDE-Chen"
    ),
    "must be one of"
  )
})

test_that("unused dots warn in cleaned density interfaces", {

  .skip_if_missing_raw_fits(c("bcg_meta-analysis", "bcg_meta-regression2"))

  expect_warning(
    plot(
      load_fit("bcg_meta-analysis", validate = FALSE),
      parameter       = "mu",
      density_method  = "qCMDE",
      density_control = list(n_points = 20, max_samples = 20),
      iwmde_n_points  = 20
    ),
    "Unused argument.*iwmde_n_points"
  )
  expect_warning(
    marginal_means(
      load_fit("bcg_meta-regression2", validate = FALSE),
      n_samples          = 1000,
      density_method     = "qCMDE",
      density_control    = list(n_points = 20, max_samples = 20),
      iwmde_max_samples  = 20
    ),
    "Unused argument.*iwmde_max_samples"
  )
})

test_that("density_control validates public density settings", {

  valid <- .density_control_normalize(
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 20,
      max_samples          = 30,
      normalization_points = 40,
      normalization_prob   = .95,
      display_grid         = "uniform"
    )
  )

  expect_equal(valid[["n_points"]], 20)
  expect_equal(valid[["max_samples"]], 30)
  expect_equal(valid[["normalization_points"]], 40)
  expect_equal(valid[["normalization_prob"]], .95)
  expect_equal(valid[["display_grid"]], "uniform")
  expect_error(
    .density_control_normalize("qCMDE", list(unknown = 1)),
    "unrecognized"
  )
  expect_error(
    .density_control_normalize("qCMDE", list(20)),
    "fully named"
  )
  expect_error(
    .density_control_normalize(
      "qCMDE",
      structure(list(20, 30), names = c("n_points", "n_points"))
    ),
    "duplicate"
  )
  expect_error(
    .density_control_normalize("KDE", list(n_points = 20)),
    "only used"
  )
  expect_error(
    .density_control_normalize("IWMDE", list(normalization_points = 20)),
    "only available"
  )
  expect_error(
    .density_control_normalize("qCMDE", list(normalization_prob = 0)),
    "higher than 0"
  )
  expect_error(
    .density_control_normalize("qCMDE", list(normalization_prob = NA_real_)),
    "cannot contain"
  )
})

test_that("IWMDE row thinning is deterministic and equally spaced", {

  rows <- seq(2L, 80L, by = 2L)
  set.seed(123)
  old_seed <- .Random.seed

  out <- .iwmde_select_active_rows(rows = rows, max_samples = 10)

  expect_equal(out, rows[unique(round(seq(1, length(rows), length.out = 10)))])
  expect_identical(.Random.seed, old_seed)
  expect_equal(.iwmde_select_active_rows(rows = rows, max_samples = 10), out)
})

test_that("IWMDE condition rows keep one-row shape and ignore missing values", {

  context <- list(
    posterior_samples = data.frame(mu = 0),
    flat_prior_list   = list()
  )
  rows <- .iwmde_parameter_condition_rows(
    context = context,
    parameter_spec = list(
      conditional      = c("missing_a", NA, "missing_b"),
      conditional_rule = "OR"
    )
  )

  expect_length(rows, 1L)
  expect_false(rows)
})

test_that("IWMDE marginal-mean specs preserve child condition metadata", {

  condition_key <- paste(c("OR", 2L, "mu_alloc", "mu_intercept"), collapse = "\r")
  condition_event <- structure(
    list(
      conditional      = c("mu_intercept", "mu_alloc"),
      conditional_rule = "OR",
      condition_key    = condition_key
    ),
    class = "BayesTools_condition_event"
  )
  prior_density_context <- structure(
    list(condition_key = condition_key),
    class = "prior_density_context"
  )
  samples <- structure(
    stats::rnorm(20),
    linear_weights              = c(mu_intercept = 1, mu_alloc = 1),
    conditional                 = "stale_parent_condition",
    conditional_rule            = "AND",
    effective_conditional       = c("mu_intercept", "mu_alloc"),
    effective_conditional_rule  = "OR",
    condition_key               = condition_key,
    condition_event             = condition_event,
    resolved_condition_event    = condition_event,
    prior_density_context       = prior_density_context
  )
  marginal_means_object <- list(
    inference        = list(
      conditional = list(mu_alloc = list(alternate = samples))
    ),
    parameters       = "mu_alloc",
    term_map         = data.frame(
      term             = "alloc",
      parameter        = "mu_alloc",
      label            = "alloc",
      stringsAsFactors = FALSE
    ),
    conditional_rule = "AND"
  )

  specs <- .iwmde_marginal_means_specs(
    marginal_means_object = marginal_means_object,
    parameter             = "alloc",
    type                  = "conditional",
    levels                = "alternate"
  )
  spec <- specs[["alloc: alternate"]]

  expect_equal(spec[["conditional"]], c("mu_intercept", "mu_alloc"))
  expect_equal(spec[["conditional_rule"]], "OR")
  expect_equal(spec[["condition_key"]], condition_key)
  expect_identical(spec[["condition_event"]], condition_event)
  expect_identical(spec[["resolved_condition_event"]], condition_event)
  expect_identical(spec[["prior_density_context"]], prior_density_context)
  expect_equal(.iwmde_parameter_condition_key(spec), condition_key)
  expect_equal(
    .iwmde_posterior_metadata(samples, "mu_alloc", "alternate")[["condition_key"]],
    condition_key
  )

  stale_samples <- samples
  attr(stale_samples, "effective_conditional") <- NULL
  attr(stale_samples, "effective_conditional_rule") <- NULL
  attr(stale_samples, "condition_key") <- NULL
  marginal_means_object[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]] <- stale_samples

  expect_error(
    .iwmde_marginal_means_specs(
      marginal_means_object = marginal_means_object,
      parameter             = "alloc",
      type                  = "conditional",
      levels                = "alternate"
    ),
    "metadata is incomplete"
  )
})

test_that("IWMDE diagnostics compute without q-grid normalization", {

  .skip_if_missing_raw_fits("bcg_meta-analysis")

  out <- plot_iwmde_diagnostics(
    object         = load_fit("bcg_meta-analysis", validate = FALSE),
    parameters     = c("mu", "tau"),
    n_points       = 20,
    max_samples    = 50,
    density_method = "IWMDE",
    plot           = FALSE,
    as_data        = TRUE
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

test_that("IWMDE uses moment-matched recommended weights", {

  .skip_if_missing_raw_fits(c("bcg_meta-analysis", "konstantopoulos2011_3lvl"))

  meta_out <- plot_iwmde_diagnostics(
    object         = load_fit("bcg_meta-analysis", validate = FALSE),
    parameters     = c("mu", "tau"),
    n_points       = 20,
    max_samples    = 50,
    density_method = "IWMDE",
    plot           = FALSE,
    as_data        = TRUE
  )
  rho_out <- plot_iwmde_diagnostics(
    object         = load_fit("konstantopoulos2011_3lvl", validate = FALSE),
    parameters     = "rho",
    n_points       = 20,
    max_samples    = 30,
    density_method = "IWMDE",
    plot           = FALSE,
    as_data        = TRUE
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
    as_data               = TRUE
  )

  .expect_iwmde_ok(out, names(out))
  expect_named(out, c("alloc: alternate", "alloc: random", "alloc: systematic"))
  expect_s3_class(mm[["source_object"]], "brma")
})

test_that("IWMDE diagnostics compute for conditional estimated marginal means", {

  .skip_if_missing_raw_fits("dat.lehmann2018_RoBMA_mods")

  fit <- load_fit("dat.lehmann2018_RoBMA_mods", validate = FALSE)
  mm  <- marginal_means(fit, n_samples = 1000)
  out <- plot_iwmde_marginal_means_diagnostics(
    object                = fit,
    parameter             = "Preregistered",
    type                  = "conditional",
    levels                = "Not Pre-Registered",
    marginal_means_object = mm,
    n_points              = 20,
    max_samples           = 20,
    plot                  = FALSE,
    as_data               = TRUE
  )

  .expect_iwmde_ok(out, names(out))
  expect_named(out, "Preregistered: Not Pre-Registered")
  expect_true(out[["Preregistered: Not Pre-Registered"]][["diagnostics"]][["n_total"]] <
                length(mm[["inference"]][["averaged"]][["mu_Preregistered"]][["Not Pre-Registered"]]))
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
    as_data               = TRUE
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

test_that("IWMDE marginal means use attached source objects", {

  .skip_if_missing_raw_fits(c("bcg_meta-regression2", "bcg_meta-regression2b"))

  fit    <- load_fit("bcg_meta-regression2", validate = FALSE)
  fit_b  <- load_fit("bcg_meta-regression2b", validate = FALSE)
  mm     <- marginal_means(fit, n_samples = 1000)
  mm_old <- mm
  mm_old[["source_object"]] <- NULL

  out <- plot_iwmde_marginal_means_diagnostics(
    object                = fit_b,
    parameter             = "alloc",
    marginal_means_object = mm,
    n_points              = 20,
    max_samples           = 20,
    plot                  = FALSE,
    as_data               = TRUE
  )
  .expect_iwmde_ok(out, names(out))

  expect_error(
    plot_iwmde_marginal_means_diagnostics(
      object                = fit,
      marginal_means_object = mm_old,
      plot                  = FALSE,
      as_data               = TRUE
    ),
    "source fitted brma object"
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
    as_data     = TRUE
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
    max_samples = 20
  ))
})
