test_that("density warnings name their observed thresholds and sample control", {

  diagnostic <- list(estimator = "iwmde", relative_mcse = .05,
    sampling_relative_mcse = .05, ess = 50)
  expect_identical(.iwmde_diagnostics_bf_warning(diagnostic), c(
    paste0("IWMDE relative MCSE is ", .iwmde_percent(.05),
      " (warning threshold ", .iwmde_percent(.05), ")."),
    paste0("IWMDE finite-population row-sampling relative error is ", .iwmde_percent(.05),
      " (warning threshold ", .iwmde_percent(.05), ")."),
    paste0("IWMDE effective sample size is ", .iwmde_count(50),
      " (warning threshold ", .iwmde_count(100), ").")))
  expect_identical(.iwmde_fixed_sample_precision_action(list()), paste0(
    " Higher precision can be obtained by increasing 'samples' in the ",
    "'density_control' argument."))
  expect_identical(.iwmde_percent(Inf), "Inf%")
  expect_identical(.iwmde_percent(-Inf), "-Inf%")
  expect_identical(.iwmde_count(Inf), "Inf")
  expect_identical(.iwmde_count(-Inf), "-Inf")
  expect_identical(.iwmde_percent(NA_real_), "NA")
  expect_identical(.iwmde_count(NA_real_), "NA")
})


test_that("ordinate failure messages carry one inspect sentence", {

  # gregexpr() reports no match as the single position -1
  inspect_count <- function(text) sum(gregexpr("[Ii]nspect", text)[[1L]] > 0L)
  expect_identical(inspect_count("target is unsupported."), 0L)
  expect_identical(inspect_count("Inspect diagnostics."), 1L)
  expect_identical(inspect_count("Inspect diagnostics. Inspect again."), 2L)
  rejected <- .hypothesis_brma_iwmde_ordinate_failure_message("qCMDE", "mu = 0",
    list(status = "ok"), paste0("qCMDE lattice resolution error estimate is unbounded. ",
      "Increase 'normalization_points' in the 'density_control' argument and inspect ",
      "density_diagnostics()"))
  expect_identical(rejected, paste0(
    "qCMDE posterior ordinate for 'mu = 0' was rejected by diagnostics: ",
    "qCMDE lattice resolution error estimate is unbounded. Increase ",
    "'normalization_points' in the 'density_control' argument and inspect ",
    "density_diagnostics()."))
  underflow <- .hypothesis_brma_iwmde_ordinate_failure_message("qCMDE", "mu = 0",
    list(status = "ok"), .iwmde_ordinate_numerical_reason("underflow"),
    numerical_status = "underflow")
  expect_identical(inspect_count(underflow), 1L)
  expect_match(underflow, " Use a region hypothesis when a point ordinate is not representable.$")
  # A reason without its own hint still gets the remedy.
  plain <- .hypothesis_brma_iwmde_ordinate_failure_message("qCMDE", "mu = 0",
    list(status = "unsupported"), "target is unsupported")
  expect_identical(inspect_count(plain), 1L)
})


test_that("density messages print trimmed values and floor tiny percentages", {

  expect_identical(.iwmde_percent(.025), "2.5%")
  expect_identical(.iwmde_percent(.5), "50%")
  expect_identical(.iwmde_percent(0), "0%")
  expect_identical(.iwmde_percent(2e-6), "0.0002%")
  expect_identical(.iwmde_percent(6.7e-301), "<1e-4%")
  expect_identical(.iwmde_count(100), "100")
  expect_identical(.iwmde_count(67.8), "67.8")
  expect_identical(.iwmde_count(0), "0")
  expect_identical(.iwmde_count(5e-300), "5e-300")
  diagnostics <- list(estimator = "q_grid_cmde", ordinate_relative_change = .127,
    lattice_resolution_relative_error = 6.7e-301, truncation_ordinate_bound = .001,
    normalization_truncation_status = "bound")
  expect_identical(
    .iwmde_diagnostics_stability_description(diagnostics, "q_grid_cmde", .128),
    paste0(" posterior ordinate error estimate is 12.8% (grid change 12.7%; ",
           "lattice resolution <1e-4%; tail truncation 0.1%, upper bound)"))
})


test_that("density conditions retain their specific classes and common family", {

  ordinate <- .iwmde_ordinate_unavailable_condition("Unavailable", list())
  unavailable <- tryCatch(.density_diagnostics_stop_unavailable(), error = identity)
  plot <- .plot_brma_iwmde_unavailable_error(list(), "qCMDE", "mu")
  for (condition in list(ordinate, unavailable, plot)) {
    expect_s3_class(condition, "RoBMA_density_error")
    expect_s3_class(condition, "error")
  }
  expect_s3_class(ordinate, "RoBMA_density_ordinate_error")
  expect_s3_class(unavailable, "RoBMA_density_diagnostics_unavailable")
  expect_s3_class(plot, "RoBMA_density_plot_error")
  expect_identical(plot[["parameter"]], "mu")
  expect_identical(plot[["density_method"]], "qCMDE")
  expect_identical(plot[["density_diagnostics"]], list())
})


test_that("plot density remedies refer to their available diagnostic field", {

  diagnostics <- list(estimator = "q_grid_cmde", plot_scale_relative_mcse = .01,
    bulk_max_relative_mcse = .01, bulk_min_ess = 120, bulk_max_weight_share = .05,
    max_relative_mcse = .01, min_finite_terms = 600L, min_ess = 120,
    max_weight_share = .05, active_mass = 1, final_normalization_integral = 1,
    normalization_relative_error = 0, normalization_mass_ratio = 1,
    max_ordinate_relative_change = 0, max_lattice_resolution_relative_error = 0,
    max_normalizer_relative_change = 0, max_quadrature_relative_change = .06,
    n_estimator_rows = 600L, n_active_state_keys = 1L)
  samples <- structure(list(), iwmde_diagnostics = list(list(status = "ok",
    diagnostics = diagnostics)))
  condition <- .plot_brma_iwmde_unavailable_error(samples, "qCMDE", "mu")
  expect_match(conditionMessage(condition),
    "inspect the returned condition's 'density_diagnostics' field", fixed = TRUE)
  expect_false(grepl("density_diagnostics()", conditionMessage(condition), fixed = TRUE))
  expect_identical(condition[["density_diagnostics"]][[1L]][["diagnostics"]], diagnostics)
})
