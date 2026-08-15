context("zplot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))

# list cached fits lazily
skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

zplot_smoke_samples <- test_profile_value(100L, 1000L)

.test_as_zplot <- function(..., max_samples = 1000) {

  as_zplot(..., max_samples = max_samples)
}

.test_zplot <- function(..., summary_max_samples = 1000, max_samples = 1000) {

  zplot(..., summary_max_samples = summary_max_samples, max_samples = max_samples)
}

.test_plot_zplot <- function(..., max_samples = 1000) {

  plot(..., max_samples = max_samples)
}

.test_lines_zplot <- function(..., max_samples = 1000) {

  lines(..., max_samples = max_samples)
}

.zplot_test_area <- function(df) {

  x <- df[["x"]]
  y <- df[["y"]]

  return(sum(diff(x) * (y[-length(y)] + y[-1L]) / 2))
}

.zplot_test_edr <- function(predictive, significance_level) {

  S        <- nrow(predictive[["mu"]])
  K        <- ncol(predictive[["mu"]])
  sei      <- predictive[["sei"]]
  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- .root_sum_squares(predictive[["tau_within"]], sei_mat)

  return(rowMeans(
    stats::pnorm(
      matrix(significance_level * sei, nrow = S, ncol = K, byrow = TRUE),
      mean       = predictive[["mu"]],
      sd         = total_sd,
      lower.tail = FALSE
    ) +
      stats::pnorm(
        matrix(-significance_level * sei, nrow = S, ncol = K, byrow = TRUE),
        mean       = predictive[["mu"]],
        sd         = total_sd,
        lower.tail = TRUE
      )
  ))
}


.zplot_test_conditional_variance <- function(Q, V) {

  return(diag(Q - Q %*% solve(Q + V, Q)))
}


test_that("zplot Gaussian conditional variance matches covariance identities", {

  Q <- matrix(c(
    1.2, 0.4,
    0.4, 0.8
  ), nrow = 2, byrow = TRUE)
  V <- matrix(c(
    0.5, 0.15,
    0.15, 0.7
  ), nrow = 2, byrow = TRUE)

  expect_equal(
    .zplot_gaussian_conditional_variance(Q, V),
    .zplot_test_conditional_variance(Q, V),
    tolerance = 1e-12
  )

  singular_V <- matrix(0.3, nrow = 2, ncol = 2)
  expect_equal(
    .zplot_gaussian_conditional_variance(Q, singular_V),
    .zplot_test_conditional_variance(Q, singular_V),
    tolerance = 1e-12
  )

  tau_within  <- matrix(c(0.2, 0.4, 0.3, 0.5), nrow = 2, byrow = TRUE)
  tau_between <- matrix(c(0.6, 0.7, 0.4, 0.8), nrow = 2, byrow = TRUE)
  vi          <- c(0.1, 0.2)
  observed <- .zplot_multilevel_conditional_variance(
    tau_within  = tau_within,
    tau_between = tau_between,
    vi          = vi,
    cluster     = c(1L, 1L)
  )
  expected <- matrix(NA_real_, nrow = 2, ncol = 2)
  for (s in seq_len(nrow(expected))) {
    Q_s <- diag(tau_within[s, ]^2) + tcrossprod(tau_between[s, ])
    expected[s, ] <- .zplot_test_conditional_variance(Q_s, diag(vi))
  }

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("zplot rejects GLMM fits early", {

  name <- "bcg_BMA.glmm"
  skip_if_missing_fits(name)

  expect_error(
    .test_as_zplot(fits[[name]]),
    "normal outcome models"
  )
})


test_that("zplot creates reusable objects and plots directly", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit <- fits[[name]]
  zp  <- .test_as_zplot(fit, max_samples = zplot_smoke_samples)

  expect_s3_class(zp, "zplot_brma")
  expect_named(zp[["zplot"]], c("estimates", "data"))
  expect_identical(
    zp[["zplot"]][["data"]][["conditioning_depth"]],
    "marginal"
  )
  expect_true(.is_ggplot(.test_zplot(
    fit,
    plot_type           = "ggplot",
    summary_max_samples = zplot_smoke_samples,
    max_samples         = zplot_smoke_samples
  )))
  expect_true(.is_ggplot(.test_zplot(
    zp,
    plot_type   = "ggplot",
    max_samples = zplot_smoke_samples
  )))
})

test_that("zplot allows descriptive known-V brma.mv displays", {

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)

  zc <- .test_as_zplot(
    fits[[name]],
    max_samples = zplot_smoke_samples
  )

  expect_s3_class(zc, "zplot_brma")
  expect_true(all(is.finite(zc[["zplot"]][["estimates"]][["EDR"]])))
  expect_true(all(is.finite(zc[["zplot"]][["estimates"]][["weights"]])))
  expect_equal(
    length(zc[["zplot"]][["data"]][["z"]]),
    nobs(fits[[name]])
  )

  density <- .test_lines_zplot(
    zc,
    as_data     = TRUE,
    max_samples = zplot_smoke_samples,
    plot_ci     = FALSE,
    length.out  = 25
  )
  expect_true(all(is.finite(density[["y"]])))
})


test_that("zplot known-V brma.mv uses marginal predictive components by default", {

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)

  fit_brma          <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
  selected_rows     <- seq_len(min(25L, nrow(posterior_samples)))
  posterior_samples <- posterior_samples[selected_rows, , drop = FALSE]

  expected_mu <- predict.brma(
    object             = fit_brma,
    type               = "terms",
    bias_adjusted      = TRUE,
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  expected_components <- .brma_mv_heterogeneity_components(
    object            = fit_brma,
    posterior_samples = posterior_samples
  )
  components <- .zplot_predictive_components(
    object            = fit_brma,
    posterior_samples = posterior_samples,
    extrapolate       = TRUE
  )

  expect_equal(components[["mu"]], as.matrix(expected_mu))
  expect_equal(
    components[["tau_within"]],
    .total_brma_mv_heterogeneity_samples(expected_components)
  )
  expect_equal(components[["sei"]], .outcome_data_sei(fit_brma))
})


test_that("zplot conditioning depths select the intended replication target", {

  name <- "konstantopoulos2011_3lvl"
  skip_if_missing_fits(name)

  fit_brma          <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
  max_samples       <- min(25L, nrow(posterior_samples))
  selected_rows     <- .thin_sample_rows(nrow(posterior_samples), max_samples)
  if (!is.null(selected_rows)) {
    posterior_samples <- posterior_samples[selected_rows, , drop = FALSE]
  }

  marginal <- .zplot_predictive_components(
    object             = fit_brma,
    posterior_samples  = posterior_samples,
    extrapolate        = TRUE,
    conditioning_depth = "marginal"
  )
  cluster <- .zplot_predictive_components(
    object             = fit_brma,
    posterior_samples  = posterior_samples,
    extrapolate        = TRUE,
    conditioning_depth = "cluster"
  )
  estimate <- .zplot_predictive_components(
    object             = fit_brma,
    posterior_samples  = posterior_samples,
    extrapolate        = TRUE,
    conditioning_depth = "estimate"
  )

  tau_result <- .evaluate.brma.tau(
    fit               = fit_brma[["fit"]],
    scale_data        = fit_brma[["data"]][["scale"]],
    scale_formula     = NULL,
    scale_priors      = fit_brma[["priors"]][["scale"]],
    is_scale          = FALSE,
    is_multilevel     = TRUE,
    K                 = nobs(fit_brma),
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(fit_brma[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(fit_brma[["priors"]])
  )

  expect_equal(marginal[["tau_within"]], tau_result[["tau_total"]])
  expect_equal(cluster[["tau_within"]], tau_result[["tau_within"]])
  expect_equal(
    marginal[["tau_within"]]^2,
    cluster[["tau_within"]]^2 + tau_result[["tau_between"]]^2
  )
  expected_estimate_variance <- matrix(
    NA_real_,
    nrow = nrow(posterior_samples),
    ncol = nobs(fit_brma)
  )
  vi            <- .outcome_data_sei(fit_brma)^2
  block_indices <- .get_multilevel_block_indices(
    fit_brma[["data"]][["outcome"]][["cluster"]]
  )
  for (s in seq_len(nrow(posterior_samples))) {
    for (idx in block_indices) {
      Q_s <- diag(
        tau_result[["tau_within"]][s, idx]^2,
        nrow = length(idx),
        ncol = length(idx)
      ) +
        tcrossprod(tau_result[["tau_between"]][s, idx])
      expected_estimate_variance[s, idx] <-
        .zplot_test_conditional_variance(Q_s, diag(vi[idx]))
    }
  }
  expect_equal(
    estimate[["tau_within"]]^2,
    expected_estimate_variance,
    tolerance = 1e-12
  )
  expect_true(all(estimate[["tau_within"]] > 0))
  expect_true(all(estimate[["tau_within"]] <= marginal[["tau_within"]]))
  expect_equal(marginal[["mu"]], as.matrix(predict.brma(
    fit_brma,
    type               = "terms",
    bias_adjusted      = TRUE,
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )))
  expect_equal(cluster[["mu"]], as.matrix(predict.brma(
    fit_brma,
    type               = "cluster",
    bias_adjusted      = TRUE,
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )))
  expect_equal(estimate[["mu"]], as.matrix(predict.brma(
    fit_brma,
    type               = "estimate",
    bias_adjusted      = TRUE,
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )))

  significance_level <- stats::qnorm(0.975)
  zc_marginal         <- .test_as_zplot(
    fit_brma,
    significance_level = significance_level,
    max_samples        = max_samples
  )
  zc_cluster          <- .test_as_zplot(
    fit_brma,
    significance_level = significance_level,
    max_samples        = max_samples,
    conditioning_depth = "cluster"
  )
  zc_estimate         <- .test_as_zplot(
    fit_brma,
    significance_level = significance_level,
    max_samples        = max_samples,
    conditioning_depth = "estimate"
  )

  expect_equal(
    zc_marginal[["zplot"]][["estimates"]][["EDR"]],
    .zplot_test_edr(marginal, significance_level)
  )
  expect_equal(
    zc_cluster[["zplot"]][["estimates"]][["EDR"]],
    .zplot_test_edr(cluster, significance_level)
  )
  expect_equal(
    zc_estimate[["zplot"]][["estimates"]][["EDR"]],
    .zplot_test_edr(estimate, significance_level)
  )
  expect_identical(
    zc_marginal[["zplot"]][["data"]][["conditioning_depth"]],
    "marginal"
  )
  expect_identical(
    zc_cluster[["zplot"]][["data"]][["conditioning_depth"]],
    "cluster"
  )
  expect_identical(
    zc_estimate[["zplot"]][["data"]][["conditioning_depth"]],
    "estimate"
  )

  legacy_zplot <- zc_cluster
  legacy_zplot[["zplot"]][["data"]][["conditioning_depth"]] <- NULL
  expect_identical(.zplot_stored_conditioning_depth(legacy_zplot), "cluster")
})


test_that("zplot estimate depth integrates ordinary conditional uncertainty", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit_brma          <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
  posterior_samples <- posterior_samples[
    seq_len(min(25L, nrow(posterior_samples))),
    ,
    drop = FALSE
  ]
  estimate <- .zplot_predictive_components(
    object             = fit_brma,
    posterior_samples  = posterior_samples,
    extrapolate        = TRUE,
    conditioning_depth = "estimate"
  )
  tau_result <- .zplot_tau_samples(
    object            = fit_brma,
    posterior_samples = posterior_samples
  )
  tau2     <- tau_result[["tau_within"]]^2
  vi       <- .outcome_data_sei(fit_brma)^2
  expected <- tau2 * matrix(vi, nrow = nrow(tau2), ncol = ncol(tau2),
                            byrow = TRUE) /
    sweep(tau2, 2L, vi, "+")

  expect_equal(estimate[["tau_within"]]^2, expected, tolerance = 1e-12)
  expect_true(all(estimate[["tau_within"]] > 0))
})


test_that("zplot estimate depth integrates known-V and random-formula uncertainty", {

  fit_names <- c("brma.mv_block_mvn", "brma.mv_block_mvn_random")
  skip_if_missing_fits(fit_names)

  for (name in fit_names) {
    fit_brma          <- fits[[name]]
    posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
    posterior_samples <- posterior_samples[
      seq_len(min(5L, nrow(posterior_samples))),
      ,
      drop = FALSE
    ]
    estimate <- .zplot_predictive_components(
      object             = fit_brma,
      posterior_samples  = posterior_samples,
      extrapolate        = TRUE,
      conditioning_depth = "estimate"
    )
    V        <- .known_v_materialize(.data_known_v_data(fit_brma[["data"]]))
    K        <- nobs(fit_brma)
    expected <- matrix(NA_real_, nrow = nrow(posterior_samples), ncol = K)

    if (.is_random(fit_brma)) {
      random_vcov <- .brma_mv_random_effects_marginal_vcov(
        object            = fit_brma,
        posterior_samples = posterior_samples,
        diagonal_only     = FALSE,
        data              = fit_brma[["data"]],
        new_levels        = "error"
      )[["samples"]]
      for (s in seq_len(nrow(expected))) {
        expected[s, ] <- .zplot_test_conditional_variance(
          random_vcov[s, , ],
          V
        )
      }
    } else {
      tau_result <- .zplot_tau_samples(
        object            = fit_brma,
        posterior_samples = posterior_samples
      )
      for (s in seq_len(nrow(expected))) {
        Q_s <- diag(tau_result[["tau_within"]][s, ]^2, nrow = K, ncol = K)
        expected[s, ] <- .zplot_test_conditional_variance(Q_s, V)
      }
    }

    expect_equal(
      estimate[["tau_within"]]^2,
      expected,
      tolerance = 1e-10,
      info = name
    )
  }
})


test_that("zplot validates cluster conditioning", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  expect_error(
    .test_as_zplot(fits[[name]], conditioning_depth = "cluster"),
    "only available for multilevel models"
  )
})


# ============================================================================ #
# Test: Simple Meta-Analysis Zplot
# ============================================================================ #

test_that("zplot for simple meta-analysis renders base and ggplot output", {

  name     <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  zc       <- .test_as_zplot(fit_brma)

  # --------------------------------------------------
  # Visual tests
  # --------------------------------------------------

  expect_vdiffr_snapshot("zplot_simple_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base"))
  })

  expect_vdiffr_snapshot(
    "zplot_simple_ggplot",
    suppressMessages(.test_plot_zplot(zc, plot_type = "ggplot"))
  )
})


# ============================================================================ #
# Test: Zplot Customization
# ============================================================================ #

test_that("zplot customization snapshots are stable", {

  skip_if_not_full_visuals("Customization snapshots are visual-gallery coverage.")

  name     <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  zc       <- .test_as_zplot(fits[[name]])

  # --------------------------------------------------
  # Test custom styles
  # --------------------------------------------------

  expect_vdiffr_snapshot("zplot_custom_base", function() {
    suppressMessages(.test_plot_zplot(
      zc, plot_type = "base",
      plot_fit = TRUE, plot_ci = TRUE,
      lwd = 2, lty = 2,           # line args
      dots_hist = list(col = "lightblue"),
      main = "Custom Zplot"
    ))
  })

  # For ggplot, we pass specific list args for hist/lines as per .get_dots_* in implementation
  # or global args that get mapped
  expect_vdiffr_snapshot(
    "zplot_custom_ggplot",
    suppressMessages(.test_plot_zplot(
      zc, plot_type = "ggplot",
      dots_hist = list(fill = "lightblue", color = "blue"),
      dots_thresholds = list(color = "red", linetype = "dashed"),
      main = "Custom Zplot GGplot"
    ))

  )

  # --------------------------------------------------
  # Test components only (hist / lines)
  # --------------------------------------------------

  expect_vdiffr_snapshot("zplot_hist_only_base", function() {
    suppressMessages(hist(zc, plot_type = "base", main = "Hist Only"))
  })

  expect_vdiffr_snapshot("zplot_lines_only_base", function() {
    # lines() adds to existing plot usually, but here we test the function
    # so we create an empty plot and add lines
    plot(0, 0, type = "n", xlim = c(-6, 6), ylim = c(0, 0.5), main = "Lines Only")
    .test_lines_zplot(zc, plot_type = "base", col = "purple")
  })

})


# ============================================================================ #
# Test: Meta-Regression Zplot
# ============================================================================ #

test_that("zplot for meta-regression renders base output", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_regression_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "Meta-Regression Zplot"))
  })

})

# ============================================================================ #
# Test: Selection Models Zplot
# ============================================================================ #

test_that("zplot for positive-direction selection model renders base output", {

  name <- "dat.lehmann2018-3PSM"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_selection_pos_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "Selection Model (Pos) Zplot"))
  })

})

test_that("zplot for negative-direction selection model renders base output", {

  skip_if_not_full_visuals("Negative-direction selection zplot is gallery coverage.")

  name <- "dat.lehmann2018-3PSM_neg"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_selection_neg_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "Selection Model (Neg) Zplot"))
  })

})

test_that("zplot handles RoBMA bias-mixture branches", {

  name <- "dat.lehmann2018_RoBMA"
  skip_if_not(name %in% names(fits), "RoBMA cached fit not available.")

  max_samples        <- test_profile_value(100L, 1000L)
  integration_points <- test_profile_value(401L, 2001L)
  fit               <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit[["fit"]])
  if (nrow(posterior_samples) > max_samples) {
    selected_ind      <- round(seq(from = 1, to = nrow(posterior_samples), length.out = max_samples))
    posterior_samples <- posterior_samples[selected_ind, , drop = FALSE]
  }
  selection         <- .zplot_selection_context(
    object            = fit,
    posterior_samples = posterior_samples,
    is_weightfunction = .is_weightfunction(fit)
  )
  weighted_rows <- which(!selection[["use_normal"]])
  skip_if_not(length(weighted_rows) > 0, "No weightfunction posterior rows in cached fit.")

  selection_args <- .zplot_selection_args(
    selection = selection,
    row       = weighted_rows[1],
    estimate  = 1,
    n         = 3
  )
  active_cuts <- selection[["n_bins"]] - 1L

  expect_equal(nrow(selection_args[["omega"]]), 3)
  expect_equal(ncol(selection_args[["omega"]]), selection[["n_bins"]])
  expect_length(selection_args[["crit_yi"]], active_cuts)

  zc <- .test_as_zplot(fit, max_samples = max_samples)
  expect_true(all(is.finite(zc[["zplot"]][["estimates"]][["EDR"]])))
  expect_true(all(zc[["zplot"]][["estimates"]][["EDR"]] >= 0))
  expect_true(all(zc[["zplot"]][["estimates"]][["EDR"]] <= 1))
  expect_true(all(is.finite(zc[["zplot"]][["estimates"]][["weights"]])))

  fitted_density <- .test_lines_zplot(
    zc, as_data = TRUE, max_samples = max_samples, plot_ci = FALSE,
    extrapolate = FALSE, length.out = 25
  )
  extrapolated_density <- .test_lines_zplot(
    zc, as_data = TRUE, max_samples = max_samples, plot_ci = FALSE,
    extrapolate = TRUE, length.out = 25
  )

  expect_true(all(is.finite(unlist(fitted_density[c("y", "y_lCI", "y_uCI")]))))
  expect_true(all(is.finite(unlist(extrapolated_density[c("y", "y_lCI", "y_uCI")]))))

  fitted_area       <- .zplot_test_area(.test_lines_zplot(
    zc,
    as_data     = TRUE,
    max_samples = max_samples,
    plot_ci     = FALSE,
    extrapolate = FALSE,
    from        = -20,
    to          = 20,
    length.out  = integration_points
  ))
  extrapolated_area <- .zplot_test_area(.test_lines_zplot(
    zc,
    as_data     = TRUE,
    max_samples = max_samples,
    plot_ci     = FALSE,
    extrapolate = TRUE,
    from        = -20,
    to          = 20,
    length.out  = integration_points
  ))
  expected_area     <- mean(zc[["zplot"]][["estimates"]][["weights"]])

  expect_equal(fitted_area, 1, tolerance = 0.01)
  expect_equal(extrapolated_area, expected_area, tolerance = 0.01)
  expect_gt(extrapolated_area, fitted_area)
})

# ============================================================================ #
# Test: PET Models Zplot
# ============================================================================ #

test_that("zplot for positive-direction PET model renders base output", {

  name <- "dat.lehmann2018-PET"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_PET_pos_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "PET Model (Pos) Zplot"))
  })

})

test_that("zplot for negative-direction PET model renders base output", {

  skip_if_not_full_visuals("Negative-direction PET zplot is gallery coverage.")

  name <- "dat.lehmann2018-PET_neg"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_PET_neg_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "PET Model (Neg) Zplot"))
  })

})

# ============================================================================ #
# Test: Multilevel Models Zplot
# ============================================================================ #

test_that("zplot for multilevel model renders base output", {

  skip_if_not_full_visuals("Multilevel zplot is gallery coverage.")

  name <- "konstantopoulos2011_3lvl"
  skip_if_missing_fits(name)

  # Retain the human-reviewed snapshot of the former cluster-conditional
  # default; the new marginal default is certified structurally above.
  zc   <- .test_as_zplot(fits[[name]], conditioning_depth = "cluster")

  expect_vdiffr_snapshot("zplot_multilevel_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "Multilevel Model Zplot"))
  })

})



