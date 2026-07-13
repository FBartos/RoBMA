context("brma.mv visual diagnostics")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-visuals.R"))

skip_if_no_fits()

mv_visual_fit_names <- c(
  "brma.mv_v14_konstantopoulos2011_cs",
  "brma.mv_v14_assink2016_nested",
  "brma.mv_v14_ishak2007_har",
  "brma.mv_v14_begg1989_study_treatment"
)
mv_visual_fixed_fit_name <- "brma.mv_block_mvn_fixed_random_null"

mv_visual_render_fit_names <- c(
  "brma.mv_v14_ishak2007_har"
)

mv_visual_regplot_cases <- list(
  list(
    name = "brma.mv_v14_ishak2007_har",
    mod  = "time_factor"
  ),
  list(
    name = "brma.mv_v14_begg1989_study_treatment",
    mod  = "trt"
  )
)

fits <- lazy_fits(c(mv_visual_fixed_fit_name, mv_visual_fit_names),
                  validate = FALSE)

.mv_visual_full <- function() {

  toupper(Sys.getenv("ROBMA_TEST_FULL_VISUALS", unset = "")) %in%
    c("TRUE", "1", "YES")
}

.mv_visual_should_render <- function(name) {

  name %in% mv_visual_render_fit_names || .mv_visual_full()
}

.expect_mv_plot_points <- function(plot_data, fit, label) {

  expect_true(is.list(plot_data), info = label)
  expect_true(is.data.frame(plot_data[["points"]]), info = label)
  expect_equal(nrow(plot_data[["points"]]), nobs(fit), info = label)
  expect_true(
    all(c("x", "y") %in% names(plot_data[["points"]])),
    info = label
  )
  expect_true(all(is.finite(plot_data[["points"]][["x"]])), info = label)
  expect_true(all(is.finite(plot_data[["points"]][["y"]])), info = label)
}

.expect_mv_components <- function(plot_data, components, label) {

  expect_true(all(components %in% names(plot_data)), info = label)
}


test_that("v14 brma.mv funnel diagnostics return data and render", {

  skip_if_missing_fits(mv_visual_fit_names)

  expected_components <- c(
    "points", "funnel", "funnel_edge1", "funnel_edge2",
    "background", "x_range", "y_range", "xlim", "ylim"
  )

  for (name in mv_visual_fit_names) {
    fit_brma <- fits[[name]]

    outcome_data <- funnel(
      fit_brma,
      residual    = FALSE,
      as_data     = TRUE,
      max_samples = 100
    )
    residual_data <- funnel(
      fit_brma,
      residual    = TRUE,
      type        = "rstandard",
      as_data     = TRUE,
      max_samples = 100
    )

    .expect_mv_components(outcome_data, expected_components, name)
    .expect_mv_components(residual_data, expected_components, name)
    .expect_mv_plot_points(outcome_data, fit_brma, paste(name, "outcome funnel"))
    .expect_mv_plot_points(residual_data, fit_brma, paste(name, "residual funnel"))

    if (.mv_visual_should_render(name)) {
      .with_temp_plot_device(expect_silent(funnel(
        fit_brma,
        residual    = FALSE,
        max_samples = 100
      )))
      expect_true(.is_ggplot(funnel(
        fit_brma,
        residual    = FALSE,
        plot_type   = "ggplot",
        max_samples = 100
      )))
    }
  }
})

test_that("v14 brma.mv Q-Q diagnostics use rstandard residuals", {

  skip_if_missing_fits(mv_visual_fit_names)

  expected_components <- c(
    "x", "y", "points", "refline", "envelope",
    "xlim", "ylim", "xlab", "ylab"
  )

  for (name in mv_visual_fit_names) {
    fit_brma <- fits[[name]]
    qq_data <- suppressWarnings(qqnorm(fit_brma, as_data = TRUE))
    qq_rstd <- suppressWarnings(qqnorm(
      fit_brma,
      type        = "rstandard",
      as_data     = TRUE,
      max_samples = 100
    ))

    .expect_mv_components(qq_data, expected_components, name)
    .expect_mv_components(qq_rstd, expected_components, name)
    .expect_mv_plot_points(qq_data, fit_brma, paste(name, "Q-Q plot"))
    .expect_mv_plot_points(qq_rstd, fit_brma, paste(name, "rstandard Q-Q plot"))
    expect_equal(qq_data[["x"]], sort(qq_data[["x"]]), info = name)
    expect_equal(qq_data[["y"]], sort(qq_data[["y"]]), info = name)

    if (.mv_visual_should_render(name)) {
      .with_temp_plot_device(expect_silent(suppressWarnings(qqnorm(
        fit_brma,
        type        = "rstandard",
        max_samples = 100
      ))))
      expect_true(.is_ggplot(suppressWarnings(qqnorm(
        fit_brma,
        type        = "rstandard",
        plot_type   = "ggplot",
        max_samples = 100
      ))))
    }
  }
})

test_that("v14 brma.mv forest adapters return data and render", {

  skip_if_missing_fits(mv_visual_fit_names)
  skip_if_not_installed("metafor")

  for (name in mv_visual_fit_names) {
    fit_brma    <- fits[[name]]
    forest_data <- as_metafor_forest(fit_brma, addpred = FALSE)

    expect_true(inherits(forest_data, "metafor_forest.brma"), info = name)
    expect_equal(nrow(forest_data[["studies"]]), nobs(fit_brma), info = name)
    expect_equal(length(forest_data[["forest_args"]][["x"]]), nobs(fit_brma),
                 info = name)
    expect_true(all(is.finite(forest_data[["studies"]][["yi"]])), info = name)
    expect_true(all(is.finite(forest_data[["studies"]][["sei"]])), info = name)
    expect_true(all(is.finite(forest_data[["pooled"]][["estimate"]])),
                info = name)
    expect_true(all(is.finite(forest_data[["pooled"]][["ci.lb"]])),
                info = name)
    expect_true(all(is.finite(forest_data[["pooled"]][["ci.ub"]])),
                info = name)
    expect_false(forest_data[["addpred"]], info = name)

    if (.mv_visual_should_render(name)) {
      .with_temp_plot_device(expect_silent(metafor::forest(
        fit_brma,
        addpred = FALSE
      )))
    }
  }
})

test_that("fixed-effect brma.mv visual diagnostics use zero heterogeneity", {

  skip_if_missing_fits(mv_visual_fixed_fit_name)
  skip_if_not_installed("metafor")

  fit_brma <- fits[[mv_visual_fixed_fit_name]]

  funnel_data <- funnel(
    fit_brma,
    sampling_heterogeneity = TRUE,
    as_data                = TRUE,
    max_samples            = 100
  )
  qq_data <- suppressWarnings(qqnorm(
    fit_brma,
    type        = "rstandard",
    as_data     = TRUE,
    max_samples = 100
  ))
  forest_data <- as_metafor_forest(fit_brma, addpred = TRUE)

  .expect_mv_components(
    funnel_data,
    c(
      "points", "funnel", "funnel_edge1", "funnel_edge2",
      "background", "x_range", "y_range", "xlim", "ylim"
    ),
    mv_visual_fixed_fit_name
  )
  .expect_mv_components(
    qq_data,
    c("x", "y", "points", "refline", "envelope", "xlim", "ylim"),
    mv_visual_fixed_fit_name
  )
  .expect_mv_plot_points(funnel_data, fit_brma, "fixed brma.mv funnel")
  .expect_mv_plot_points(qq_data, fit_brma, "fixed brma.mv Q-Q")
  expect_true(inherits(forest_data, "metafor_forest.brma"))
  expect_equal(nrow(forest_data[["studies"]]), nobs(fit_brma))
  expect_true(all(is.finite(forest_data[["pooled"]][["estimate"]])))

  .with_temp_plot_device(expect_silent(funnel(
    fit_brma,
    sampling_heterogeneity = TRUE,
    max_samples            = 100
  )))
  expect_true(.is_ggplot(funnel(
    fit_brma,
    sampling_heterogeneity = TRUE,
    plot_type              = "ggplot",
    max_samples            = 100
  )))
  .with_temp_plot_device(expect_silent(suppressWarnings(qqnorm(
    fit_brma,
    type        = "rstandard",
    max_samples = 100
  ))))
  expect_true(.is_ggplot(suppressWarnings(qqnorm(
    fit_brma,
    type        = "rstandard",
    plot_type   = "ggplot",
    max_samples = 100
  ))))
  .with_temp_plot_device(expect_silent(metafor::forest(
    fit_brma,
    addpred = TRUE
  )))
})

test_that("v14 brma.mv regplot diagnostics cover moderator models", {

  skip_if_missing_fits(vapply(
    mv_visual_regplot_cases,
    function(case) case[["name"]],
    character(1)
  ))

  expected_components <- c(
    "points", "pred", "ci", "pi", "si",
    "xlim", "ylim", "xlab", "ylab", "mod_type"
  )

  for (case in mv_visual_regplot_cases) {
    fit_brma <- fits[[case[["name"]]]]
    label    <- paste(case[["name"]], "regplot")

    regplot_data <- regplot(
      fit_brma,
      mod         = case[["mod"]],
      pi          = TRUE,
      si          = TRUE,
      as_data     = TRUE,
      max_samples = 100
    )

    .expect_mv_components(regplot_data, expected_components, label)
    .expect_mv_plot_points(regplot_data, fit_brma, label)
    expected_band_rows <- if (regplot_data[["mod_type"]] == "continuous") {
      2 * nrow(regplot_data[["pred"]])
    } else {
      nrow(regplot_data[["pred"]])
    }
    expect_equal(nrow(regplot_data[["ci"]]), expected_band_rows, info = label)
    expect_equal(nrow(regplot_data[["pi"]]), expected_band_rows, info = label)
    expect_equal(nrow(regplot_data[["si"]]), expected_band_rows, info = label)

    .with_temp_plot_device(expect_silent(regplot(
      fit_brma,
      mod         = case[["mod"]],
      pi          = TRUE,
      si          = TRUE,
      max_samples = 100
    )))
    expect_true(.is_ggplot(regplot(
      fit_brma,
      mod         = case[["mod"]],
      pi          = TRUE,
      si          = TRUE,
      plot_type   = "ggplot",
      max_samples = 100
    )))
  }
})

test_that("v14 brma.mv zplot diagnostics return reusable data", {

  skip_if_missing_fits(mv_visual_fit_names)

  for (name in mv_visual_fit_names) {
    fit_brma <- fits[[name]]
    zc       <- as_zplot(fit_brma, max_samples = 100)
    density  <- lines(
      zc,
      as_data     = TRUE,
      max_samples = 100,
      plot_ci     = FALSE,
      from        = -80,
      to          = 80,
      length.out  = 25
    )

    expect_true(inherits(zc, "zplot_brma"), info = name)
    expect_equal(length(zc[["zplot"]][["data"]][["z"]]), nobs(fit_brma),
                 info = name)
    expect_true(all(is.finite(zc[["zplot"]][["estimates"]][["EDR"]])),
                info = name)
    expect_true(all(is.finite(zc[["zplot"]][["estimates"]][["weights"]])),
                info = name)
    expect_true(all(c("x", "y", "y_lCI", "y_uCI") %in% names(density)),
                info = name)
    expect_true(all(is.finite(density[["y"]])), info = name)
  }
})

test_that("v14 brma.mv radial diagnostics are explicitly unsupported", {

  skip_if_missing_fits(mv_visual_fit_names)

  for (name in mv_visual_fit_names) {
    fit_brma <- fits[[name]]

    expect_error(
      radial(fit_brma, as_data = TRUE),
      "Radial plots",
      fixed = TRUE,
      info  = name
    )
  }
})
