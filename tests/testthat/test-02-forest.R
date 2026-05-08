context("Forest plot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))

test_that("RoBMA does not export a competing forest generic", {

  expect_false("forest" %in% getNamespaceExports("RoBMA"))
  expect_false("metafor" %in% names(getNamespaceImports("RoBMA")))
})

# list cached fits lazily
skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

.test_forest <- function(x, ..., addpred = TRUE) {

  suppressWarnings(metafor::forest(x, addpred = addpred, ...))
}

.forest_snapshot_name <- function(name) {

  name <- gsub("[^[:alnum:]]+", "_", name)
  name <- gsub("_+$", "", name)

  return(paste0("forest_", name))
}

.forest_core_names <- c(
  "bcg_meta-analysis",
  "bcg_meta-regression",
  "konstantopoulos2011_3lvl",
  "bcg_glmm",
  "dat.lehmann2018-PET",
  "dat.lehmann2018-PEESE",
  "dat.lehmann2018-3PSM",
  "dat.lehmann2018_BMA.norm",
  "bcg_BMA.glmm",
  "dat.lehmann2018_RoBMA"
)

.forest_extended_names <- c(
  "bcg_meta-regression2",
  "bangertdrowns2004_location-scale",
  "konstantopoulos2011_3lvl2",
  "bcg_glmm_reg",
  "nielweise2008_glmm",
  "dat.lehmann2018-PET_neg",
  "dat.lehmann2018-PETreg",
  "dat.lehmann2018-PEESE_neg",
  "dat.lehmann2018-PEESEreg",
  "dat.lehmann2018-3PSM_neg",
  "dat.lehmann2018-3PSMreg",
  "dat.lehmann2018_BMA.norm_mods",
  "nielweise2008_BMA.glmm",
  "dat.lehmann2018_RoBMA_mods",
  "dat.lehmann2018_RoBMA_3lvl_mods_scale"
)

.forest_bcg_count_columns <- function(name) {

  dat <- info[[name]][["metafor"]][["data"]]
  return(as.matrix(dat[, c("tpos", "tneg", "cpos", "cneg")]))
}

.forest_auxiliary_columns <- function(fit) {

  columns <- cbind(
    seq_len(nobs(fit)),
    round(.outcome_data_sei(fit), 2)
  )
  colnames(columns) <- c("No.", "SE")

  return(columns)
}


# ============================================================================ #
# Test: Adapter Data
# ============================================================================ #

test_that("as_metafor_forest returns metafor-ready vector arguments", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit <- fits[[name]]
  out <- as_metafor_forest(fit, addpred = TRUE)

  expect_s3_class(out, "metafor_forest.brma")
  expect_equal(nrow(out[["studies"]]), nobs(fit))
  expect_equal(
    out[["studies"]][["yi"]],
    .outcome_data_yi(fit),
    tolerance = 1e-12
  )
  expect_equal(
    out[["studies"]][["sei"]],
    .outcome_data_sei(fit),
    tolerance = 1e-12
  )
  expect_equal(out[["forest_args"]][["x"]], out[["studies"]][["yi"]])
  expect_equal(out[["forest_args"]][["sei"]], out[["studies"]][["sei"]])
  expect_equal(attributes(out[["forest_args"]][["x"]])[["measure"]], .measure(fit))
  expect_equal(
    attributes(out[["forest_args"]][["x"]])[["slab"]],
    out[["studies"]][["slab"]]
  )
  expect_true(
    all(c("estimate", "ci.lb", "ci.ub", "pi.lb", "pi.ub") %in%
      names(out[["pooled"]]))
  )
  expect_true(
    all(is.finite(unlist(out[["pooled"]][
      1,
      c("estimate", "ci.lb", "ci.ub", "pi.lb", "pi.ub")
    ])))
  )
  expect_equal(out[["addpoly_args"]][["x"]], out[["pooled"]][["estimate"]])
  expect_true(is.finite(attributes(out[["addpoly_args"]][["pi.lb"]])[["se"]]))
  expect_identical(out[["addpoly_args"]][["rows"]], -1)
  expect_identical(out[["addpoly_args"]][["predstyle"]], "bar")
  expect_identical(out[["addpoly_args"]][["mlab"]], "Pooled Effect")
  expect_true(out[["default_ylim"]])
  expect_true(out[["forest_args"]][["ylim"]][1] <= -3)
})

test_that("as_metafor_forest preserves metafor auxiliary column arguments", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit          <- fits[[name]]
  columns      <- .forest_auxiliary_columns(fit)
  column_xpos  <- c(-4, -3)
  study_slab   <- paste("Trial", seq_len(nobs(fit)))

  out <- as_metafor_forest(
    fit,
    addpred   = TRUE,
    slab      = study_slab,
    ilab      = columns,
    ilab.lab  = c("No.", "SE"),
    ilab.xpos = column_xpos,
    ilab.pos  = 2,
    xlim      = c(-5, 4),
    atransf   = exp
  )

  expect_identical(out[["forest_args"]][["slab"]], study_slab)
  expect_identical(out[["forest_args"]][["ilab"]], columns)
  expect_identical(out[["forest_args"]][["ilab.lab"]], c("No.", "SE"))
  expect_identical(out[["forest_args"]][["ilab.xpos"]], column_xpos)
  expect_identical(out[["forest_args"]][["ilab.pos"]], 2)
  expect_identical(out[["forest_args"]][["xlim"]], c(-5, 4))
  expect_identical(out[["addpoly_args"]][["atransf"]], exp)

  out_ylim <- as_metafor_forest(fit, ylim = c(-5, 20))

  expect_identical(out_ylim[["forest_args"]][["ylim"]], c(-5, 20))
  expect_false(out_ylim[["default_ylim"]])
})

test_that("as_metafor_forest supports representative fitted object families", {

  skip_if_missing_fits(.forest_core_names)

  for (name in .forest_core_names) {
    fit <- fits[[name]]
    out <- as_metafor_forest(fit, addpred = TRUE)

    expect_s3_class(out, "metafor_forest.brma", info = name)
    expect_equal(nrow(out[["studies"]]), nobs(fit), info = name)
    expect_equal(length(out[["forest_args"]][["x"]]), nobs(fit), info = name)
    expect_equal(length(out[["forest_args"]][["sei"]]), nobs(fit), info = name)
    expect_true(all(is.finite(out[["studies"]][["yi"]])), info = name)
    expect_true(all(is.finite(out[["studies"]][["sei"]])), info = name)
    expect_true(all(is.finite(out[["pooled"]][["estimate"]])), info = name)
    expect_true(all(is.finite(out[["pooled"]][["ci.lb"]])), info = name)
    expect_true(all(is.finite(out[["pooled"]][["ci.ub"]])), info = name)
    expect_true(all(is.finite(out[["pooled"]][["pi.lb"]])), info = name)
    expect_true(all(is.finite(out[["pooled"]][["pi.ub"]])), info = name)
  }
})


# ============================================================================ #
# Test: S3 Dispatch
# ============================================================================ #

test_that("forest methods dispatch through metafor", {

  skip_if_not_installed("metafor")

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit <- fits[[name]]
  out <- as_metafor_forest(fit, addpred = TRUE)

  expect_false("forest" %in% getNamespaceExports("RoBMA"))
  expect_identical(metafor::forest(fit, addpred = TRUE, as_data = TRUE), out)
  .with_temp_plot_device(expect_silent(metafor::forest(fit, addpred = TRUE)))
  .with_temp_plot_device(expect_silent(metafor::forest(out, predstyle = "bar")))
  .with_temp_plot_device(expect_silent({
    do.call(metafor::forest.default, out[["forest_args"]])
    do.call(metafor::addpoly.default, out[["addpoly_args"]])
  }))
  .with_temp_plot_device(expect_silent(
    metafor::forest(
      fit,
      addfit = FALSE,
      slab   = paste0("Study ", seq_len(nobs(fit)))
    )
  ))

  .with_temp_plot_device(expect_silent(metafor::forest(
    as_metafor_forest(
      fit,
      addpred   = TRUE,
      ilab      = .forest_bcg_count_columns(name),
      ilab.lab  = c("T+", "T-", "C+", "C-"),
      ilab.xpos = c(-8.5, -7.5, -6.5, -5.5),
      xlim      = c(-10, 4)
    ),
    predstyle = "bar"
  )))
})

test_that("metafor forest validates brma-specific arguments", {

  skip_if_not_installed("metafor")

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit <- fits[[name]]

  expect_error(
    metafor::forest(fit, addfit = FALSE, addpred = TRUE),
    "requires 'addfit = TRUE'"
  )
  expect_error(
    metafor::forest(fit, predstyle = "bar"),
    "require 'addpred = TRUE'"
  )
  expect_error(
    as_metafor_forest(fit, predstyle = "bar"),
    "controlled by the RoBMA forest adapter"
  )
})

test_that("metafor forest adapter validates adapter-specific arguments", {

  skip_if_not_installed("metafor")

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit <- fits[[name]]

  expect_error(
    metafor::forest(as_metafor_forest(fit), predstyle = "bar"),
    "created with 'addpred = TRUE'"
  )
  expect_error(
    metafor::forest(as_metafor_forest(fit), addpred = TRUE),
    "controlled by the RoBMA forest adapter"
  )
})


# ============================================================================ #
# Test: Visual Coverage
# ============================================================================ #

test_that("Forest plot for simple meta-analysis matches metafor structure", {

  skip_if_not_installed("metafor")

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  expect_vdiffr_snapshot("forest_simple_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::forest(fit_metafor, main = "metafor", addpred = TRUE, atransf = exp)
    metafor::forest(fit_brma,    main = "brma",    addpred = TRUE, atransf = exp)
  })

  expect_vdiffr_snapshot("forest_simple_prediction_bar", function() {
    metafor::forest(
      as_metafor_forest(fit_brma, addpred = TRUE),
      predstyle = "bar",
      atransf   = exp
    )
  })

  expect_vdiffr_snapshot("forest_simple_prediction_shade", function() {
    metafor::forest(
      as_metafor_forest(fit_brma, addpred = TRUE),
      predstyle = "shade",
      atransf   = exp
    )
  })

  expect_vdiffr_snapshot("forest_simple_ilab_columns", function() {
    dat <- info[[name]][["metafor"]][["data"]]
    metafor::forest(
      as_metafor_forest(
        fit_brma,
        addpred   = TRUE,
        slab      = paste(dat[["author"]], dat[["year"]]),
        ilab      = .forest_bcg_count_columns(name),
        ilab.lab  = c("T+", "T-", "C+", "C-"),
        ilab.xpos = c(-8.5, -7.5, -6.5, -5.5),
        xlim      = c(-10, 4),
        alim      = c(-3, 2),
        at        = log(c(0.05, 0.25, 1, 4)),
        atransf   = exp,
        header    = c("Study", "RR [95% CI]")
      ),
      predstyle = "bar"
    )
  })
})

test_that("Forest plot renders representative fitted object families", {

  skip_if_not_installed("metafor")
  skip_if_missing_fits(.forest_core_names)

  for (name in setdiff(.forest_core_names, "bcg_meta-analysis")) {
    local({
      name_local <- name
      fit_brma   <- fits[[name_local]]

      expect_vdiffr_snapshot(.forest_snapshot_name(name_local), function() {
        .test_forest(fit_brma, main = name_local)
      })
    })
  }
})

test_that("Forest plot renders extended visual gallery", {

  skip_if_not_installed("metafor")
  skip_if_not_full_visuals("Forest variants across secondary model objects.")
  skip_if_missing_fits(.forest_extended_names)

  for (name in .forest_extended_names) {
    local({
      name_local <- name
      fit_brma   <- fits[[name_local]]

      expect_vdiffr_snapshot(.forest_snapshot_name(name_local), function() {
        .test_forest(fit_brma, main = name_local)
      })
    })
  }
})
