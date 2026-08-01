context("Forest plot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))

test_that("dense metafor diamonds retain their four visual corners", {

  file <- tempfile(fileext = ".svg")
  on.exit(unlink(file), add = TRUE)
  writeLines(
    "<polygon points='1,2 3,4 3,4 5,6 5,6 7,8 7,8 1,2 ' style='fill: #000000;' />",
    file
  )

  .canonicalize_dense_diamond_polygons(file, vertices_per_edge = 2L)

  expect_identical(
    readLines(file, warn = FALSE),
    "<polygon points='1,2 3,4 5,6 7,8 ' style='fill: #000000;' />"
  )
})


test_that("dense polygon canonicalization preserves changed geometry", {

  file <- tempfile(fileext = ".svg")
  on.exit(unlink(file), add = TRUE)
  polygon <- paste0(
    "<polygon points='",
    "1,2 2,3 3,4 ",
    "3,4 4.5,9 5,6 ",
    "5,6 6,7 7,8 ",
    "7,8 4,5 1,2 ",
    "' style='fill: #000000;' />"
  )
  writeLines(polygon, file)

  .canonicalize_dense_diamond_polygons(file, vertices_per_edge = 3L)

  expect_identical(readLines(file, warn = FALSE), polygon)
})


test_that("RoBMA does not export a competing forest generic", {

  expect_false("forest" %in% getNamespaceExports("RoBMA"))
  expect_false("metafor" %in% names(getNamespaceImports("RoBMA")))
})


test_that("forest level normalization accepts documented percent values", {

  expect_equal(.forest_normalize_level(95), .95)
  expect_equal(.forest_normalize_level(.95), .95)
  expect_equal(.forest_normalize_level(1), .01)
  expect_error(
    .forest_normalize_level(100),
    "greater than 0 and less than 100",
    fixed = TRUE
  )
})


test_that("forest detects compiled random-slope ambiguity", {

  dat <- data.frame(
    yi    = c(.1, .2, .3),
    x     = c(0, 1, 2),
    study = c("s1", "s2", "s3")
  )
  intercept <- brma.mv(
    yi                        = yi,
    V                         = diag(c(.04, .05, .06)),
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  slope <- brma.mv(
    yi                        = yi,
    V                         = diag(c(.04, .05, .06)),
    random                    = ~ us(1 + x | study),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  slope_only <- brma.mv(
    yi                        = yi,
    V                         = diag(c(.04, .05, .06)),
    random                    = ~ us(0 + x | study),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_false(.forest_prediction_design_is_ambiguous(intercept))
  expect_true(.forest_prediction_design_is_ambiguous(slope))
  expect_true(.forest_prediction_design_is_ambiguous(slope_only))
})

# list cached fits lazily
skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

.test_forest <- function(x, ..., addpred = FALSE) {

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
  expect_equal(as.numeric(out[["forest_args"]][["x"]]), out[["studies"]][["yi"]])
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
      c("estimate", "ci.lb", "ci.ub")
    ])))
  )
  expect_true(all(is.finite(unlist(out[["prediction"]]))))
  expect_true(all(is.na(out[["pooled"]][c("pi.lb", "pi.ub")])))
  expect_equal(
    out[["addpoly_args"]][["x"]],
    out[["prediction"]][["estimate"]]
  )
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
    out <- as_metafor_forest(fit, addpred = FALSE)

    expect_true(inherits(out, "metafor_forest.brma"), info = name)
    expect_equal(nrow(out[["studies"]]), nobs(fit), info = name)
    expect_equal(length(out[["forest_args"]][["x"]]), nobs(fit), info = name)
    expect_equal(length(out[["forest_args"]][["sei"]]), nobs(fit), info = name)
    expect_true(all(is.finite(out[["studies"]][["yi"]])), info = name)
    expect_true(all(is.finite(out[["studies"]][["sei"]])), info = name)
    expect_true(all(is.finite(out[["pooled"]][["estimate"]])), info = name)
    expect_true(all(is.finite(out[["pooled"]][["ci.lb"]])), info = name)
    expect_true(all(is.finite(out[["pooled"]][["ci.ub"]])), info = name)
    expect_null(out[["prediction"]], info = name)
    expect_false(out[["addpred"]], info = name)
  }
})

test_that("forest prediction uses one explicit moderator row throughout", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  fit     <- fits[[name]]
  newdata <- data.frame(ablat = 70, year = 1980)
  expected_terms <- .forest_summary_row(predict(
    fit,
    newdata = newdata,
    type    = "terms",
    quiet   = TRUE
  ))
  set.seed(91)
  expected_effect <- .forest_summary_row(predict(
    fit,
    newdata = newdata,
    type    = "estimate",
    quiet   = TRUE
  ))
  set.seed(91)
  out <- as_metafor_forest(
    fit,
    newdata = newdata,
    addpred = TRUE
  )

  expect_equal(out[["prediction"]][["estimate"]], expected_terms[["estimate"]])
  expect_equal(out[["prediction"]][["ci.lb"]], expected_terms[["ci.lb"]])
  expect_equal(out[["prediction"]][["ci.ub"]], expected_terms[["ci.ub"]])
  expect_equal(out[["prediction"]][["pi.lb"]], expected_effect[["ci.lb"]])
  expect_equal(out[["prediction"]][["pi.ub"]], expected_effect[["ci.ub"]])
  expect_equal(out[["addpoly_args"]][["x"]], expected_terms[["estimate"]])
  expect_identical(out[["addpoly_args"]][["mlab"]], "Predicted Effect")
  expect_error(
    as_metafor_forest(fit, addpred = TRUE),
    "requires one explicit 'newdata' row",
    fixed = TRUE
  )
  expect_error(
    as_metafor_forest(
      fit,
      newdata = rbind(newdata, newdata),
      addpred = TRUE
    ),
    "exactly one 'newdata' row",
    fixed = TRUE
  )
})

test_that("forest known-R prediction requires supported explicit levels", {

  name <- "brma.mv_block_mvn_known_R"
  skip_if_missing_fits(name)

  fit     <- fits[[name]]
  newdata <- data.frame(x = 0, z = -1, study = "s1")
  set.seed(92)
  out <- as_metafor_forest(fit, newdata = newdata, addpred = TRUE)

  expect_true(all(is.finite(unlist(out[["prediction"]]))))
  expect_identical(out[["addpoly_args"]][["mlab"]], "Predicted Effect")

  newdata[["study"]] <- "unseen"
  expect_error(
    as_metafor_forest(fit, newdata = newdata, addpred = TRUE),
    "R_new"
  )
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
  set.seed(1)
  out_dispatch <- metafor::forest(fit, addpred = TRUE, as_data = TRUE)
  set.seed(1)
  expect_identical(out_dispatch, as_metafor_forest(fit, addpred = TRUE))
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

})


test_that("Forest prediction-shade gallery is stable", {

  skip_if_not_full_visuals(
    "Metafor's high-resolution prediction gradient is release-only coverage."
  )
  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

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
