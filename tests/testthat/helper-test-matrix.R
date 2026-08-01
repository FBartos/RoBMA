test_tier <- function() {

  tier <- "core"
  if (is_certification_profile()) {
    tier <- c(tier, "extended")
  }

  return(tier)
}

visual_test_tier <- function() {

  tier <- test_tier()
  if (is_certification_profile()) {
    tier <- c(tier, "visual-gallery")
  }

  return(tier)
}

test_profile_value <- function(standard, certification) {

  if (is_certification_profile()) {
    return(certification)
  }

  return(standard)
}

case_value <- function(case, name, default = NULL) {

  if (!name %in% names(case)) {
    return(default)
  }

  value <- case[[name]]
  if (is.list(value)) {
    value <- value[[1]]
  }
  if (length(value) == 0 || (length(value) == 1 && is.na(value))) {
    return(default)
  }

  return(value)
}

case_name <- function(case) {

  return(case_value(case, "name"))
}

case_label <- function(case) {

  label <- case_value(case, "label", default = case_name(case))
  return(label)
}

case_has_check <- function(case, check) {

  return(check %in% case_value(case, "checks", character()))
}

filter_cases <- function(cases, tier = test_tier()) {

  if ("tier" %in% names(cases)) {
    cases <- cases[cases[["tier"]] %in% tier, , drop = FALSE]
  }

  if ("name" %in% names(cases)) {
    catalog_names <- fit_catalog()[["name"]]
    cached_case   <- cases[["name"]] %in% catalog_names
    active_case   <- cases[["name"]] %in% active_fit_catalog()[["name"]]
    cases <- cases[!cached_case | active_case, , drop = FALSE]
  }

  return(cases)
}

for_each_case <- function(cases, callback, tier = test_tier()) {

  cases <- filter_cases(cases, tier = tier)
  if (nrow(cases) == 0L) {
    return(invisible(NULL))
  }

  for (i in seq_len(nrow(cases))) {
    local({
      case <- cases[i, , drop = FALSE]
      callback(case)
    })
  }

  return(invisible(NULL))
}

test_that_case <- function(description, case, code) {

  testthat::test_that(paste0(description, " [", case_label(case), "]"), {
    code
  })
}

hatvalue_cases <- function() {

  data.frame(
    name = c(
      "bcg_meta-analysis",
      "bcg_meta-regression",
      "bcg_meta-regression4",
      "konstantopoulos2011_3lvl",
      "konstantopoulos2011_3lvl2",
      "dat.lehmann2018-PET",
      "dat.lehmann2018-PETreg",
      "dat.lehmann2018-PET_neg",
      "dat.lehmann2018-PEESE",
      "dat.lehmann2018-PEESEreg",
      "dat.lehmann2018-PEESE_neg"
    ),
    label = c(
      "normal simple",
      "normal meta-regression",
      "normal interaction",
      "normal multilevel",
      "normal multilevel meta-regression",
      "PET",
      "PET meta-regression",
      "PET negative",
      "PEESE",
      "PEESE meta-regression",
      "PEESE negative"
    ),
    tolerance = rep(0.05, 11),
    tier = c("core", "core", "core", "core", "extended", "core", "core", "extended",
             "core", "core", "extended"),
    stringsAsFactors = FALSE
  )
}

vif_cases <- function() {

  out <- data.frame(
    name = c(
      "bcg_meta-regression",
      "bcg_meta-regression2",
      "bcg_meta-regression3",
      "bcg_meta-regression4",
      "bangertdrowns2004_location-scale",
      "konstantopoulos2011_3lvl2",
      "dat.lehmann2018-PETreg",
      "dat.lehmann2018-3PSMreg",
      "bcg_glmm_reg"
    ),
    label = c(
      "normal continuous",
      "normal factor",
      "normal interaction year",
      "normal interaction factor",
      "normal location-scale",
      "normal multilevel meta-regression",
      "PET meta-regression",
      "selection meta-regression",
      "GLMM meta-regression"
    ),
    tolerance = rep(0.10, 9),
    tier = c("core", "extended", "extended", "core", "core", "extended", "core", "core", "core"),
    stringsAsFactors = FALSE
  )
  out[["btt"]] <- I(list(NULL, NULL, NULL, NULL, NULL, NULL, list(3, 2), NULL, NULL))

  return(out)
}

dfbetas_metafor_cases <- function() {

  out <- data.frame(
    name = c(
      "bcg_meta-analysis",
      "bcg_meta-regression",
      "bcg_meta-regression4",
      "konstantopoulos2011_3lvl",
      "konstantopoulos2011_3lvl2",
      "dat.lehmann2018-PET",
      "dat.lehmann2018-PETreg",
      "dat.lehmann2018-PET_neg",
      "dat.lehmann2018-PEESE",
      "dat.lehmann2018-PEESEreg",
      "dat.lehmann2018-PEESE_neg"
    ),
    label = c(
      "normal simple",
      "normal meta-regression",
      "normal interaction",
      "normal multilevel",
      "normal multilevel meta-regression",
      "PET",
      "PET meta-regression",
      "PET negative",
      "PEESE",
      "PEESE meta-regression",
      "PEESE negative"
    ),
    tolerance = rep(0.10, 11),
    oracle = c("equal", "equal", "structure", "equal", "equal", "mapped", "mapped", "mapped",
               "mapped", "mapped", "mapped"),
    tier = c("core", "core", "core", "core", "extended", "core", "core", "extended",
             "core", "core", "extended"),
    stringsAsFactors = FALSE
  )
  out[["skip_rows"]] <- I(list(integer(), c(4, 6, 13), integer(), integer(),
                               integer(), integer(), integer(), integer(),
                               integer(), integer(), integer()))
  out[["metafor_cols"]] <- I(list(NULL, NULL, NULL, NULL, NULL, 1, c(1, 3), 1,
                                  1, c(1, 3), 1))
  out[["brma_cols"]]    <- I(list(NULL, NULL, NULL, NULL, NULL, 1, c(1, 2), 1,
                                  1, c(1, 2), 1))

  return(out)
}

summary_heterogeneity_cases <- function() {

  data.frame(
    name = c(
      "bcg_meta-analysis",
      "konstantopoulos2011_3lvl",
      "bcg_meta-regression",
      "bangertdrowns2004_location-scale",
      "dat.lehmann2018-PET",
      "dat.lehmann2018-PET_neg",
      "dat.lehmann2018-PEESE",
      "dat.lehmann2018-PEESEreg",
      "dat.lehmann2018-PEESE_neg",
      "dat.lehmann2018-3PSM",
      "dat.lehmann2018-3PSM_neg",
      "nielweise2008_glmm"
    ),
    label = c(
      "normal simple",
      "normal multilevel",
      "normal meta-regression",
      "normal location-scale",
      "PET",
      "PET negative",
      "PEESE",
      "PEESE meta-regression",
      "PEESE negative",
      "selection",
      "selection negative",
      "GLMM"
    ),
    kind = c("standard", "multilevel", "standard", "scale", "standard",
             "standard", "standard", "standard", "standard", "selection",
             "selection", "standard"),
    tolerance = c(0.05, 0.05, 0.10, 0.05, 0.05, 0.05, 0.05, 0.05,
                  0.05, 0.05, 0.05, 0.10),
    h2_tolerance = c(0.05, NA, 0.20, 0.20, 0.05, 0.05, 0.05, 0.05,
                     0.05, NA, NA, 0.10),
    i2_tolerance = c(0.05, NA, 0.10, 0.05, 0.05, 0.05, 0.05, 0.05,
                     0.05, NA, NA, 0.20),
    tier = c("core", "core", "core", "core", "core", "extended", "core",
             "core", "extended", "core", "extended", "core"),
    stringsAsFactors = FALSE
  )
}

influence_metafor_cases <- function() {

  out <- data.frame(
    name = c(
      "bcg_meta-analysis",
      "bcg_meta-regression",
      "bcg_meta-regression4",
      "konstantopoulos2011_3lvl",
      "konstantopoulos2011_3lvl2"
    ),
    label = c(
      "normal simple",
      "normal meta-regression",
      "normal interaction",
      "normal multilevel",
      "normal multilevel meta-regression"
    ),
    oracle = c("equal", "equal_without_covratio", "rank", "cooks_finite", "cooks_finite"),
    tier = c("core", "core", "core", "core", "extended"),
    stringsAsFactors = FALSE
  )
  out[["skip_rows"]] <- I(list(integer(), c(4, 6, 13), c(4, 5, 6, 8, 10),
                               integer(), integer()))

  return(out)
}

residual_metafor_cases <- function() {

  out <- data.frame(
    name = c(
      "bcg_meta-analysis",
      "bcg_meta-regression",
      "bcg_meta-regression4",
      "bangertdrowns2004_location-scale",
      "konstantopoulos2011_3lvl",
      "konstantopoulos2011_3lvl2",
      "dat.lehmann2018-3PSM",
      "dat.lehmann2018-3PSMreg",
      "dat.lehmann2018-3PSM_neg",
      "dat.lehmann2018-PET",
      "dat.lehmann2018-PETreg",
      "dat.lehmann2018-PET_neg",
      "dat.lehmann2018-PEESE",
      "dat.lehmann2018-PEESEreg",
      "dat.lehmann2018-PEESE_neg",
      "nielweise2008_glmm",
      "bcg_glmm_reg"
    ),
    label = c(
      "normal simple",
      "normal meta-regression",
      "normal interaction",
      "normal location-scale",
      "normal multilevel",
      "normal multilevel meta-regression",
      "selection",
      "selection meta-regression",
      "selection negative",
      "PET",
      "PET meta-regression",
      "PET negative",
      "PEESE",
      "PEESE meta-regression",
      "PEESE negative",
      "GLMM",
      "GLMM meta-regression"
    ),
    kind = c(
      "standard",
      "regression",
      "interaction",
      "standard",
      "multilevel",
      "multilevel_no_loo",
      "selection_pos",
      "selection_reg",
      "selection_neg",
      "pet",
      "pet_reg",
      "pet",
      "peese",
      "peese_reg",
      "peese",
      "glmm",
      "glmm_reg"
    ),
    tolerance = c(0.05, 0.10, 0.10, 0.05, 0.05, 0.05, 0.05,
                  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                  0.05, 0.05, 0.10),
    tier = c("core", "core", "core", "core", "core", "extended", "core",
             "extended", "extended", "core", "core", "extended", "core",
             "core", "extended", "core", "core"),
    stringsAsFactors = FALSE
  )
  out[["rstudent"]] <- I(list(
    "equal", "rank", NULL, NULL, "rank", NULL, "selection_pos", NULL,
    "selection_neg", "equal", NULL, "equal", "equal", NULL, "equal",
    NULL, NULL
  ))

  return(out)
}

prediction_metafor_cases <- function() {

  data.frame(
    name = c(
      "bcg_meta-analysis",
      "bcg_meta-regression",
      "bcg_meta-regression4",
      "bangertdrowns2004_location-scale",
      "konstantopoulos2011_3lvl",
      "konstantopoulos2011_3lvl2",
      "bcg_glmm",
      "bcg_glmm_reg",
      "dat.lehmann2018-3PSM",
      "dat.lehmann2018-3PSMreg",
      "dat.lehmann2018-3PSM_neg",
      "dat.lehmann2018-PET",
      "dat.lehmann2018-PETreg",
      "dat.lehmann2018-PET_neg",
      "dat.lehmann2018-PEESE",
      "dat.lehmann2018-PEESEreg",
      "dat.lehmann2018-PEESE_neg"
    ),
    label = c(
      "normal simple",
      "normal meta-regression",
      "normal interaction",
      "normal location-scale",
      "normal multilevel",
      "normal multilevel meta-regression",
      "GLMM",
      "GLMM meta-regression",
      "selection",
      "selection meta-regression",
      "selection negative",
      "PET",
      "PET meta-regression",
      "PET negative",
      "PEESE",
      "PEESE meta-regression",
      "PEESE negative"
    ),
    kind = c(
      "simple",
      "regression",
      "interaction",
      "scale",
      "multilevel",
      "multilevel",
      "glmm",
      "glmm_reg",
      "selection",
      "selection",
      "selection",
      "pet",
      "pet_reg",
      "pet",
      "peese",
      "peese_reg",
      "peese"
    ),
    tolerance = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.10, 0.05,
                  0.15, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                  0.05, 0.05, 0.05),
    tau_tolerance = c(0.05, NA, NA, 0.05, 0.05, 0.05, 0.05,
                      0.10, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                      0.05, 0.05, 0.05),
    fitted_tolerance = c(NA, 0.05, NA, 0.05, NA, NA, NA, NA, NA,
                         NA, NA, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10),
    tier = c("core", "core", "core", "core", "core", "extended", "core",
             "core", "core", "extended", "extended", "core", "core", "extended",
             "core", "core", "extended"),
    stringsAsFactors = FALSE
  )
}

prediction_newdata_metafor_cases <- function() {

  data.frame(
    name = c(
      "bcg_meta-regression",
      "bcg_meta-regression2",
      "bangertdrowns2004_location-scale",
      "konstantopoulos2011_3lvl2",
      "bcg_glmm_reg",
      "dat.lehmann2018-PETreg",
      "dat.lehmann2018-PEESEreg",
      "dat.lehmann2018-3PSMreg"
    ),
    label = c(
      "normal continuous moderators",
      "normal factor moderator",
      "normal location-scale",
      "normal multilevel moderator",
      "GLMM factor moderator",
      "PET moderator",
      "PEESE moderator",
      "selection moderator"
    ),
    kind = c(
      "normal_continuous",
      "normal_factor",
      "location_scale",
      "multilevel_mod",
      "glmm_factor",
      "pet_reg",
      "peese_reg",
      "selection_reg"
    ),
    tolerance = c(0.10, 0.10, 0.05, 0.05, 0.10, 0.07, 0.07, 0.05),
    tau_tolerance = c(NA, NA, 0.05, NA, NA, NA, NA, NA),
    tier = c("core", "core", "core", "extended", "core", "core", "core",
             "extended"),
    stringsAsFactors = FALSE
  )
}

marginal_means_cases <- function() {

  data.frame(
    name = c(
      "bcg_meta-regression",
      "bcg_meta-regression2",
      "bcg_meta-regression2b",
      "dat.lehmann2018_BMA.norm_mods",
      "dat.lehmann2018_RoBMA_mods"
    ),
    label = c(
      "continuous moderator",
      "factor moderator",
      "factor moderator transformed",
      "BMA moderator",
      "RoBMA moderator"
    ),
    parameter = c(
      "year",
      "alloc",
      "alloc",
      "Preregistered",
      "Preregistered"
    ),
    tier = c("core", "core", "visual-gallery", "core", "core"),
    stringsAsFactors = FALSE
  )
}

marginal_means_interaction_plot_cases <- function() {

  data.frame(
    name = c(
      "bcg_meta-regression3",
      "bcg_meta-regression4",
      "dat.lehmann2018_RoBMA_mods2"
    ),
    label = c(
      "factor by continuous interaction",
      "factor by factor interaction",
      "RoBMA factor by factor interaction"
    ),
    parameter = c(
      "alloc:year",
      "alloc:year_before1969",
      "Preregistered:Gender"
    ),
    tier = c("core", "core", "core"),
    stringsAsFactors = FALSE
  )
}

.announce_existing_visual_snapshots <- function() {

  snapshotter <- getOption("testthat.snapshotter")
  if (is.null(snapshotter) || !snapshotter$is_active() ||
      is.null(snapshotter$file) || !nzchar(snapshotter$file)) {
    return(invisible(FALSE))
  }

  snapshot_dir <- testthat::test_path("_snaps", snapshotter$file)
  if (!dir.exists(snapshot_dir)) {
    return(invisible(FALSE))
  }

  snapshots <- list.files(snapshot_dir)
  snapshots <- snapshots[!grepl(".new.", snapshots, fixed = TRUE)]
  for (snapshot in snapshots) {
    # Equivalent to announce_snapshot_file(), without requiring edition 3.
    snapshotter$announce_file_snapshot(snapshot)
  }

  return(invisible(TRUE))
}

skip_if_not_full_visuals <- function(reason = NULL) {

  if (!is_certification_profile()) {
    # Conditional file snapshots must be announced before skipping or testthat
    # treats their committed baselines as obsolete. Full-visual runs still
    # exercise every snapshot and therefore retain normal stale-file cleanup.
    .announce_existing_visual_snapshots()

    detail <- if (is.null(reason)) "" else paste0(" ", reason)
    testthat::skip(paste0(
      "Skipping extended visual gallery by default.",
      detail,
      " Run the certification profile to include it."
    ))
  }
}
