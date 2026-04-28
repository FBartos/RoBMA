context("VIF")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

expect_vif_table <- function(x, n_terms = NULL) {

  expect_s3_class(x, "data.frame")
  expect_true(all(c("term", "df", "GVIF", "GVIF^(1/(2*df))") %in% names(x)))
  expect_true(all(nzchar(x[["term"]])))
  expect_true(all(x[["df"]] >= 1))
  expect_true(all(is.finite(x[["GVIF"]])))
  expect_true(all(is.finite(x[["GVIF^(1/(2*df))"]])))
  expect_true(all(x[["GVIF"]] >= 1 - sqrt(.Machine$double.eps)))
  expect_true(all(x[["GVIF^(1/(2*df))"]] >= 1 - sqrt(.Machine$double.eps)))

  if (!is.null(n_terms)) {
    expect_equal(nrow(x), n_terms)
  }
}

brma_term_btt <- function(fit_brma) {

  X      <- RoBMA:::.get_model_matrix(fit_brma)
  assign <- attr(X, "assign")

  if (is.null(assign)) {
    stop("Model matrix is missing the 'assign' attribute.", call. = FALSE)
  }

  term_indices <- sort(unique(assign[assign != 0]))
  btt          <- lapply(term_indices, function(i) which(assign == i))

  return(btt)
}

metafor_vif_value <- function(fit_metafor, btt) {

  out       <- metafor::vif(fit_metafor, btt = btt)
  component <- if (!is.null(out[["beta"]])) out[["beta"]] else out
  row       <- component[["vif"]][[1]]

  return(c(
    GVIF = unname(row[["vif"]]),
    GSIF = unname(row[["sif"]])
  ))
}

metafor_vif_table <- function(fit_brma, fit_metafor, btt = NULL) {

  brma_vif <- vif(fit_brma, posterior_correlation = FALSE)[["vif"]]

  if (is.null(btt)) {
    btt <- brma_term_btt(fit_brma)
  }

  expected <- do.call(rbind, lapply(btt, function(x) metafor_vif_value(fit_metafor, x)))

  return(data.frame(
    term              = brma_vif[["term"]],
    df                = brma_vif[["df"]],
    GVIF              = unname(expected[, "GVIF"]),
    "GVIF^(1/(2*df))" = unname(expected[, "GSIF"]),
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  ))
}

expect_vif_matches_metafor <- function(name, tolerance = 0.10, btt = NULL) {

  fit_brma    <- fits[[name]]
  fit_metafor <- info[[name]][["metafor"]]

  brma_vif    <- vif(fit_brma, posterior_correlation = FALSE)[["vif"]]
  metafor_vif <- metafor_vif_table(fit_brma, fit_metafor, btt = btt)

  expect_vif_table(brma_vif, nrow(metafor_vif))
  expect_equal(brma_vif[["df"]], metafor_vif[["df"]],
               info = paste("brma VIF degrees of freedom should match metafor for", name))
  expect_equal(brma_vif[["GVIF"]], metafor_vif[["GVIF"]], tolerance = tolerance,
               info = paste("brma GVIF should match metafor for", name))
  expect_equal(brma_vif[["GVIF^(1/(2*df))"]], metafor_vif[["GVIF^(1/(2*df))"]],
               tolerance = tolerance,
               info = paste("brma adjusted GVIF should match metafor for", name))
}


# ============================================================================ #
# Test: Normal Meta-Regression VIF
# ============================================================================ #

test_that("VIF for normal meta-regression models matches metafor", {

  model_names <- c(
    "bcg_meta-regression",
    "bcg_meta-regression2",
    "bcg_meta-regression3",
    "bcg_meta-regression4",
    "bangertdrowns2004_location-scale",
    "konstantopoulos2011_3lvl2"
  )
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    expect_vif_matches_metafor(name)
  }
})


# ============================================================================ #
# Test: Publication-Bias Meta-Regression VIF
# ============================================================================ #

test_that("VIF for publication-bias meta-regression models matches metafor", {

  model_names <- c("dat.lehmann2018-PETreg", "dat.lehmann2018-3PSMreg")
  skip_if_missing_fits(model_names)

  # metafor PET reference includes sqrt(vi) before the moderator, while brma
  # appends PET after user-specified moderators.
  expect_vif_matches_metafor("dat.lehmann2018-PETreg", btt = list(3, 2))
  expect_vif_matches_metafor("dat.lehmann2018-3PSMreg")
})


# ============================================================================ #
# Test: GLMM Meta-Regression VIF
# ============================================================================ #

test_that("VIF for GLMM meta-regression models matches metafor", {

  model_names <- c("bcg_glmm_reg")
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    expect_vif_matches_metafor(name)
  }
})


# ============================================================================ #
# Test: Posterior Correlation Diagnostics
# ============================================================================ #

test_that("VIF returns optional posterior correlation diagnostics", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  vif_with_post <- vif(fit_brma)
  post_cor      <- vif_with_post[["posterior_correlation"]]

  expect_s3_class(vif_with_post, "vif.brma")
  expect_vif_table(vif_with_post[["vif"]], 2)
  expect_type(post_cor, "double")
  expect_equal(nrow(post_cor), ncol(post_cor))
  expect_equal(unname(diag(post_cor)), rep(1, nrow(post_cor)), tolerance = 1e-12)
  expect_equal(post_cor, t(post_cor), tolerance = 1e-12)
  expect_true(all(is.finite(post_cor)))

  vif_without_post <- vif(fit_brma, posterior_correlation = FALSE)
  expect_null(vif_without_post[["posterior_correlation"]])
})


# ============================================================================ #
# Test: Intercept-Only Models
# ============================================================================ #

test_that("VIF errors for models without moderators", {

  model_names <- c(
    "bcg_meta-analysis",
    "konstantopoulos2011_3lvl",
    "nielweise2008_glmm",
    "dat.lehmann2018-PET",
    "dat.lehmann2018-3PSM"
  )
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    expect_error(
      vif(fits[[name]]),
      "only meaningful for models with moderators",
      info = paste("vif should error without moderators for", name)
    )
  }
})
