context("Metafor coefficient references")

source(testthat::test_path("common-functions.R"))

skip_if_no_fits()
skip_if_not_installed("metafor")
skip_if_not_installed("metadat")

coef_fit_names <- c(
  "bcg_meta-analysis",
  "bcg_meta-regression",
  "bcg_meta-regression2",
  "bcg_meta-regression3",
  "bcg_meta-regression4",
  "bangertdrowns2004_location-scale",
  "konstantopoulos2011_3lvl",
  "konstantopoulos2011_3lvl2",
  "bcg_glmm",
  "bcg_glmm_reg",
  "nielweise2008_glmm",
  "dat.lehmann2018-PET",
  "dat.lehmann2018-PET_neg",
  "dat.lehmann2018-PETreg",
  "dat.lehmann2018-PEESE",
  "dat.lehmann2018-3PSM",
  "dat.lehmann2018-4PSM",
  "dat.lehmann2018-3PSM_neg",
  "dat.lehmann2018-3PSMreg"
)

fits <- lazy_fits(coef_fit_names, validate = TRUE)
info <- lazy_infos(coef_fit_names, validate = TRUE)

summary_mean <- function(name, row) {

  return(fits[[name]]$summary[row, "Mean"])
}

metafor_ref <- function(name) {

  return(info[[name]][["metafor"]])
}

expect_summary_mean_matches <- function(name, row, expected, tolerance) {

  expect_equal(
    summary_mean(name, row),
    expected,
    tolerance = tolerance,
    info      = paste0("Summary mean '", row, "' matches metafor for '", name, "'")
  )
}

coef_cases <- list(
  list("bcg_meta-analysis", "mu", function(m) m$beta[[1]], 0.05),
  list("bcg_meta-analysis", "tau", function(m) sqrt(m$tau2), 0.05),
  list("bcg_meta-regression", "(mu) ablat", function(m) m$beta[[2]], 0.10),
  list("bcg_meta-regression", "(mu) year", function(m) m$beta[[3]], 0.10),
  list("bcg_meta-regression2", "(mu) intercept", function(m) m$beta[[1]], 0.05),
  list("bcg_meta-regression2", "(mu) alloc[random]", function(m) m$beta[[2]], 0.15),
  list("bcg_meta-regression2", "(mu) alloc[systematic]", function(m) m$beta[[3]], 0.05),
  list("bcg_meta-regression3", "(mu) year", function(m) m$beta[[4]], 0.15),
  list("bcg_meta-regression3", "(mu) alloc[systematic]:year", function(m) m$beta[[5]], 0.15),
  list("bcg_meta-regression3", "(mu) alloc[random]:year", function(m) m$beta[[6]], 0.15),
  list("bcg_meta-regression4", "(mu) alloc[systematic]", function(m) m$beta[[3]], 0.15),
  list("bangertdrowns2004_location-scale", "(mu) intercept", function(m) m$beta[[1]], 0.05),
  list("bangertdrowns2004_location-scale", "(mu) meta[1]", function(m) m$beta[[2]], 0.05),
  list("bangertdrowns2004_location-scale", "(mu) ni100", function(m) m$beta[[3]], 0.05),
  list("bangertdrowns2004_location-scale", "(log_tau) ni100", function(m) 0.5 * m$alpha[[2]], 0.15),
  list("konstantopoulos2011_3lvl", "mu", function(m) m$beta[[1]], 0.01),
  list("konstantopoulos2011_3lvl", "tau", function(m) sqrt(m$tau2), 0.01),
  list("konstantopoulos2011_3lvl", "rho", function(m) m$rho, 0.05),
  list("konstantopoulos2011_3lvl2", "(mu) intercept", function(m) m$beta[[1]], 0.10),
  list("konstantopoulos2011_3lvl2", "(mu) vi", function(m) m$beta[[2]], 0.10),
  list("konstantopoulos2011_3lvl2", "tau", function(m) sqrt(m$tau2), 0.10),
  list("konstantopoulos2011_3lvl2", "rho", function(m) m$rho, 0.10),
  list("bcg_glmm", "mu", function(m) m$beta[[1]], 0.05),
  list("bcg_glmm", "tau", function(m) sqrt(m$tau2), 0.05),
  list("bcg_glmm_reg", "(mu) intercept", function(m) m$beta[[1]], 0.05),
  list("bcg_glmm_reg", "(mu) alloc[random]", function(m) m$beta[[2]], 0.15),
  list("bcg_glmm_reg", "(mu) alloc[systematic]", function(m) m$beta[[3]], 0.10),
  list("bcg_glmm_reg", "tau", function(m) sqrt(m$tau2), 0.10),
  list("nielweise2008_glmm", "mu", function(m) m$beta[[1]], 0.05),
  list("nielweise2008_glmm", "tau", function(m) sqrt(m$tau2), 0.10),
  list("dat.lehmann2018-PET", "mu", function(m) m$beta[[1]], 0.05),
  list("dat.lehmann2018-PET", "tau", function(m) sqrt(m$tau2), 0.05),
  list("dat.lehmann2018-PET", "PET", function(m) m$beta[[2]], 0.20),
  list("dat.lehmann2018-PET_neg", "mu", function(m) m$beta[[1]], 0.05),
  list("dat.lehmann2018-PET_neg", "tau", function(m) sqrt(m$tau2), 0.05),
  list("dat.lehmann2018-PET_neg", "PET", function(m) -m$beta[[2]], 0.20),
  list("dat.lehmann2018-PETreg", "(mu) intercept", function(m) m$beta[[1]], 0.05),
  list("dat.lehmann2018-PETreg", "(mu) Preregistered[Pre-Registered]", function(m) m$beta[[3]], 0.05),
  list("dat.lehmann2018-PETreg", "PET", function(m) m$beta[[2]], 0.20),
  list("dat.lehmann2018-PETreg", "tau", function(m) sqrt(m$tau2), 0.05),
  list("dat.lehmann2018-PEESE", "mu", function(m) m$beta[[1]], 0.05),
  list("dat.lehmann2018-PEESE", "tau", function(m) sqrt(m$tau2), 0.05),
  list("dat.lehmann2018-PEESE", "PEESE", function(m) m$beta[[2]], 0.05),
  list("dat.lehmann2018-3PSM", "mu", function(m) m$beta[[1]], 0.01),
  list("dat.lehmann2018-3PSM", "tau", function(m) sqrt(m$tau2), 0.01),
  list("dat.lehmann2018-3PSM", "omega[0.025,1]", function(m) m$delta[[2]], 0.10),
  list("dat.lehmann2018-4PSM", "mu", function(m) m$beta[[1]], 0.05),
  list("dat.lehmann2018-4PSM", "tau", function(m) sqrt(m$tau2), 0.05),
  list("dat.lehmann2018-4PSM", "omega[0.025,0.5]", function(m) m$delta[[2]], 0.20),
  list("dat.lehmann2018-4PSM", "omega[0.5,1]", function(m) m$delta[[3]], 0.20),
  list("dat.lehmann2018-3PSM_neg", "mu", function(m) m$beta[[1]], 0.01),
  list("dat.lehmann2018-3PSM_neg", "tau", function(m) sqrt(m$tau2), 0.01),
  list("dat.lehmann2018-3PSM_neg", "omega[0.025,1]", function(m) m$delta[[2]], 0.10),
  list("dat.lehmann2018-3PSMreg", "(mu) intercept", function(m) m$beta[[1]], 0.05),
  list("dat.lehmann2018-3PSMreg", "(mu) Preregistered[Pre-Registered]", function(m) m$beta[[2]], 0.05),
  list("dat.lehmann2018-3PSMreg", "tau", function(m) sqrt(m$tau2), 0.05),
  list("dat.lehmann2018-3PSMreg", "omega[0.025,1]", function(m) m$delta[[2]], 0.10)
)

test_that("cached summary means match metafor fit references", {

  for (case in coef_cases) {
    name      <- case[[1]]
    row       <- case[[2]]
    expected  <- case[[3]](metafor_ref(name))
    tolerance <- case[[4]]

    expect_summary_mean_matches(name, row, expected, tolerance)
  }
})

test_that("location-scale intercept matches metafor on the log-tau scale", {

  name <- "bangertdrowns2004_location-scale"

  expect_equal(
    log(summary_mean(name, "(log_tau) exp(intercept)")),
    0.5 * metafor_ref(name)$alpha[[1]],
    tolerance = 0.15,
    info      = paste0("Log-tau intercept matches metafor for '", name, "'")
  )
})

test_that("scaled predictor oracle for BCG meta-regression matches JAGS estimates", {

  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(
    measure = "RR",
    ai      = tpos,
    bi      = tneg,
    ci      = cpos,
    di      = cneg,
    data    = dat.bcg
  )

  fit_metafor <- metafor::rma(yi, vi, mods = ~ scale(ablat) + scale(year), data = dat)
  fit_brma    <- BayesTools::JAGS_estimates_table(fits[["bcg_meta-regression"]]$fit)

  expect_equal(fit_metafor$beta[[1]], fit_brma["(mu) intercept", "Mean"], tolerance = 0.10)
  expect_equal(fit_metafor$beta[[2]], fit_brma["(mu) ablat", "Mean"], tolerance = 0.10)
  expect_equal(fit_metafor$beta[[3]], fit_brma["(mu) year", "Mean"], tolerance = 0.10)
})
