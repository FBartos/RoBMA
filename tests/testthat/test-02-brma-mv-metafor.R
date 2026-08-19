context("brma.mv metafor references")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-metafor.R"))

skip_if_no_fits()
skip_if_not_installed("metafor")

mv_fixed_metafor_fit_name <- "brma.mv_block_mvn_fixed_random_null"
mv_known_r_metafor_fit_name <- "brma.mv_block_mvn_known_R"
mv_metafor_fit_names <- c(
  "brma.mv_v14_konstantopoulos2011_cs",
  "brma.mv_v14_assink2016_nested",
  "brma.mv_v14_ishak2007_har",
  "brma.mv_v14_begg1989_study_treatment"
)

.mv_active_metafor_fit_names <- function(names = mv_metafor_fit_names) {

  names <- intersect(names, active_fit_catalog()[["name"]])
  if (length(names) == 0L) {
    testthat::skip("No fixtures for this metafor parity check are active.")
  }

  return(names)
}

fits <- lazy_fits(c(mv_fixed_metafor_fit_name, mv_known_r_metafor_fit_name,
                    mv_metafor_fit_names),
                  validate = FALSE)
info <- lazy_infos(c(mv_fixed_metafor_fit_name, mv_known_r_metafor_fit_name,
                     mv_metafor_fit_names),
                   validate = FALSE)

.mv_metafor <- function(name) {

  info[[name]][["metafor"]]
}

.mv_summary_mean <- function(name, section, row) {

  table <- summary(
    fits[[name]],
    include_mcmc_diagnostics = FALSE
  )[[section]]

  table[row, "Mean"]
}

.mv_heterogeneity_mean <- function(name, component, row) {

  heterogeneity <- summary_heterogeneity(
    fits[[name]],
    component = component
  )

  heterogeneity[["estimates"]][row, "Mean"]
}

.expect_close_abs <- function(observed, expected, tolerance, label) {

  expect_true(
    is.finite(observed) && is.finite(expected) &&
      abs(observed - expected) <= tolerance,
    info = paste0(
      label,
      " | observed = ", signif(observed, 6),
      ", expected = ", signif(expected, 6),
      ", abs diff = ", signif(abs(observed - expected), 6),
      ", tolerance = ", tolerance
    )
  )
}

.metafor_normal_prior_beta <- function(fit_metafor, mean, sd) {

  beta                 <- as.numeric(fit_metafor[["beta"]])
  likelihood_precision <- solve(fit_metafor[["vb"]])
  prior_mean           <- rep(mean, length.out = length(beta))
  prior_sd             <- rep(sd, length.out = length(beta))
  prior_precision      <- diag(1 / prior_sd^2)

  return(as.numeric(solve(
    likelihood_precision + prior_precision,
    likelihood_precision %*% beta + prior_precision %*% prior_mean
  )))
}

.metafor_fixed_prediction <- function(fit_metafor,
                                      beta = fit_metafor[["beta"]]) {

  as.vector(fit_metafor[["X"]] %*% beta)
}

.metafor_marginal_residual <- function(fit_metafor,
                                       beta = fit_metafor[["beta"]]) {

  as.numeric(
    fit_metafor[["yi"]] - .metafor_fixed_prediction(fit_metafor, beta)
  )
}

.ishak_prior_aware_beta <- function(fit_metafor) {

  prior <- fits[["brma.mv_v14_ishak2007_har"]][["priors"]][["mods"]][[
    "time_factor"
  ]]
  return(.metafor_normal_prior_beta(
    fit_metafor = fit_metafor,
    mean         = prior[["parameters"]][["mean"]],
    sd           = prior[["parameters"]][["sd"]]
  ))
}

.brma_fixed_prediction_mean <- function(fit_brma) {

  colMeans(as.matrix(predict(
    fit_brma,
    newdata = NULL,
    type    = "terms",
    quiet   = TRUE
  )))
}

.mv_sample_means <- function(samples) {

  colMeans(as.matrix(samples))
}

.mv_posterior_subset <- function(fit_brma, n = 600) {

  posterior_samples <- as.matrix(fit_brma[["fit"]][["mcmc"]])
  rows <- .nested_srs_rows(seq_len(nrow(posterior_samples)), n)

  posterior_samples[rows, , drop = FALSE]
}

test_that("fixed-effect brma.mv with random = NULL matches metafor", {

  skip_if_missing_fits(mv_fixed_metafor_fit_name)

  name        <- mv_fixed_metafor_fit_name
  fit_brma    <- fits[[name]]
  fit_metafor <- .mv_metafor(name)

  brma_pred    <- .brma_fixed_prediction_mean(fit_brma)
  metafor_pred <- .metafor_fixed_prediction(fit_metafor)
  brma_blup    <- .mv_sample_means(blup(fit_brma))

  expect_false(.is_random(fit_brma))
  .expect_close_abs(
    observed  = .mv_summary_mean(name, "estimates", 1L),
    expected  = as.numeric(fit_metafor[["beta"]]),
    tolerance = 0.05,
    label     = paste(name, "fixed effect")
  )
  .expect_close_abs(
    observed  = max(abs(brma_pred - metafor_pred)),
    expected  = 0,
    tolerance = 0.05,
    label     = paste(name, "fixed fitted values")
  )
  expect_equal(
    unname(fitted(fit_brma)),
    unname(brma_pred),
    tolerance = 1e-12
  )
  expect_equal(
    unname(fitted(fit_brma, conditioning_depth = "estimate")),
    unname(brma_blup),
    tolerance = 1e-12
  )
  expect_equal(
    unname(brma_blup),
    unname(.mv_sample_means(true_effects(fit_brma))),
    tolerance = 1e-12
  )
})


test_that("known-R brma.mv rstandard and VIF track metafor references", {

  skip_if_missing_fits(mv_known_r_metafor_fit_name)

  name        <- mv_known_r_metafor_fit_name
  fit_brma    <- fits[[name]]
  fit_metafor <- .mv_metafor(name)

  brma_standard    <- rstandard(fit_brma)
  metafor_standard <- metafor::rstandard.rma.mv(fit_metafor)

  expect_residual_table(brma_standard, nobs(fit_brma), info = name)
  expect_equal(nrow(brma_standard), length(metafor_standard[["z"]]), info = name)
  expect_true(
    stats::cor(brma_standard[["z"]], metafor_standard[["z"]],
               method = "spearman", use = "complete.obs") > 0.85,
    info = paste(name, "rstandard z rank correlation")
  )
  .expect_close_abs(
    observed  = max(abs(brma_standard[["z"]] - metafor_standard[["z"]])),
    expected  = 0,
    tolerance = 0.30,
    label     = paste(name, "maximum rstandard z difference")
  )

  brma_vif    <- vif(fit_brma, posterior_correlation = FALSE)[["vif"]]
  metafor_vif <- metafor_vif_table(fit_brma, fit_metafor)

  expect_vif_table(brma_vif, nrow(metafor_vif), info = name)
  expect_equal(brma_vif[["df"]], metafor_vif[["df"]], info = name)
  expect_equal(
    brma_vif[["GVIF"]],
    metafor_vif[["GVIF"]],
    tolerance = 0.20,
    info      = paste(name, "known-R GVIF")
  )
  expect_equal(
    brma_vif[["GVIF^(1/(2*df))"]],
    metafor_vif[["GVIF^(1/(2*df))"]],
    tolerance = 0.20,
    info      = paste(name, "known-R adjusted GVIF")
  )
})


test_that("v14 brma.mv fixed effects match metafor references", {

  fit_names <- .mv_active_metafor_fit_names()
  skip_if_missing_fits(fit_names)

  cases <- list(
    list(
      name      = "brma.mv_v14_konstantopoulos2011_cs",
      rows      = "intercept",
      expected  = function(m) as.numeric(m[["beta"]]),
      tolerance = 0.05
    ),
    list(
      name      = "brma.mv_v14_assink2016_nested",
      rows      = "intercept",
      expected  = function(m) as.numeric(m[["beta"]]),
      tolerance = 0.05
    ),
    list(
      name      = "brma.mv_v14_ishak2007_har",
      rows      = paste0("time_factor[", 1:4, "]"),
      expected  = .ishak_prior_aware_beta,
      tolerance = 0.40
    ),
    list(
      name      = "brma.mv_v14_begg1989_study_treatment",
      rows      = c("intercept", "trt[BMT]"),
      expected  = function(m) as.numeric(m[["beta"]]),
      tolerance = 0.05
    )
  )
  cases <- Filter(function(case) case[["name"]] %in% fit_names, cases)

  for (case in cases) {
    name     <- case[["name"]]
    observed <- vapply(
      case[["rows"]],
      function(row) .mv_summary_mean(name, "estimates_mods", row),
      numeric(1)
    )
    expected <- case[["expected"]](.mv_metafor(name))

    for (i in seq_along(observed)) {
      .expect_close_abs(
        observed  = observed[[i]],
        expected  = expected[[i]],
        tolerance = case[["tolerance"]],
        label     = paste(name, case[["rows"]][[i]], "fixed effect")
      )
    }
  }
})

test_that("v14 brma.mv heterogeneity components match metafor references", {

  fit_names <- .mv_active_metafor_fit_names()
  skip_if_missing_fits(fit_names)

  cases <- list(
    list(
      name      = "brma.mv_v14_konstantopoulos2011_cs",
      component = "district",
      row       = "tau",
      expected  = function(m) sqrt(m[["tau2"]]),
      tolerance = 0.05
    ),
    list(
      name      = "brma.mv_v14_assink2016_nested",
      component = "study/esid",
      row       = "sd_total",
      expected  = function(m) sqrt(sum(m[["sigma2"]])),
      tolerance = 0.05
    ),
    list(
      name      = "brma.mv_v14_assink2016_nested",
      component = "study",
      row       = "tau",
      expected  = function(m) sqrt(m[["sigma2"]][[1]]),
      tolerance = 0.05
    ),
    list(
      name      = "brma.mv_v14_assink2016_nested",
      component = "esid_study",
      row       = "tau",
      expected  = function(m) sqrt(m[["sigma2"]][[2]]),
      tolerance = 0.05
    ),
    list(
      name      = "brma.mv_v14_assink2016_nested",
      component = "study/esid",
      row       = "var_prop(study)",
      expected  = function(m) m[["sigma2"]][[1]] / sum(m[["sigma2"]]),
      tolerance = 0.08
    ),
    list(
      name      = "brma.mv_v14_assink2016_nested",
      component = "study/esid",
      row       = "var_prop(esid_study)",
      expected  = function(m) m[["sigma2"]][[2]] / sum(m[["sigma2"]]),
      tolerance = 0.08
    ),
    list(
      name      = "brma.mv_v14_ishak2007_har",
      component = "study",
      row       = "sd(time[1])",
      expected  = function(m) sqrt(m[["tau2"]][[1]]),
      tolerance = 0.75
    ),
    list(
      name      = "brma.mv_v14_ishak2007_har",
      component = "study",
      row       = "sd(time[2])",
      expected  = function(m) sqrt(m[["tau2"]][[2]]),
      tolerance = 0.75
    ),
    list(
      name      = "brma.mv_v14_ishak2007_har",
      component = "study",
      row       = "sd(time[3])",
      expected  = function(m) sqrt(m[["tau2"]][[3]]),
      tolerance = 0.75
    ),
    list(
      name      = "brma.mv_v14_ishak2007_har",
      component = "study",
      row       = "sd(time[4])",
      expected  = function(m) sqrt(m[["tau2"]][[4]]),
      tolerance = 1.10
    ),
    list(
      name      = "brma.mv_v14_begg1989_study_treatment",
      component = "study",
      row       = "tau",
      expected  = function(m) sqrt(m[["sigma2"]][[1]]),
      tolerance = 0.06
    ),
    list(
      name      = "brma.mv_v14_begg1989_study_treatment",
      component = "treatment",
      row       = "tau",
      expected  = function(m) sqrt(m[["tau2"]]),
      tolerance = 0.06
    )
  )
  cases <- Filter(function(case) case[["name"]] %in% fit_names, cases)

  for (case in cases) {
    name     <- case[["name"]]
    observed <- .mv_heterogeneity_mean(name, case[["component"]], case[["row"]])
    expected <- case[["expected"]](.mv_metafor(name))

    .expect_close_abs(
      observed  = observed,
      expected  = expected,
      tolerance = case[["tolerance"]],
      label     = paste(name, case[["component"]], case[["row"]])
    )
  }
})

test_that("v14 brma.mv heterogeneity component selectors expose expected names", {

  fit_names <- .mv_active_metafor_fit_names()
  skip_if_missing_fits(fit_names)

  expected <- list(
    brma.mv_v14_konstantopoulos2011_cs        = "district",
    brma.mv_v14_assink2016_nested             = c("study/esid", "esid_study", "study"),
    brma.mv_v14_ishak2007_har                 = "study",
    brma.mv_v14_begg1989_study_treatment      = c("study", "treatment")
  )

  for (name in fit_names) {
    fit_brma <- fits[[name]]
    out_all  <- summary_heterogeneity(fit_brma, component = "all")

    if (length(expected[[name]]) == 1L) {
      expect_equal(names(out_all), "estimates", info = name)
    } else {
      expect_equal(names(out_all), expected[[name]], info = name)
      total <- summary_heterogeneity(fit_brma, component = "total")
      expect_equal(names(total), "estimates", info = paste(name, "total"))
    }

    for (component in expected[[name]]) {
      selected <- summary_heterogeneity(fit_brma, component = component)
      expect_equal(names(selected), "estimates",
                   info = paste(name, component))
    }
    expect_error(
      summary_heterogeneity(fit_brma, component = "missing"),
      "Unknown|Available",
      info = name
    )
  }
})

test_that("v14 brma.mv random-covariance parameters match metafor references", {

  fit_names <- .mv_active_metafor_fit_names(c(
    "brma.mv_v14_konstantopoulos2011_cs",
    "brma.mv_v14_ishak2007_har",
    "brma.mv_v14_begg1989_study_treatment"
  ))
  skip_if_missing_fits(fit_names)

  cases <- list(
    list(
      name      = "brma.mv_v14_konstantopoulos2011_cs",
      row       = "cor",
      expected  = function(m) m[["rho"]],
      tolerance = 0.12
    ),
    list(
      name      = "brma.mv_v14_ishak2007_har",
      row       = "cor",
      expected  = function(m) m[["rho"]],
      tolerance = 0.12
    ),
    list(
      name      = "brma.mv_v14_begg1989_study_treatment",
      row       = "study: cor",
      expected  = function(m) m[["rho"]],
      tolerance = 1e-12
    )
  )
  cases <- Filter(function(case) case[["name"]] %in% fit_names, cases)

  for (case in cases) {
    name     <- case[["name"]]
    observed <- .mv_summary_mean(name, "estimates_random", case[["row"]])
    expected <- case[["expected"]](.mv_metafor(name))

    .expect_close_abs(
      observed  = observed,
      expected  = expected,
      tolerance = case[["tolerance"]],
      label     = paste(name, case[["row"]])
    )
  }
})

test_that("v14 brma.mv fixed fitted values and residuals match metafor references", {

  fit_names <- .mv_active_metafor_fit_names()
  skip_if_missing_fits(fit_names)

  tolerances <- c(
    brma.mv_v14_konstantopoulos2011_cs        = 0.05,
    brma.mv_v14_assink2016_nested             = 0.05,
    brma.mv_v14_ishak2007_har                 = 0.40,
    brma.mv_v14_begg1989_study_treatment      = 0.05
  )

  for (name in fit_names) {
    fit_brma    <- fits[[name]]
    fit_metafor <- .mv_metafor(name)
    tolerance   <- tolerances[[name]]

    reference_beta <- if (identical(name, "brma.mv_v14_ishak2007_har")) {
      .ishak_prior_aware_beta(fit_metafor)
    } else {
      fit_metafor[["beta"]]
    }
    brma_pred    <- .brma_fixed_prediction_mean(fit_brma)
    metafor_pred <- .metafor_fixed_prediction(fit_metafor, reference_beta)
    brma_resid   <- residuals(
      fit_brma,
      type               = "outcome",
      conditioning_depth = "marginal"
    )
    metafor_resid <- .metafor_marginal_residual(
      fit_metafor,
      reference_beta
    )

    expect_equal(length(brma_pred), length(metafor_pred))
    expect_equal(length(brma_resid), length(metafor_resid))

    .expect_close_abs(
      observed  = max(abs(brma_pred - metafor_pred)),
      expected  = 0,
      tolerance = tolerance,
      label     = paste(name, "maximum fixed fitted-value difference")
    )
    .expect_close_abs(
      observed  = max(abs(brma_resid - metafor_resid)),
      expected  = 0,
      tolerance = tolerance,
      label     = paste(name, "maximum marginal residual difference")
    )
  }
})


test_that("explicit brma.mv new-effect predictions track metafor targets", {

  fit_names <- .mv_active_metafor_fit_names(c(
    "brma.mv_v14_konstantopoulos2011_cs",
    "brma.mv_v14_assink2016_nested",
    "brma.mv_v14_ishak2007_har"
  ))
  skip_if_missing_fits(fit_names)

  predict_rma <- getS3method(
    "predict",
    "rma",
    envir = asNamespace("metafor")
  )
  cases <- list(
    list(
      name = "brma.mv_v14_konstantopoulos2011_cs",
      newdata = function(dat) dat[1L, c("school", "district"), drop = FALSE],
      reference = function(fit, dat) {
        predict_rma(fit, tau2.levels = dat[["school"]][[1L]])
      },
      tolerance = c(mean = .10, lower = .15, upper = .15)
    ),
    list(
      name = "brma.mv_v14_assink2016_nested",
      newdata = function(dat) data.frame(.prediction_row = 1L),
      reference = function(fit, dat) predict_rma(fit),
      tolerance = c(mean = .10, lower = .20, upper = .20)
    ),
    list(
      name = "brma.mv_v14_ishak2007_har",
      newdata = function(dat) {
        dat[1L, c("time", "time_factor", "study"), drop = FALSE]
      },
      reference = function(fit, dat) {
        predict_rma(
          fit,
          newmods     = fit[["X"]][1L, , drop = FALSE],
          tau2.levels = dat[["time"]][[1L]]
        )
      },
      tolerance = c(mean = .75, lower = 1.75, upper = 1.75)
    )
  )
  cases <- Filter(function(case) case[["name"]] %in% fit_names, cases)

  for (case in cases) {
    name     <- case[["name"]]
    fit_brma <- fits[[name]]
    fit_meta <- .mv_metafor(name)
    dat      <- info[[name]][["data"]]
    newdata  <- case[["newdata"]](dat)
    reference <- case[["reference"]](fit_meta, newdata)

    set.seed(81)
    observed <- summary(predict(
      fit_brma,
      newdata = newdata,
      type    = "estimate",
      quiet   = TRUE
    ))

    .expect_close_abs(
      observed  = observed[1L, "Mean"],
      expected  = reference[["pred"]][[1L]],
      tolerance = case[["tolerance"]][["mean"]],
      label     = paste(name, "new-effect center")
    )
    .expect_close_abs(
      observed  = observed[1L, "0.025"],
      expected  = reference[["pi.lb"]][[1L]],
      tolerance = case[["tolerance"]][["lower"]],
      label     = paste(name, "new-effect lower bound")
    )
    .expect_close_abs(
      observed  = observed[1L, "0.975"],
      expected  = reference[["pi.ub"]][[1L]],
      tolerance = case[["tolerance"]][["upper"]],
      label     = paste(name, "new-effect upper bound")
    )
  }
})

test_that("v14 brma.mv ranef components track metafor references", {

  fit_names <- .mv_active_metafor_fit_names()
  skip_if_missing_fits(fit_names)

  cases <- list(
    list(
      name              = "brma.mv_v14_konstantopoulos2011_cs",
      component         = "district",
      metafor_component = "~school | district",
      tolerance         = 0.08
    ),
    list(
      name              = "brma.mv_v14_assink2016_nested",
      component         = "study",
      metafor_component = "study",
      tolerance         = 0.08
    ),
    list(
      name              = "brma.mv_v14_assink2016_nested",
      component         = "esid_study",
      metafor_component = "study/esid",
      tolerance         = 0.08
    ),
    list(
      name              = "brma.mv_v14_ishak2007_har",
      component         = "study",
      metafor_component = "~time | study",
      tolerance         = 2.25,
      rank              = 0.98
    ),
    list(
      name              = "brma.mv_v14_begg1989_study_treatment",
      component         = "study",
      metafor_component = "study",
      tolerance         = 0.08
    ),
    list(
      name              = "brma.mv_v14_begg1989_study_treatment",
      component         = "treatment",
      metafor_component = "~trt | study",
      tolerance         = 0.08
    )
  )
  cases <- Filter(function(case) case[["name"]] %in% fit_names, cases)

  for (case in cases) {
    name        <- case[["name"]]
    fit_brma    <- fits[[name]]
    fit_metafor <- .mv_metafor(name)
    posterior   <- .mv_posterior_subset(fit_brma)

    observed <- .mv_sample_means(ranef(
      fit_brma,
      component          = case[["component"]],
      .posterior_samples = posterior
    ))
    expected <- metafor::ranef(
      fit_metafor
    )[[case[["metafor_component"]]]][["intrcpt"]]

    expect_equal(length(observed), length(expected), info = name)
    expect_true(all(is.finite(observed)), info = name)
    expect_true(all(is.finite(expected)), info = name)
    .expect_close_abs(
      observed  = max(abs(observed - expected)),
      expected  = 0,
      tolerance = case[["tolerance"]],
      label     = paste(name, case[["component"]], "ranef maximum difference")
    )
    if (!is.null(case[["rank"]])) {
      expect_true(
        stats::cor(observed, expected, method = "spearman") > case[["rank"]],
        info = paste(name, case[["component"]], "ranef rank correlation")
      )
    }
  }
})

test_that("v14 brma.mv ranef, blup, and true_effects use consistent targets", {

  fit_names <- .mv_active_metafor_fit_names()
  skip_if_missing_fits(fit_names)

  expected_components <- list(
    brma.mv_v14_konstantopoulos2011_cs        = "district",
    brma.mv_v14_assink2016_nested             = c("esid_study", "study"),
    brma.mv_v14_ishak2007_har                 = "study",
    brma.mv_v14_begg1989_study_treatment      = c("study", "treatment")
  )

  for (name in fit_names) {
    fit_brma  <- fits[[name]]
    posterior <- .mv_posterior_subset(fit_brma, n = 300)

    terms    <- as.matrix(predict(
      fit_brma,
      type               = "terms",
      quiet              = TRUE,
      .posterior_samples = posterior
    ))
    estimate <- as.matrix(predict(
      fit_brma,
      type               = "blup",
      quiet              = TRUE,
      .posterior_samples = posterior
    ))
    total <- as.matrix(ranef(
      fit_brma,
      component          = "total",
      expand             = TRUE,
      .posterior_samples = posterior
    ))
    components <- ranef(
      fit_brma,
      simplify           = FALSE,
      expand             = TRUE,
      .posterior_samples = posterior
    )

    expect_equal(names(components), expected_components[[name]], info = name)
    expect_equal(
      unname(total),
      unname(Reduce(`+`, lapply(components, as.matrix))),
      tolerance = 1e-12,
      info      = paste(name, "total random effect")
    )
    expect_equal(
      unname(total),
      unname(estimate - terms),
      tolerance = 1e-12,
      info      = paste(name, "estimate minus terms")
    )
    expect_equal(
      unname(as.matrix(blup(fit_brma, .posterior_samples = posterior))),
      unname(as.matrix(true_effects(fit_brma, .posterior_samples = posterior))),
      tolerance = 1e-12,
      info      = paste(name, "blup true_effects alias")
    )
    expect_error(
      ranef(fit_brma, component = "missing",
            .posterior_samples = posterior),
      "Unknown random-effect component"
    )
  }
})

test_that("v14 brma.mv hatvalues track metafor leverages", {

  fit_names <- .mv_active_metafor_fit_names()
  skip_if_missing_fits(fit_names)

  for (name in fit_names) {
    fit_brma    <- fits[[name]]
    fit_metafor <- .mv_metafor(name)

    brma_hat    <- suppressWarnings(hatvalues(fit_brma, max_samples = 200))
    metafor_hat <- as.numeric(metafor::hatvalues.rma.mv(fit_metafor))

    expect_equal(length(brma_hat), length(metafor_hat))
    expect_true(all(is.finite(brma_hat)))
    expect_true(all(is.finite(metafor_hat)))
    expect_true(
      stats::cor(brma_hat, metafor_hat, method = "spearman") > 0.85,
      info = paste(name, "hatvalue rank correlation")
    )
    .expect_close_abs(
      observed  = max(abs(brma_hat - metafor_hat)),
      expected  = 0,
      tolerance = 0.15,
      label     = paste(name, "maximum hatvalue difference")
    )
  }
})
