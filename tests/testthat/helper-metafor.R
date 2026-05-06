.test_object <- function(name) {

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (exists(name, envir = frames[[i]], inherits = FALSE)) {
      return(get(name, envir = frames[[i]], inherits = FALSE))
    }
  }

  return(get(name, envir = parent.frame(), inherits = TRUE))
}

.test_fits <- function() {

  return(.test_object("fits"))
}

.test_info <- function() {

  return(.test_object("info"))
}

.case_fit <- function(case) {

  return(.test_fits()[[case_name(case)]])
}

.case_metafor <- function(case) {

  return(.test_info()[[case_name(case)]][["metafor"]])
}

.metafor_dfbetas <- function(fit_metafor) {

  if (inherits(fit_metafor, "rma.mv")) {
    return(as.data.frame(metafor::dfbetas.rma.mv(fit_metafor)))
  }

  return(as.data.frame(metafor::dfbetas.rma.uni(fit_metafor)))
}

.expect_residuals_equal <- function(name, actual, expected, tolerance, label) {

  testthat::expect_equal(
    actual,
    as.vector(expected),
    tolerance = tolerance,
    info      = paste(name, label, "matches metafor")
  )
}

.expect_residual_table_equal <- function(name, actual, expected, tolerance,
                                         cols = c("resid", "se", "z"),
                                         rows = NULL, z_tolerance = tolerance,
                                         label = "rstandard") {

  if (is.null(rows)) {
    rows <- seq_len(nrow(actual))
  }

  for (col in cols) {
    col_tolerance <- if (col == "z") z_tolerance else tolerance
    testthat::expect_equal(
      actual[[col]][rows],
      expected[[col]][rows],
      tolerance = col_tolerance,
      info      = paste(name, label, col, "matches metafor")
    )
  }
}

.expect_residual_rstudent_rank <- function(name, brma_rstudent, metafor_rstudent,
                                           resid_cor = 0.90, se_cor = NULL,
                                           z_cor = 0.90) {

  testthat::expect_true(all(is.finite(brma_rstudent$resid)),
                        info = paste(name, "LOO-PIT residual means finite"))
  testthat::expect_true(all(is.finite(brma_rstudent$se)) && all(brma_rstudent$se > 0),
                        info = paste(name, "LOO predictive SEs finite and positive"))
  testthat::expect_true(
    stats::cor(brma_rstudent$resid, metafor_rstudent$resid,
               method = "spearman", use = "complete.obs") > resid_cor,
    info = paste(name, "LOO-PIT residual ranks align with metafor")
  )
  if (!is.null(se_cor)) {
    testthat::expect_true(
      stats::cor(brma_rstudent$se, metafor_rstudent$se,
                 method = "spearman", use = "complete.obs") > se_cor,
      info = paste(name, "LOO predictive SE ranks align with metafor")
    )
  }
  testthat::expect_true(
    stats::cor(brma_rstudent$z, metafor_rstudent$z,
               method = "spearman", use = "complete.obs") > z_cor,
    info = paste(name, "LOO-PIT z ranks align with metafor")
  )
}

expect_residuals_match_metafor <- function(case) {

  name        <- case_name(case)
  fit_brma    <- .case_fit(case)
  fit_metafor <- .case_metafor(case)
  kind        <- case_value(case, "kind", "standard")
  tolerance   <- case_value(case, "tolerance", 0.05)

  if (kind %in% c("standard", "regression", "multilevel", "multilevel_no_loo",
                  "pet", "pet_reg", "peese", "peese_reg")) {
    .expect_residuals_equal(name, residuals(fit_brma), residuals(fit_metafor),
                            tolerance, "outcome residuals")
    .expect_residuals_equal(name, residuals(fit_brma, type = "pearson"),
                            residuals(fit_metafor, type = "pearson"),
                            tolerance, "pearson residuals")
    .expect_residual_table_equal(name, rstandard(fit_brma), rstandard(fit_metafor),
                                 tolerance)
  }

  if (kind == "regression") {
    .expect_residual_table_equal(
      name,
      rstandard(fit_brma, conditioning_depth = "estimate"),
      rstandard(fit_metafor, type = "conditional"),
      tolerance,
      label = "conditional rstandard"
    )
  }

  if (kind == "interaction") {
    .expect_residual_table_equal(
      name,
      rstandard(fit_brma),
      rstandard(fit_metafor),
      tolerance,
      cols = c("resid", "se")
    )

    metafor_conditional <- suppressWarnings(rstandard(fit_metafor, type = "conditional"))
    brma_conditional    <- rstandard(fit_brma, conditioning_depth = "estimate")
    valid <- is.finite(metafor_conditional$se) & is.finite(metafor_conditional$z)

    .expect_residual_table_equal(
      name,
      brma_conditional,
      metafor_conditional,
      tolerance,
      cols = c("resid", "se"),
      rows = which(valid),
      label = "conditional rstandard"
    )
    testthat::expect_true(all(abs(brma_conditional$resid[!valid]) < 1e-8),
                          info = paste(name, "saturated residuals are zero"))
    testthat::expect_true(all(brma_conditional$se[!valid] < 1e-8),
                          info = paste(name, "saturated SEs collapse"))
    testthat::expect_true(all(abs(brma_conditional$z[!valid]) < 1e-8),
                          info = paste(name, "saturated z values are zero"))
  }

  if (kind %in% c("multilevel", "multilevel_no_loo")) {
    metafor_conditional <- metafor_rstandard_conditional_mv(fit_metafor)
    brma_conditional    <- rstandard(fit_brma, conditioning_depth = "estimate")
    brma_theta          <- colMeans(as.matrix(predict(fit_brma, type = "estimate", quiet = TRUE)))
    metafor_theta       <- fit_metafor[["yi"]] - metafor_conditional$resid

    testthat::expect_equal(as.vector(brma_theta), as.vector(metafor_theta),
                           tolerance = 0.11,
                           info = paste(name, "estimate-level BLUP oracle"))
    .expect_residual_table_equal(
      name,
      brma_conditional,
      metafor_conditional,
      tolerance = 0.11,
      cols = "resid",
      label = "multilevel conditional rstandard"
    )
    .expect_residual_table_equal(
      name,
      brma_conditional,
      metafor_conditional,
      tolerance = 0.05,
      cols = "se",
      label = "multilevel conditional rstandard"
    )
  }

  if (kind %in% c("selection_pos", "selection_reg", "selection_neg")) {
    .expect_residuals_equal(name, residuals(fit_brma), residuals(fit_metafor),
                            tolerance, "outcome residuals")
    testthat::expect_error(
      residuals(fit_brma, type = "pearson"),
      "Pearson residuals are not available for selection models",
      info = paste(name, "pearson residuals are rejected")
    )
    testthat::expect_error(
      residuals(fit_brma, type = "rstandard"),
      "Standardized residuals.*are not available for selection models",
      info = paste(name, "rstandard residuals are rejected")
    )
  }

  if (kind %in% c("glmm", "glmm_reg")) {
    .expect_residuals_equal(name, residuals(fit_brma), residuals(fit_metafor),
                            tolerance, "outcome residuals")
    testthat::expect_error(
      residuals(fit_brma, type = "pearson"),
      "Pearson residuals are only available for normal outcome models",
      info = paste(name, "pearson residuals are rejected")
    )
    testthat::expect_error(
      residuals(fit_brma, type = "rstandard"),
      "Standardized residuals.*are only available for normal outcome models",
      info = paste(name, "rstandard residuals are rejected")
    )
  }

  rstudent_kind <- case_value(case, "rstudent", NULL)
  if (identical(rstudent_kind, "equal")) {
    .expect_residual_table_equal(
      name,
      suppressWarnings(rstudent(fit_brma)),
      rstudent(fit_metafor),
      tolerance,
      z_tolerance = if (kind == "pet") 0.10 else tolerance,
      label = "rstudent"
    )
  } else if (identical(rstudent_kind, "rank")) {
    .expect_residual_rstudent_rank(
      name,
      suppressWarnings(rstudent(fit_brma)),
      rstudent(fit_metafor),
      resid_cor = if (kind == "multilevel") 0.65 else 0.90,
      se_cor    = if (kind == "multilevel") 0.80 else NULL,
      z_cor     = if (kind == "multilevel") 0.65 else 0.90
    )
  } else if (identical(rstudent_kind, "selection_pos") ||
             identical(rstudent_kind, "selection_neg")) {
    brma_resid    <- residuals(fit_brma)
    brma_rstudent <- suppressWarnings(rstudent(fit_brma))
    if (rstudent_kind == "selection_pos") {
      testthat::expect_true(brma_rstudent$z[4] > 4,
                            info = paste(name, "large positive outlier"))
    } else {
      testthat::expect_true(brma_rstudent$z[4] < -4,
                            info = paste(name, "large negative outlier"))
    }
    testthat::expect_true(all(abs(brma_rstudent$z[-4]) < 4),
                          info = paste(name, "single large outlier"))
    testthat::expect_true(stats::cor(brma_resid, brma_rstudent$z) > 0.9,
                          info = paste(name, "rstudent and residuals align"))
    testthat::expect_true(all(is.finite(brma_rstudent$se)),
                          info = paste(name, "selection LOO SEs finite"))
    testthat::expect_true(all(brma_rstudent$se > 0),
                          info = paste(name, "selection LOO SEs positive"))
  } else if (identical(rstudent_kind, "glmm_align")) {
    brma_resid    <- residuals(fit_brma)
    brma_rstudent <- suppressWarnings(rstudent(fit_brma, type = "estimate"))
    testthat::expect_true(stats::cor(brma_resid, brma_rstudent$z) > 0.9,
                          info = paste(name, "rstudent and residuals align"))
  }
}

.sample_means <- function(x) {

  return(summary(x)[, "Mean"])
}

.sample_mean <- function(x, row, column = "Mean") {

  return(summary(x)[row, column])
}

.expect_blup_tracks_observed <- function(name, fit_brma, theta) {

  testthat::expect_true(
    stats::cor(fit_brma$data$outcome$yi, theta, method = "spearman") > 0.8,
    info = paste(name, "BLUPs correlate with observed estimates")
  )
}

.expect_blup_shrinks <- function(name, fit_brma, theta, mu) {

  yi <- fit_brma$data$outcome$yi
  testthat::expect_true(
    mean(abs(theta - mu)) < mean(abs(yi - mu)),
    info = paste(name, "BLUPs shrink toward pooled estimate")
  )
  .expect_blup_tracks_observed(name, fit_brma, theta)
}

expect_prediction_matches_metafor <- function(case) {

  name        <- case_name(case)
  fit_brma    <- .case_fit(case)
  fit_metafor <- .case_metafor(case)
  kind        <- case_value(case, "kind", "simple")
  tolerance   <- case_value(case, "tolerance", 0.05)
  tau_tol     <- case_value(case, "tau_tolerance", tolerance)
  fitted_tol  <- case_value(case, "fitted_tolerance", tolerance)

  if (kind == "simple") {
    mu_predict <- .sample_mean(predict(fit_brma, type = "terms", newdata = TRUE), "mu")
    mu_wrapper <- .sample_mean(pooled_effect(fit_brma), "mu")
    testthat::expect_equal(mu_predict, mu_wrapper,
                           info = paste(name, "pooled_effect wrapper"))
    testthat::expect_equal(mu_predict, fit_metafor$beta[[1]], tolerance = tolerance,
                           info = paste(name, "pooled effect matches metafor"))

    tau_predict <- .sample_mean(predict(fit_brma, type = "terms.scale", newdata = TRUE), "tau")
    tau_wrapper <- .sample_mean(pooled_heterogeneity(fit_brma), "tau")
    testthat::expect_equal(tau_predict, tau_wrapper,
                           info = paste(name, "pooled_heterogeneity wrapper"))
    testthat::expect_equal(tau_predict, sqrt(fit_metafor$tau2), tolerance = tau_tol,
                           info = paste(name, "tau matches metafor"))

    theta_predict <- .sample_means(predict(fit_brma, type = "effect"))
    theta_blup    <- .sample_means(blup(fit_brma))
    theta_true    <- .sample_means(true_effects(fit_brma))
    testthat::expect_equal(theta_predict, theta_blup,
                           info = paste(name, "predict effect matches blup"))
    testthat::expect_equal(theta_blup, theta_true,
                           info = paste(name, "true_effects match blup"))
    testthat::expect_equal(theta_blup, metafor::blup(fit_metafor)$pred,
                           tolerance = 0.05,
                           info = paste(name, "BLUPs match metafor"))
    .expect_blup_tracks_observed(name, fit_brma, theta_true)
    return(invisible(TRUE))
  }

  if (kind %in% c("regression", "interaction")) {
    metafor_fitted <- predict(fit_metafor)$pred
    brma_fitted    <- .sample_means(predict(fit_brma, type = "terms"))

    if (kind == "interaction") {
      testthat::expect_true(stats::cor(brma_fitted, metafor_fitted) > 0.90,
                            info = paste(name, "fitted values align"))
    } else {
      testthat::expect_equal(brma_fitted, metafor_fitted, tolerance = fitted_tol,
                             info = paste(name, "fitted values match metafor"))
    }

    testthat::expect_equal(
      .sample_mean(pooled_effect(fit_brma), "mu"),
      mean(metafor_fitted),
      tolerance = tolerance,
      info      = paste(name, "pooled effect matches mean metafor fitted")
    )

    theta_brma    <- .sample_means(blup(fit_brma))
    theta_metafor <- metafor::blup(fit_metafor)$pred
    if (kind == "interaction") {
      testthat::expect_true(stats::cor(theta_brma, theta_metafor) > 0.90,
                            info = paste(name, "BLUPs align with metafor"))
    } else {
      testthat::expect_equal(theta_brma, theta_metafor, tolerance = 0.05,
                             info = paste(name, "BLUPs match metafor"))
    }
    .expect_blup_tracks_observed(name, fit_brma, theta_brma)
    return(invisible(TRUE))
  }

  if (kind == "scale") {
    metafor_fitted <- predict(fit_metafor)$pred
    brma_fitted    <- .sample_means(predict(fit_brma, type = "terms"))
    testthat::expect_equal(brma_fitted, metafor_fitted, tolerance = fitted_tol,
                           info = paste(name, "location predictions"))

    metafor_tau <- sqrt(as.vector(exp(fit_metafor$Z %*% fit_metafor$alpha)))
    brma_tau    <- .sample_means(predict(fit_brma, type = "terms.scale"))
    testthat::expect_equal(brma_tau, metafor_tau, tolerance = tau_tol,
                           info = paste(name, "study tau predictions"))
    testthat::expect_equal(.sample_mean(pooled_heterogeneity(fit_brma), "tau"),
                           mean(metafor_tau), tolerance = tau_tol,
                           info = paste(name, "pooled tau"))

    theta_brma <- .sample_means(blup(fit_brma))
    testthat::expect_equal(theta_brma, metafor::blup(fit_metafor)$pred,
                           tolerance = 0.10,
                           info = paste(name, "BLUPs match metafor"))
    .expect_blup_tracks_observed(name, fit_brma, theta_brma)
    return(invisible(TRUE))
  }

  if (kind == "multilevel") {
    brma_mu  <- .sample_mean(pooled_effect(fit_brma), "mu")
    brma_tau <- .sample_mean(pooled_heterogeneity(fit_brma), "tau")

    testthat::expect_equal(brma_mu, fit_metafor$beta[[1]], tolerance = tolerance,
                           info = paste(name, "pooled effect"))
    testthat::expect_equal(brma_tau, sqrt(fit_metafor$tau2), tolerance = tau_tol,
                           info = paste(name, "total tau"))

    theta <- .sample_means(blup(fit_brma))
    testthat::expect_equal(length(theta), nrow(fit_brma$data$outcome),
                           info = paste(name, "BLUP length"))
    .expect_blup_shrinks(name, fit_brma, theta, brma_mu)
    return(invisible(TRUE))
  }

  if (kind %in% c("glmm", "glmm_reg")) {
    brma_mu  <- .sample_mean(pooled_effect(fit_brma), "mu")
    brma_tau <- .sample_mean(pooled_heterogeneity(fit_brma), "tau")
    testthat::expect_equal(brma_mu, fit_metafor$beta[[1]], tolerance = tolerance,
                           info = paste(name, "pooled effect"))
    testthat::expect_equal(brma_tau, sqrt(fit_metafor$tau2), tolerance = tau_tol,
                           info = paste(name, "pooled tau"))

    theta <- .sample_means(blup(fit_brma))
    testthat::expect_equal(length(theta), nrow(fit_brma$data$outcome),
                           info = paste(name, "BLUP length"))
    return(invisible(TRUE))
  }

  if (kind == "selection") {
    brma_mu  <- .sample_mean(pooled_effect(fit_brma), "mu")
    brma_tau <- .sample_mean(pooled_heterogeneity(fit_brma), "tau")
    testthat::expect_equal(brma_mu, fit_metafor$beta[[1]], tolerance = tolerance,
                           info = paste(name, "pooled effect"))
    testthat::expect_equal(brma_tau, sqrt(fit_metafor$tau2), tolerance = tau_tol,
                           info = paste(name, "pooled tau"))

    theta <- .sample_means(blup(fit_brma))
    testthat::expect_equal(length(theta), nrow(fit_brma$data$outcome),
                           info = paste(name, "BLUP length"))
    .expect_blup_shrinks(name, fit_brma, theta, brma_mu)
    return(invisible(TRUE))
  }

  if (kind %in% c("pet", "pet_reg", "peese", "peese_reg")) {
    is_regression_bias <- kind %in% c("pet_reg", "peese_reg")
    bias_column        <- if (kind %in% c("pet", "pet_reg")) "sqrt(vi)" else "vi"

    if (is_regression_bias) {
      bias_newmods <- matrix(
        0,
        nrow     = nrow(fit_metafor$X),
        ncol     = 1,
        dimnames = list(NULL, bias_column)
      )
      metafor_newmods_bias_adjusted <- cbind(
        bias_newmods,
        fit_metafor$X[, -c(1, 2), drop = FALSE]
      )
      metafor_mu <- mean(predict(fit_metafor, newmods = metafor_newmods_bias_adjusted)$pred)
    } else {
      metafor_newmods_bias_adjusted <- rep(0, nrow(fit_brma$data$outcome))
      metafor_mu <- predict(fit_metafor, newmods = 0)$pred
    }

    brma_mu <- .sample_mean(pooled_effect(fit_brma), "mu")
    testthat::expect_equal(brma_mu, metafor_mu, tolerance = tolerance,
                           info = paste(name, "bias-adjusted pooled effect"))
    testthat::expect_equal(.sample_mean(pooled_heterogeneity(fit_brma), "tau"),
                           sqrt(fit_metafor$tau2), tolerance = tau_tol,
                           info = paste(name, "pooled tau"))

    metafor_fitted <- predict(fit_metafor)$pred
    brma_fitted    <- .sample_means(predict(fit_brma, type = "terms"))
    testthat::expect_equal(brma_fitted, metafor_fitted, tolerance = fitted_tol,
                           info = paste(name, "unadjusted fitted values"))

    theta <- .sample_means(blup(fit_brma))
    testthat::expect_equal(theta, metafor::blup(fit_metafor)$pred,
                           tolerance = 0.05,
                           info = paste(name, "BLUPs match metafor"))
    .expect_blup_tracks_observed(name, fit_brma, theta)

    adjusted_metafor <- predict(fit_metafor, newmods = metafor_newmods_bias_adjusted)$pred
    adjusted_brma    <- .sample_means(predict(fit_brma, type = "terms", bias_adjusted = TRUE))
    testthat::expect_equal(adjusted_brma, adjusted_metafor, tolerance = 0.05,
                           info = paste(name, "bias-adjusted fitted values"))
    return(invisible(TRUE))
  }

  stop(paste0("Unknown prediction case kind '", kind, "'."), call. = FALSE)
}

.metafor_linear_predictor <- function(fit_metafor, formula, newdata) {

  x_new <- stats::model.matrix(formula, newdata)
  return(as.vector(x_new %*% fit_metafor[["beta"]]))
}

.newdata_factor <- function(fit_brma, variable, values = NULL) {

  levels_variable <- levels(fit_brma[["data"]][["mods"]][[variable]])
  if (is.null(values)) {
    values <- levels_variable
  }

  return(factor(values, levels = levels_variable))
}

.newdata_prediction_case_data <- function(kind, fit_brma) {

  if (kind == "normal_continuous") {
    return(data.frame(
      ablat = c(10, 40, 70),
      year  = c(1950, 1968, 1980)
    ))
  }

  if (kind == "normal_factor") {
    alloc <- .newdata_factor(
      fit_brma,
      "alloc",
      values = c("random", "alternate", "systematic")
    )
    return(data.frame(alloc = alloc))
  }

  if (kind == "location_scale") {
    meta_levels <- levels(fit_brma[["data"]][["mods"]][["meta"]])
    return(data.frame(
      meta  = factor(c(meta_levels[1], meta_levels[2], meta_levels[1]),
                     levels = meta_levels),
      ni100 = c(0.5, 1.0, 2.0)
    ))
  }

  if (kind == "multilevel_mod") {
    return(data.frame(vi = c(0.02, 0.08, 0.20)))
  }

  if (kind == "glmm_factor") {
    alloc <- .newdata_factor(fit_brma, "alloc")
    return(data.frame(alloc = alloc))
  }

  if (kind == "pet_reg") {
    preregistered <- .newdata_factor(fit_brma, "Preregistered")
    return(data.frame(
      Preregistered = preregistered,
      sei           = c(0.10, 0.20),
      vi            = c(0.01, 0.04)
    ))
  }

  if (kind == "peese_reg") {
    preregistered <- .newdata_factor(fit_brma, "Preregistered")
    return(data.frame(
      Preregistered = preregistered,
      vi            = c(0.01, 0.04)
    ))
  }

  if (kind == "selection_reg") {
    preregistered <- .newdata_factor(fit_brma, "Preregistered")
    return(data.frame(Preregistered = preregistered))
  }

  stop(paste0("Unknown newdata prediction case kind '", kind, "'."), call. = FALSE)
}

.metafor_newdata_location <- function(kind, fit_metafor, newdata,
                                      bias_adjusted = FALSE) {

  if (kind == "normal_continuous") {
    return(.metafor_linear_predictor(fit_metafor, ~ ablat + year, newdata))
  }

  if (kind %in% c("normal_factor", "glmm_factor")) {
    return(.metafor_linear_predictor(fit_metafor, ~ alloc, newdata))
  }

  if (kind == "location_scale") {
    return(.metafor_linear_predictor(fit_metafor, ~ meta + ni100, newdata))
  }

  if (kind == "multilevel_mod") {
    return(.metafor_linear_predictor(fit_metafor, ~ vi, newdata))
  }

  if (kind %in% c("pet_reg", "peese_reg")) {
    metafor_newdata <- newdata
    if (bias_adjusted) {
      metafor_newdata[["vi"]] <- 0
    }
    formula <- if (kind == "pet_reg") {
      ~ sqrt(vi) + Preregistered
    } else {
      ~ vi + Preregistered
    }
    return(.metafor_linear_predictor(
      fit_metafor,
      formula,
      metafor_newdata
    ))
  }

  if (kind == "selection_reg") {
    return(.metafor_linear_predictor(fit_metafor, ~ Preregistered, newdata))
  }

  stop(paste0("Unknown newdata prediction case kind '", kind, "'."), call. = FALSE)
}

.metafor_newdata_tau <- function(kind, fit_metafor, newdata) {

  if (kind == "location_scale") {
    z_new <- stats::model.matrix(~ ni100, newdata)
    return(sqrt(as.vector(exp(z_new %*% fit_metafor[["alpha"]]))))
  }

  return(NULL)
}

expect_newdata_prediction_matches_metafor <- function(case) {

  name        <- case_name(case)
  fit_brma    <- .case_fit(case)
  fit_metafor <- .case_metafor(case)
  kind        <- case_value(case, "kind")
  tolerance   <- case_value(case, "tolerance", 0.05)
  tau_tol     <- case_value(case, "tau_tolerance", tolerance)
  newdata     <- .newdata_prediction_case_data(kind, fit_brma)

  brma_terms <- .sample_means(predict(
    fit_brma,
    newdata = newdata,
    type    = "terms",
    quiet   = TRUE
  ))
  metafor_terms <- .metafor_newdata_location(kind, fit_metafor, newdata)

  testthat::expect_equal(
    brma_terms,
    metafor_terms,
    tolerance = tolerance,
    info      = paste(name, "newdata location predictions match metafor")
  )

  if (kind == "pet_reg") {
    brma_adjusted <- .sample_means(predict(
      fit_brma,
      newdata       = newdata,
      type          = "terms",
      bias_adjusted = TRUE,
      quiet         = TRUE
    ))
    metafor_adjusted <- .metafor_newdata_location(
      kind          = kind,
      fit_metafor   = fit_metafor,
      newdata       = newdata,
      bias_adjusted = TRUE
    )
    testthat::expect_equal(
      brma_adjusted,
      metafor_adjusted,
      tolerance = tolerance,
      info      = paste(name, "bias-adjusted newdata predictions match metafor")
    )
  }

  metafor_tau <- .metafor_newdata_tau(kind, fit_metafor, newdata)
  if (!is.null(metafor_tau)) {
    brma_tau <- .sample_means(predict(
      fit_brma,
      newdata = newdata,
      type    = "terms.scale",
      quiet   = TRUE
    ))
    testthat::expect_equal(
      brma_tau,
      metafor_tau,
      tolerance = tau_tol,
      info      = paste(name, "newdata scale predictions match metafor")
    )
  }
}

expect_hatvalues_match_metafor <- function(case) {

  name        <- case_name(case)
  fit_brma    <- .case_fit(case)
  fit_metafor <- .case_metafor(case)
  tolerance   <- case_value(case, "tolerance", 0.05)

  testthat::expect_equal(
    unname(hatvalues(fit_brma)),
    as.vector(hatvalues(fit_metafor)),
    tolerance = tolerance,
    info      = paste(name, "hatvalues match metafor")
  )
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

expect_vif_matches_metafor <- function(case) {

  name        <- case_name(case)
  fit_brma    <- .case_fit(case)
  fit_metafor <- .case_metafor(case)
  tolerance   <- case_value(case, "tolerance", 0.10)
  btt         <- case_value(case, "btt", NULL)

  brma_vif    <- vif(fit_brma, posterior_correlation = FALSE)[["vif"]]
  metafor_vif <- metafor_vif_table(fit_brma, fit_metafor, btt = btt)

  expect_vif_table(brma_vif, nrow(metafor_vif), info = name)
  testthat::expect_equal(
    brma_vif[["df"]],
    metafor_vif[["df"]],
    info = paste(name, "VIF degrees of freedom match metafor")
  )
  testthat::expect_equal(
    brma_vif[["GVIF"]],
    metafor_vif[["GVIF"]],
    tolerance = tolerance,
    info      = paste(name, "GVIF matches metafor")
  )
  testthat::expect_equal(
    brma_vif[["GVIF^(1/(2*df))"]],
    metafor_vif[["GVIF^(1/(2*df))"]],
    tolerance = tolerance,
    info      = paste(name, "adjusted GVIF matches metafor")
  )
}

expect_dfbetas_match_metafor <- function(case) {

  name        <- case_name(case)
  fit_brma    <- .case_fit(case)
  fit_metafor <- .case_metafor(case)
  oracle      <- case_value(case, "oracle", "equal")
  tolerance   <- case_value(case, "tolerance", 0.10)

  brma_dfbetas <- dfbetas(fit_brma)
  if (oracle == "structure") {
    expect_dfbetas_table(brma_dfbetas, nobs(fit_brma), info = name)
    return(invisible(TRUE))
  }

  metafor_dfbetas <- .metafor_dfbetas(fit_metafor)
  skip_rows       <- case_value(case, "skip_rows", integer())
  metafor_cols    <- case_value(case, "metafor_cols", NULL)
  brma_cols       <- case_value(case, "brma_cols", NULL)

  if (is.null(metafor_cols)) {
    metafor_cols <- seq_len(ncol(metafor_dfbetas))
  }
  if (is.null(brma_cols)) {
    brma_cols <- metafor_cols
  }

  rows <- seq_len(nrow(metafor_dfbetas))
  rows <- setdiff(rows, skip_rows)

  for (i in seq_along(metafor_cols)) {
    testthat::expect_equal(
      metafor_dfbetas[rows, metafor_cols[[i]]],
      brma_dfbetas[rows, brma_cols[[i]]],
      tolerance = tolerance,
      info      = paste(name, "DFBETAS match metafor")
    )
  }
}

expect_summary_heterogeneity_matches_metafor <- function(case) {

  name        <- case_name(case)
  fit_brma    <- .case_fit(case)
  fit_metafor <- .case_metafor(case)
  brma_het    <- summary_heterogeneity(fit_brma)
  kind        <- case_value(case, "kind", "standard")
  tolerance   <- case_value(case, "tolerance", 0.05)
  i2_tol      <- case_value(case, "i2_tolerance", tolerance)
  h2_tol      <- case_value(case, "h2_tolerance", tolerance)

  if (kind == "multilevel") {
    testthat::expect_equal(
      brma_het$estimates["tau", "Mean"],
      sqrt(fit_metafor$tau2),
      tolerance = tolerance,
      info      = paste(name, "total tau matches metafor")
    )
    testthat::expect_equal(
      brma_het$estimates["tau [within]", "Mean"],
      sqrt(fit_metafor$tau2 * (1 - fit_metafor$rho)),
      tolerance = tolerance,
      info      = paste(name, "within tau matches metafor")
    )
    testthat::expect_equal(
      brma_het$estimates["tau [between]", "Mean"],
      sqrt(fit_metafor$tau2 * fit_metafor$rho),
      tolerance = tolerance,
      info      = paste(name, "between tau matches metafor")
    )
    testthat::expect_equal(
      brma_het$estimates["rho", "Mean"],
      fit_metafor$rho,
      tolerance = tolerance,
      info      = paste(name, "rho matches metafor")
    )

    W         <- diag(1 / fit_metafor$vi)
    X         <- model.matrix(fit_metafor)
    P         <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    typical_v <- (fit_metafor$k - fit_metafor$p) / sum(diag(P))
    total_var <- fit_metafor$tau2

    i2_total   <- 100 * total_var / (total_var + typical_v)
    i2_within  <- 100 * (total_var * (1 - fit_metafor$rho)) / (total_var + typical_v)
    i2_between <- 100 * (total_var * fit_metafor$rho) / (total_var + typical_v)

    testthat::expect_equal(brma_het$estimates["I2", "Median"], i2_total,
                           tolerance = 0.05, info = paste(name, "I2 total"))
    testthat::expect_equal(brma_het$estimates["I2 [within]", "Median"], i2_within,
                           tolerance = 0.10, info = paste(name, "I2 within"))
    testthat::expect_equal(brma_het$estimates["I2 [between]", "Median"], i2_between,
                           tolerance = 0.10, info = paste(name, "I2 between"))
    return(invisible(TRUE))
  }

  if (kind == "scale") {
    testthat::expect_equal(brma_het$estimates["I2", "Median"], mean(fit_metafor$I2),
                           tolerance = i2_tol, info = paste(name, "I2"))
    testthat::expect_equal(brma_het$estimates["H2", "Median"], mean(fit_metafor$H2),
                           tolerance = h2_tol, info = paste(name, "H2"))
    return(invisible(TRUE))
  }

  testthat::expect_equal(brma_het$estimates["tau", "Mean"], sqrt(fit_metafor$tau2),
                         tolerance = tolerance, info = paste(name, "tau"))
  testthat::expect_equal(brma_het$estimates["tau2", "Mean"], fit_metafor$tau2,
                         tolerance = tolerance, info = paste(name, "tau2"))

  if (kind != "selection") {
    testthat::expect_equal(brma_het$estimates["I2", "Median"], fit_metafor$I2,
                           tolerance = i2_tol, info = paste(name, "I2"))
    testthat::expect_equal(brma_het$estimates["H2", "Median"], fit_metafor$H2,
                           tolerance = h2_tol, info = paste(name, "H2"))
  }
}

expect_influence_matches_metafor <- function(case) {

  name        <- case_name(case)
  fit_brma    <- .case_fit(case)
  fit_metafor <- .case_metafor(case)
  oracle      <- case_value(case, "oracle", "equal")

  if (oracle == "cooks_finite") {
    testthat::expect_true(
      all(is.finite(suppressWarnings(cooks.distance(fit_brma)))),
      info = paste(name, "cooks.distance is finite")
    )
    return(invisible(TRUE))
  }

  inf_metafor <- influence(fit_metafor)
  inf_brma    <- suppressWarnings(influence(fit_brma))

  if (oracle == "rank") {
    testthat::expect_true(
      stats::cor(inf_metafor$inf$rstudent, inf_brma$inf$rstudent,
                 method = "spearman", use = "pairwise") > 0.9,
      info = paste(name, "rstudent rank matches metafor")
    )
    testthat::expect_true(all(is.finite(inf_brma$inf$dffits)),
                          info = paste(name, "dffits finite"))
    testthat::expect_true(all(is.finite(inf_brma$inf$cook.d)),
                          info = paste(name, "cook.d finite"))
    testthat::expect_true(
      stats::cor(inf_metafor$inf$cov.r, inf_brma$inf$cov.r,
                 method = "spearman", use = "pairwise") > 0.9,
      info = paste(name, "cov.r rank matches metafor")
    )
    testthat::expect_true(
      stats::cor(inf_metafor$inf$tau2.del, inf_brma$inf$tau.del,
                 method = "spearman", use = "pairwise") > 0.9,
      info = paste(name, "tau.del rank matches metafor")
    )
    testthat::expect_true(
      stats::cor(inf_metafor$inf$hat, inf_brma$inf$hat,
                 method = "spearman", use = "pairwise") > 0.9,
      info = paste(name, "hat rank matches metafor")
    )
    return(invisible(TRUE))
  }

  skip_rows <- case_value(case, "skip_rows", integer())
  rows      <- setdiff(seq_len(nrow(inf_brma$inf)), skip_rows)
  tol       <- if (oracle == "equal") 0.05 else 0.10

  testthat::expect_equal(inf_metafor$inf$rstudent[rows], inf_brma$inf$rstudent[rows],
                         tolerance = if (oracle == "equal") 0.05 else 0.15,
                         info = paste(name, "rstudent"))
  testthat::expect_true(all(is.finite(inf_brma$inf$dffits[rows])),
                        info = paste(name, "dffits finite"))
  testthat::expect_true(all(is.finite(inf_brma$inf$cook.d[rows])),
                        info = paste(name, "cook.d finite"))

  if (oracle == "equal") {
    testthat::expect_equal(inf_metafor$inf$cov.r[rows], inf_brma$inf$cov.r[rows],
                           tolerance = tol, info = paste(name, "cov.r"))
  }
  testthat::expect_equal(inf_metafor$inf$tau2.del[rows], inf_brma$inf$tau.del[rows]^2,
                         tolerance = tol, info = paste(name, "tau2.del"))
  testthat::expect_equal(inf_metafor$inf$hat[rows], inf_brma$inf$hat[rows],
                         tolerance = tol, info = paste(name, "hat"))

  if (oracle == "equal") {
    testthat::expect_equal(inf_metafor$dfbs$intrcpt, inf_brma$dfbs$mu,
                           tolerance = 0.10, info = paste(name, "dfbetas"))
    testthat::expect_equal(inf_brma$inf$dffits, unname(dffits(fit_brma)),
                           tolerance = 1e-12)
    testthat::expect_equal(inf_brma$inf$cook.d, unname(cooks.distance(fit_brma)),
                           tolerance = 1e-12)
    testthat::expect_equal(inf_metafor$inf$cov.r, unname(covratio(fit_brma)),
                           tolerance = 0.05)
  }
}
