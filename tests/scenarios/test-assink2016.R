if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
REGENERATE_SCENARIO_FILES <- FALSE
SHOW_SCENARIO_OUTPUT      <- FALSE
scenario_start("assink2016")
# testthat::test_file("tests/scenarios/test-assink2016.R")

### Description
# compare multivariate nested random-effects targets with metafor

testthat::test_that("Assink multivariate nested random-effects model", {

  set.seed(1)
  data("dat.assink2016", package = "metadat")

  V_assink <- metafor::vcalc(
    vi, cluster = study, type = deltype, obs = esid,
    rho = c(0.7, 0.5), data = dat.assink2016
  )
  V_assink_diagonal <- diag(dat.assink2016[["vi"]])

  ### Model fits ----
  fit_metafor           <- metafor::rma.mv(yi, V_assink, random = ~ 1 | study / esid, data = dat.assink2016)
  fit_metafor.cs        <- metafor::rma.mv(yi, V_assink, random = ~ esid | study, data = dat.assink2016)
  fit_metafor_diag      <- metafor::rma.mv(yi, V_assink_diagonal, random = ~ 1 | study / esid, data = dat.assink2016)
  fit_metafor_diag.cs   <- metafor::rma.mv(yi, V_assink_diagonal, random = ~ esid | study, data = dat.assink2016)
  fit_metafor_no_study  <- metafor::rma.mv(yi, V_assink, random = ~ 1 | id, data = dat.assink2016)
  fit_metafor_no_effect <- metafor::rma.mv(yi, V_assink, random = ~ 1 | esid, data = dat.assink2016)

  fit_brma.mv <- scenario_fit("fit_brma.mv", {
    tmp <- brma.mv(yi = yi, V = V_assink, measure = "SMD", random = ~ 1 | study / esid, data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_brma.mv_diag <- scenario_fit("fit_brma.mv_diag", {
    tmp <- brma.mv(yi = yi, V = V_assink_diagonal, measure = "SMD", random = ~ 1 | study / esid, data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_brma.mv_no_study <- scenario_fit("fit_brma.mv_no_study", {
    tmp <- brma.mv(yi = yi, V = V_assink, measure = "SMD", random = ~ 1 | study:esid, data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_brma.mv_no_effect <- scenario_fit("fit_brma.mv_no_effect", {
    tmp <- brma.mv(yi = yi, V = V_assink, measure = "SMD", random = ~ 1 | study, data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_brma_cluster <- scenario_fit("fit_brma_cluster", {
    tmp <- brma(yi = yi, vi = vi, cluster = study, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_brma <- scenario_fit("fit_brma", {
    tmp <- brma(yi = yi, vi = vi, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })


  ### Model summary ----
  fit_metafor_diag
  fit_metafor_diag.cs
  scenario_text("summary-fit_brma.mv_diag", summary(fit_brma.mv_diag))

  fit_metafor
  fit_metafor.cs
  scenario_text("summary-fit_brma.mv", summary(fit_brma.mv))

  fit_metafor_no_study
  scenario_text("summary-fit_brma.mv_no_study", summary(fit_brma.mv_no_study))

  break # error here:
  fit_metafor_no_effect
  scenario_text("summary-fit_brma.mv_no_effect", summary(fit_brma.mv_no_effect))



  ### Fixed effect and heterogeneity ----
  fit_summary          <- summary(fit_brma, include_mcmc_diagnostics = FALSE)
  heterogeneity        <- summary_heterogeneity(fit_brma, component = "all")
  random_effects       <- ranef(fit_brma, simplify = FALSE)[["location"]]
  random_components    <- names(random_effects)
  nested_component     <- setdiff(random_components, "study")
  allocation_component <- setdiff(names(heterogeneity), random_components)

  total_heterogeneity  <- heterogeneity[[allocation_component]][["estimates"]]
  study_heterogeneity  <- heterogeneity[["study"]][["estimates"]]
  nested_heterogeneity <- heterogeneity[[nested_component]][["estimates"]]
  study_fraction       <- "var_frac(random_total: study)"
  nested_fraction      <- paste0("var_frac(random_total: ", nested_component, ")")

  scenario_text("fit_parameter_comparison", {data.frame(
    target  = c("intercept", "total heterogeneity SD", "study-level heterogeneity SD",
                "study/esid-level heterogeneity SD", "study variance fraction", "study/esid variance fraction"),
    brma.mv = c(fit_summary[["estimates_mods"]]["intercept", "Mean"], total_heterogeneity["tau", "Mean"],
                study_heterogeneity["tau", "Mean"], nested_heterogeneity["tau", "Mean"],
                total_heterogeneity[study_fraction, "Mean"], total_heterogeneity[nested_fraction, "Mean"]),
    metafor = c(as.numeric(fit_metafor[["beta"]]), sqrt(sum(fit_metafor[["sigma2"]])),
                sqrt(fit_metafor[["sigma2"]][[1L]]), sqrt(fit_metafor[["sigma2"]][[2L]]),
                fit_metafor[["sigma2"]][[1L]] / sum(fit_metafor[["sigma2"]]),
                fit_metafor[["sigma2"]][[2L]] / sum(fit_metafor[["sigma2"]]))
  )})

  ### Predictions ----
  # The explicit row requests a marginal new true effect. The terms target is the fixed location and uses metafor's confidence interval.
  newdata <- data.frame(.prediction_row = 1L)
  set.seed(1617)
  prediction_estimate <- summary(predict(fit_brma, newdata = newdata, type = "estimate", quiet = TRUE))
  prediction_terms    <- summary(predict(fit_brma, newdata = newdata, type = "terms", quiet = TRUE))
  prediction_metafor  <- predict(fit_metafor)

  scenario_text("fit_prediction_comparison", {data.frame(
    target        = c("fixed location", "marginal new true effect"),
    brma.mv_mean  = c(prediction_terms[1L, "Mean"], prediction_estimate[1L, "Mean"]),
    brma.mv_lower = c(prediction_terms[1L, "0.025"], prediction_estimate[1L, "0.025"]),
    brma.mv_upper = c(prediction_terms[1L, "0.975"], prediction_estimate[1L, "0.975"]),
    metafor_mean  = rep(prediction_metafor[["pred"]][[1L]], 2L),
    metafor_lower = c(prediction_metafor[["ci.lb"]][[1L]], prediction_metafor[["pi.lb"]][[1L]]),
    metafor_upper = c(prediction_metafor[["ci.ub"]][[1L]], prediction_metafor[["pi.ub"]][[1L]])
  )})

  ### Random effects ----
  random_metafor <- metafor::ranef(fit_metafor, expand = TRUE)
  study_mean     <- colMeans(as.matrix(random_effects[["study"]]))
  nested_mean    <- colMeans(as.matrix(random_effects[[nested_component]]))
  study_rows     <- !duplicated(dat.assink2016[["study"]])

  study_mean       <- study_mean[study_rows]
  study_metafor    <- random_metafor[["study"]][["intrcpt"]][study_rows]
  study_average    <- (study_mean + study_metafor) / 2
  study_difference <- study_mean - study_metafor
  study_bias       <- mean(study_difference)
  study_limits     <- study_bias + c(-1, 1) * 1.96 * stats::sd(study_difference)

  scenario_plot("fit_random_study_comparison", {
    plot(study_average, study_difference, main = "Study random effects", xlab = "Mean random effect (RoBMA and metafor)",
         ylab = "RoBMA - metafor", ylim = range(c(0, study_difference, study_limits)), pch = 16)
    abline(h = 0, lty = 2)
    abline(h = study_bias, lwd = 2)
    abline(h = study_limits, lty = 3)
  })

  nested_metafor    <- random_metafor[["study/esid"]][["intrcpt"]]
  nested_average    <- (nested_mean + nested_metafor) / 2
  nested_difference <- nested_mean - nested_metafor
  nested_bias       <- mean(nested_difference)
  nested_limits     <- nested_bias + c(-1, 1) * 1.96 * stats::sd(nested_difference)

  scenario_plot("fit_random_nested_comparison", {
    plot(nested_average, nested_difference, main = "Study/esid random effects", xlab = "Mean random effect (RoBMA and metafor)",
         ylab = "RoBMA - metafor", ylim = range(c(0, nested_difference, nested_limits)), pch = 16)
    abline(h = 0, lty = 2)
    abline(h = nested_bias, lwd = 2)
    abline(h = nested_limits, lty = 3)
  })

  ### Model consistency ----
  scenario_plot("fit_posterior_effect_comparison", {
    plot(fit_brma_cluster, "mu", prior = TRUE, ylim = c(0, 7), main = "Pooled effect")
    lines(fit_brma_no_cluster, "mu", col = "blue")
    lines(fit_brma, "mu", col = "red")
    legend("topright", legend = c("hierarchical, diagonal vi", "non-hierarchical, diagonal vi", "hierarchical, Vcalc"),
           col = c("black", "blue", "red"), lty = 1, bty = "n")
  })

  scenario_plot("fit_posterior_heterogeneity_comparison", {
    plot(fit_brma_cluster, "tau", prior = TRUE, xlim = c(0, 1.1), ylim = c(0, 8), main = "Total heterogeneity SD")
    lines(fit_brma_no_cluster, "tau", col = "blue")
    lines(fit_brma, "sd_total(random_total)", component = "random", col = "red")
    legend("topright", legend = c("hierarchical, diagonal vi", "non-hierarchical, diagonal vi", "hierarchical, Vcalc"),
           col = c("black", "blue", "red"), lty = 1, bty = "n")
  })

  ### Heterogeneity allocation ----
  # These are matching study-level heterogeneity variance fractions. The rho used to construct V_assink fixes sampling covariance and is not estimated.
  # At rho = 0, the clustered model reduces to ordinary brma; at rho = 1, it reduces to a study-only brma.mv model with diagonal sampling covariance.
  set.seed(1616)
  BF_rho0_qCMDE <- suppressWarnings(hypothesis(fit_brma_cluster, "rho = 0", density_method = "qCMDE", density_control = list(samples = 12000)))
  BF_rho1_qCMDE <- suppressWarnings(hypothesis(fit_brma_cluster, "rho = 1", density_method = "qCMDE", density_control = list(samples = 12000)))

  scenario_text("fit_rho_bayes_factor_comparison", {data.frame(
    null            = c("Vcalc: no study random effect", "diagonal vi: rho = 0", "diagonal vi: rho = 1"),
    bridge_sampling = c(bf(fit_brma, fit_brma_no_study)[["bf"]],
                        bridgesampling::bf(bridge_sampler(fit_brma_cluster), bridge_sampler(fit_brma_no_cluster))[["bf"]],
                        bridgesampling::bf(bridge_sampler(fit_brma_cluster), bridge_sampler(fit_brma_study_only))[["bf"]]),
    qCMDE_density   = c(NA_real_, BF_rho0_qCMDE[["BF"]], BF_rho1_qCMDE[["BF"]]),
    qCMDE_error     = c(NA_real_, BF_rho0_qCMDE[["BF_error"]], BF_rho1_qCMDE[["BF_error"]])
  )})

  scenario_plot("fit_posterior_rho_comparison", {
    plot(fit_brma_cluster, "rho", prior = TRUE, density_method = "KDE", xlim = c(0, 1),
         xlab = "Study variance fraction (rho)", ylab = "Posterior density", main = "Heterogeneity allocation density")
    lines(fit_brma, "var_frac(random_total: study)", component = "random", density_method = "KDE", col = "red")
    legend("topright", legend = c("diagonal vi: rho", "Vcalc: study variance fraction", "uniform prior"),
           col = c("black", "red", "gray"), lty = 1, bty = "n")
  })
})
