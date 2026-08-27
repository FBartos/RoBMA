if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("white2020")
# testthat::test_file("tests/scenarios/test-white2020.R")

### Description
# compare two independent random-intercept structures and the nested models
# omitting each component with the documented metafor analyses

testthat::test_that("White study and observation random-effects model", {

  set.seed(1)
  data("dat.white2020", package = "metadat")

  dat <- metafor::escalc(measure = "ZCOR", ri = r, ni = n, data = dat.white2020)

  ### model fits ----
  fit_metafor             <- metafor::rma.mv(yi, vi, random = list(~ 1 | study_id, ~ 1 | obs), data = dat)
  fit_metafor_study       <- metafor::rma.mv(yi, vi, random = ~ 1 | study_id, data = dat)
  fit_metafor_observation <- metafor::rma.mv(yi, vi, random = ~ 1 | obs, data = dat)

  fit_brma <- scenario_fit("fit_brma", {
    temp_fit <- brma.mv(yi = yi, V = vi, ni = n, random = list(study = ~ 1 | study_id, observation = ~ 1 | obs), data = dat, measure = "ZCOR", seed = 1)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_study <- scenario_fit("fit_brma_study", {
    temp_fit <- brma.mv(yi = yi, V = vi, ni = n, random = list(study = ~ 1 | study_id), data = dat, measure = "ZCOR", seed = 1)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_observation <- scenario_fit("fit_brma_observation", {
    temp_fit <- brma.mv(yi = yi, V = vi, ni = n, random = list(observation = ~ 1 | obs), data = dat, measure = "ZCOR", seed = 1)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })

  ### model summaries ----
  fit_metafor
  scenario_text("summary-fit_brma", summary(fit_brma, include_mcmc_diagnostics = FALSE))
  fit_metafor_study
  scenario_text("summary-fit_brma-study", summary(fit_brma_study, include_mcmc_diagnostics = FALSE))
  fit_metafor_observation
  scenario_text("summary-fit_brma-observation", summary(fit_brma_observation, include_mcmc_diagnostics = FALSE))

  scenario_text("summary-heterogeneity",             summary_heterogeneity(fit_brma))
  scenario_text("summary-heterogeneity-study",       summary_heterogeneity(fit_brma_study))
  scenario_text("summary-heterogeneity-observation", summary_heterogeneity(fit_brma_observation))

  metafor_parameters <- c("intercept", study_SD = "sigma[study_id]", observation_SD = "sigma[obs]", total_random_SD = "sigma[total]")
  robma_parameters        <- c("intercept", study_SD = "sd", observation_SD = "sd", total_random_SD = "sd_total")
  robma_single_parameters <- c("intercept", study_SD = "sd", observation_SD = "sd", total_random_SD = "sd")
  robma_components        <- c(NA, "study", "observation", NA)
  robma_study_components  <- c(NA, "study", "observation", "study")
  robma_obs_components    <- c(NA, "study", "observation", "observation")
  comparison <- rbind(
    ex_m(fit_metafor, metafor_parameters),             ex_r(fit_brma, robma_parameters, component = robma_components),
    ex_m(fit_metafor_study, metafor_parameters),       ex_r(fit_brma_study, robma_single_parameters, component = robma_study_components),
    ex_m(fit_metafor_observation, metafor_parameters), ex_r(fit_brma_observation, robma_single_parameters, component = robma_obs_components)
  )
  scenario_text("metafor-comparison", data.frame(
    model           = rep(c("study+observation", "study", "observation"), each = 2L),
    implementation  = rep(c("metafor", "RoBMA"), 3L),
    comparison,
    row.names = NULL
  ))


  scenario_text("model-fit-bma", round(post_prob(fit_brma, fit_brma_study, fit_brma_observation), 3))
  scenario_text("model-fit-loo", loo_model_weights(fit_brma, fit_brma_study, fit_brma_observation))

  ### parameter plots ----
  scenario_plot("random", {
    par(mfrow = c(3, 2))

    plot(fit_brma, "sd_total", prior = TRUE)
    lines(fit_brma, "sd_total", density_method = "qCMDE", lty = 2)

    plot(fit_brma, "var_total", prior = TRUE, xlim = c(0, 0.5), ylim = c(0, 50))
    lines(fit_brma, "var_total", density_method = "qCMDE", lty = 2)

    plot(fit_brma, "var_prop(study)", prior = TRUE)
    lines(fit_brma, "var_prop(study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma, "var_prop(observation)", prior = TRUE)
    lines(fit_brma, "var_prop(observation)", density_method = "qCMDE", lty = 2)

    plot(fit_brma, "study: sd", prior = TRUE)
    lines(fit_brma, "study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000))

    plot(fit_brma, "observation: sd", prior = TRUE)
    lines(fit_brma, "observation: sd", density_method = "qCMDE", lty = 2)
  })


  ### hypothesis and model comparison ----
  scenario_text("random-component-bayes-factors", {
    density_bf <- hypothesis(
      fit_brma, c("var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1"),
      density_method  = "qCMDE", density_control = list(samples = 2000L)
    )
    data.frame(
      omitted_component = c("study", "observation"),
      qCMDE_BF           = density_bf[["BF"]],
      bridge_BF          = c(bf(fit_brma, fit_brma_observation)[["bf"]], bf(fit_brma, fit_brma_study)[["bf"]])
    )
  })


  ### random-effect comparisons ----
  scenario_plot("ranef-comparison", {
    ranef_metafor <- metafor::ranef(fit_metafor)
    ranef_brma    <- ranef(fit_brma)

    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(ranef_metafor[["study_id"]][["intrcpt"]], as.data.frame(ranef_brma$study)[["Mean"]], "Study effects")
    scenario_agreement_plot(ranef_metafor[["obs"]][["intrcpt"]], as.data.frame(ranef_brma$observation)[["Mean"]], "Study effects")
  })

  ### diagnostics ----
  scenario_plot("marginal_diagnostics",   suppressWarnings(plot_marginal_diagnostics(fit_metafor, fit_brma)))

  ### diagnostic plots ----
  set.seed(1)
  scenario_plot("funnel", funnel(fit_brma))
  scenario_plot("qqnorm", qqnorm(fit_brma))
  scenario_plot("zplot",  zplot(fit_brma))
})
