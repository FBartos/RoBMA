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
  scenario_text("summary-heterogeneity", summary_heterogeneity(fit_brma))

  brma_tau <- function(fit, component) summary_heterogeneity(fit, component = component)[["estimates"]]["tau", "Mean"]
  scenario_text("metafor-comparison", data.frame(
    model           = rep(c("study+observation", "study", "observation"), each = 2L),
    implementation  = rep(c("metafor", "RoBMA"), 3L),
    intercept        = c(as.numeric(stats::coef(fit_metafor)), as.numeric(stats::coef(fit_brma)[[1L]]),
                         as.numeric(stats::coef(fit_metafor_study)), as.numeric(stats::coef(fit_brma_study)[[1L]]),
                         as.numeric(stats::coef(fit_metafor_observation)), as.numeric(stats::coef(fit_brma_observation)[[1L]])),
    study_SD         = c(sqrt(fit_metafor[["sigma2"]][[1L]]), brma_tau(fit_brma, "study"),
                         sqrt(fit_metafor_study[["sigma2"]][[1L]]), brma_tau(fit_brma_study, "study"), NA, NA),
    observation_SD   = c(sqrt(fit_metafor[["sigma2"]][[2L]]), brma_tau(fit_brma, "observation"),
                         NA, NA, sqrt(fit_metafor_observation[["sigma2"]][[1L]]), brma_tau(fit_brma_observation, "observation")),
    total_random_SD  = c(sqrt(sum(fit_metafor[["sigma2"]])), brma_tau(fit_brma, "total"),
                         sqrt(fit_metafor_study[["sigma2"]][[1L]]), brma_tau(fit_brma_study, "study"),
                         sqrt(fit_metafor_observation[["sigma2"]][[1L]]), brma_tau(fit_brma_observation, "observation"))
  ))

  get_looic <- function(fit) loo(fit)[["estimates"]]["looic", 1L]
  scenario_text("model-fit", data.frame(
    model                   = c("study+observation", "study", "observation"),
    log_marginal_likelihood = c(logml(fit_brma), logml(fit_brma_study), logml(fit_brma_observation)),
    looic                   = c(get_looic(fit_brma), get_looic(fit_brma_study), get_looic(fit_brma_observation))
  ))

  ### hypothesis and model comparison ----
  scenario_text("random-component-bayes-factors", {
    density_bf <- hypothesis(
      fit_brma,
      c("var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1"),
      density_method  = "qCMDE",
      density_control = list(samples = 2000L)
    )
    # Both columns are BF(full / nested): the component is present vs absent.
    data.frame(
      omitted_component = c("study", "observation"),
      qCMDE_BF           = density_bf[["BF"]],
      bridge_BF          = c(bf(fit_brma, fit_brma_observation)[["bf"]], bf(fit_brma, fit_brma_study)[["bf"]])
    )
  })

  loo_performance <- function(fit) {

    fit_loo  <- loo(fit)
    pareto_k <- fit_loo[["diagnostics"]][["pareto_k"]]
    return(c(
      elpd_loo           = fit_loo[["estimates"]]["elpd_loo", 1L],
      elpd_loo_SE        = fit_loo[["estimates"]]["elpd_loo", 2L],
      p_loo              = fit_loo[["estimates"]]["p_loo", 1L],
      looic              = fit_loo[["estimates"]]["looic", 1L],
      maximum_pareto_k   = max(pareto_k),
      pareto_k_above_0.7 = sum(pareto_k > 0.7),
      pareto_k_above_1   = sum(pareto_k > 1)
    ))
  }
  scenario_text("loo-performance", rbind(
    `study+observation` = loo_performance(fit_brma),
    study               = loo_performance(fit_brma_study),
    observation         = loo_performance(fit_brma_observation)
  ))
  scenario_text("loo-model-weights", loo_model_weights(fit_brma, fit_brma_study, fit_brma_observation))

  ### random-effect comparisons ----
  ranef_metafor          <- metafor::ranef(fit_metafor)
  ranef_metafor_expanded <- metafor::ranef(fit_metafor, expand = TRUE)
  ranef_brma             <- ranef(fit_brma, simplify = FALSE)
  ranef_brma_total       <- as.data.frame(ranef(fit_brma, component = "total", expand = TRUE))[["Mean"]]

  component_means <- function(x, component) {

    output        <- as.data.frame(x)[["Mean"]]
    names(output) <- sub(paste0("^u_", component, "\\[([^]]+)\\]$"), "\\1", rownames(as.data.frame(x)))
    return(output)
  }
  ranef_brma_study       <- component_means(ranef_brma[["study"]], "study")
  ranef_brma_observation <- component_means(ranef_brma[["observation"]], "observation")
  ranef_metafor_total    <- ranef_metafor_expanded[["study_id"]][["intrcpt"]] + ranef_metafor_expanded[["obs"]][["intrcpt"]]

  agreement_summary <- function(reference, estimate) {

    difference <- estimate - reference
    return(data.frame(
      mean_difference         = mean(difference),
      root_mean_square_error  = sqrt(mean(difference^2)),
      maximum_abs_difference  = max(abs(difference)),
      correlation             = stats::cor(reference, estimate)
    ))
  }
  scenario_text("ranef-agreement", rbind(
    study       = agreement_summary(ranef_metafor[["study_id"]][["intrcpt"]], ranef_brma_study[rownames(ranef_metafor[["study_id"]])]),
    observation = agreement_summary(ranef_metafor[["obs"]][["intrcpt"]], ranef_brma_observation[rownames(ranef_metafor[["obs"]])]),
    combined    = agreement_summary(ranef_metafor_total, ranef_brma_total)
  ))

  scenario_plot("ranef-comparison", {
    par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(ranef_metafor[["study_id"]][["intrcpt"]], ranef_brma_study[rownames(ranef_metafor[["study_id"]])], "Study effects")
    scenario_agreement_plot(ranef_metafor[["obs"]][["intrcpt"]], ranef_brma_observation[rownames(ranef_metafor[["obs"]])], "Observation effects")
    scenario_agreement_plot(ranef_metafor_total, ranef_brma_total, "Combined effects")
  })
})
