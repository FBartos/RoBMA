if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("ishak2007")
# testthat::test_file("tests/scenarios/test-ishak2007.R")

### Description
# compare heterogeneous and homogeneous AR(1) random-effects models with the
# longitudinal dat.ishak2007 example from metadat and the corresponding metafor
# models. The AR model is the equal-variance restriction of HAR.

testthat::test_that("Ishak longitudinal heterogeneous AR model", {

  set.seed(1)
  data("dat.ishak2007", package = "metadat")

  dat <- reshape(dat.ishak2007, direction = "long", idvar = "study", v.names = c("yi", "vi"), varying = list(c(2, 4, 6, 8), c(3, 5, 7, 9)))
  dat <- dat[order(dat$study, dat$time), ]
  dat <- dat[!is.na(dat$yi), ]
  rownames(dat) <- NULL
  dat$time_factor <- factor(dat$time)

  # Follow the metadat example: the sampling errors use the AR(1) correlation
  # 0.97 reported by Ishak et al. (2007); HAR/AR below describe true effects.
  V <- metafor::vcalc(vi, cluster = study, time1 = time, phi = 0.97, data = dat)

  unit_information_sd <- estimate_unit_information_sd(sei = sqrt(14.3), ni = 15)
  time_prior         <- prior_factor(distribution = "normal", parameters = list(mean = 0, sd = unit_information_sd), contrast = "independent")
  random_prior_har0  <- prior_random(study = random_block(covariance = random_covariance(cor = prior("spike", list(0)))))
  convergence_checks <- set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
  qcmde_control      <- list(samples = 1000L)

  ### Model fits ----
  fit_metafor_har  <- metafor::rma.mv(yi, V, mods = ~ 0 + time_factor, random = ~ time | study, struct = "HAR", data = dat)
  fit_metafor_har0 <- metafor::rma.mv(yi, V, mods = ~ 0 + time_factor, random = ~ time | study, struct = "HAR", rho = 0, data = dat)
  fit_metafor_ar   <- metafor::rma.mv(yi, V, mods = ~ 0 + time_factor, random = ~ time | study, struct = "AR",  data = dat)

  fit_brma_har <- scenario_fit("fit_brma_har", {
    temp_fit <- brma.mv(yi = yi, V = V, measure = "GEN", mods = ~ 0 + time_factor, random = ~ har(time | study),
                        prior_effect = NULL, prior_mods = time_prior, prior_unit_information_sd = unit_information_sd, data = dat,
                        chains = 3, sample = 10000, burnin = 4000, adapt = 2000, seed = 3, silent = FALSE, convergence_checks = convergence_checks)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_har0 <- scenario_fit("fit_brma_har0", {
    temp_fit <- brma.mv(yi = yi, V = V, measure = "GEN", mods = ~ 0 + time_factor, random = ~ har(time | study),
                        prior_effect = NULL, prior_mods = time_prior, prior_unit_information_sd = unit_information_sd, prior_heterogeneity = random_prior_har0, data = dat,
                        chains = 3, sample = 10000, burnin = 4000, adapt = 2000, seed = 5, silent = FALSE, convergence_checks = convergence_checks)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_ar <- scenario_fit("fit_brma_ar", {
    temp_fit <- brma.mv(yi = yi, V = V, measure = "GEN", mods = ~ 0 + time_factor, random = ~ ar(time | study),
                        prior_effect = NULL, prior_mods = time_prior, prior_unit_information_sd = unit_information_sd, data = dat,
                        chains = 3, sample = 10000, burnin = 4000, adapt = 2000, seed = 4, silent = FALSE, convergence_checks = convergence_checks)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })

  ### Model summaries ----
  fit_metafor_har
  scenario_text("summary-har", summary(fit_brma_har))
  scenario_text("summary-heterogeneity-har", summary_heterogeneity(fit_brma_har))

  fit_metafor_ar
  scenario_text("summary-ar", summary(fit_brma_ar))
  scenario_text("summary-heterogeneity-ar", summary_heterogeneity(fit_brma_ar))

  metafor_parameters_har <- c(time_1 = "time_factor[1]", time_2 = "time_factor[2]", time_3 = "time_factor[3]", time_4 = "time_factor[4]", sd_1 = "tau[1]", sd_2 = "tau[2]", sd_3 = "tau[3]", sd_4 = "tau[4]", correlation = "rho")
  robma_parameters_har   <- c(time_1 = "time_factor[1]", time_2 = "time_factor[2]", time_3 = "time_factor[3]", time_4 = "time_factor[4]", sd_1 = "sd(time[1])", sd_2 = "sd(time[2])", sd_3 = "sd(time[3])", sd_4 = "sd(time[4])", correlation = "cor")
  robma_components_har   <- c(rep(NA_character_, 4L), rep("study", 5L))
  metafor_parameters_ar  <- c(time_1 = "time_factor[1]", time_2 = "time_factor[2]", time_3 = "time_factor[3]", time_4 = "time_factor[4]", sd_1 = "tau", sd_2 = "tau", sd_3 = "tau", sd_4 = "tau", correlation = "rho")
  robma_parameters_ar    <- c(time_1 = "time_factor[1]", time_2 = "time_factor[2]", time_3 = "time_factor[3]", time_4 = "time_factor[4]", sd_1 = "sd", sd_2 = "sd", sd_3 = "sd", sd_4 = "sd", correlation = "cor")
  robma_components_ar    <- c(rep(NA_character_, 4L), rep("study", 4L), NA_character_)
  scenario_text("metafor-comparison", data.frame(
    model          = rep(c("HAR", "AR"), each = 2L),
    implementation = rep(c("metafor", "RoBMA"), 2L),
    rbind(
      ex_m(fit_metafor_har, metafor_parameters_har), ex_r(fit_brma_har, robma_parameters_har, component = robma_components_har),
      ex_m(fit_metafor_ar, metafor_parameters_ar),   ex_r(fit_brma_ar, robma_parameters_ar, component = robma_components_ar)
    ),
    row.names = NULL
  ))

  ### Nested HAR versus AR comparison ----
  # Metafor compares the equal-variance AR restriction with HAR by an LRT;
  # RoBMA compares the same two structures by marginal and predictive fit. The
  # estimate-unit LOO weights are retained despite poor Pareto-k diagnostics.
  scenario_text("metafor-structure-comparison", anova(fit_metafor_ar, fit_metafor_har))
  scenario_text("robma-structure-comparison-marglik", round(post_prob(fit_brma_ar, fit_brma_har), 3))
  scenario_text("robma-structure-comparison-loo", loo_model_weights(fit_brma_ar, fit_brma_har))

  ### Correlation density test ----
  # Compare the HAR model against the same model with cor fixed to zero. The
  # bridge Bayes factor is the model-based reference for the density estimates.
  anova(fit_metafor_har0, fit_metafor_har)
  scenario_text("robma-correlation-density-test", {
    bf_bridge <- bf(fit_brma_har, fit_brma_har0)
    bf_kde    <- hypothesis(fit_brma_har, "cor = 0")
    bf_qcmde  <- hypothesis(fit_brma_har, "cor = 0", density_method = "qCMDE", density_control = qcmde_control)

    data.frame(bridge = bf_bridge$bf, KDE = bf_kde$BF, qCMDE = bf_qcmde$BF)
  })
  scenario_text("robma-correlation-comparison-loo", loo_model_weights(fit_brma_har0, fit_brma_har))

  ### simple plots ----
  scenario_plot("time_factor", {
    plot(fit_brma_har, "time_factor", prior = TRUE)
    lines(fit_brma_har, "time_factor", density_method = "qCMDE", lty = 2)

    lines(fit_brma_ar, "time_factor", col = c("green", "blue"))
    lines(fit_brma_ar, "time_factor", col = c("green", "blue"), lty = 2, density_method = "qCMDE")
  })

  scenario_plot("sd_common", {
    plot(fit_brma_har, "sd_common", prior = TRUE)
    lines(fit_brma_har, "sd_common", density_method = "qCMDE", lty = 2)

    lines(fit_brma_ar, "sd_common", col = "blue")
    lines(fit_brma_ar, "sd_common", col = "blue", lty = 2, density_method = "qCMDE")
  })

  ### Random-parameter plots ----
  scenario_plot("random-parameters", {
    par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))

    plot(fit_brma_har, "sd_common", prior = TRUE, main = "Common SD: HAR and AR")
    lines(fit_brma_har, "sd_common", density_method = "qCMDE", density_control = qcmde_control, lty = 2)
    lines(fit_brma_ar, "sd", col = "blue")
    lines(fit_brma_ar, "sd", density_method = "qCMDE", density_control = qcmde_control, col = "blue", lty = 2)
    legend("topright", legend = c("HAR KDE", "HAR qCMDE", "AR KDE", "AR qCMDE"), col = c("black", "black", "blue", "blue"), lty = c(1, 2, 1, 2), bty = "n", cex = 0.7)

    plot(fit_brma_har, "sd(time[1])", prior = TRUE, main = "HAR time 1 SD")
    lines(fit_brma_har, "sd(time[1])", density_method = "qCMDE", density_control = qcmde_control, lty = 2)

    plot(fit_brma_har, "sd(time[2])", prior = TRUE, main = "HAR time 2 SD")
    lines(fit_brma_har, "sd(time[2])", density_method = "qCMDE", density_control = qcmde_control, lty = 2)

    plot(fit_brma_har, "sd(time[3])", prior = TRUE, main = "HAR time 3 SD")
    lines(fit_brma_har, "sd(time[3])", density_method = "qCMDE", density_control = qcmde_control, lty = 2)

    plot(fit_brma_har, "sd(time[4])", prior = TRUE, main = "HAR time 4 SD")
    lines(fit_brma_har, "sd(time[4])", density_method = "qCMDE", density_control = qcmde_control, lty = 2)

    plot(fit_brma_har, "cor", prior = TRUE, main = "Correlation: HAR and AR")
    lines(fit_brma_har, "cor", density_method = "qCMDE", density_control = qcmde_control, lty = 2)
    lines(fit_brma_ar, "cor", col = "blue")
    lines(fit_brma_ar, "cor", density_method = "qCMDE", density_control = qcmde_control, col = "blue", lty = 2)
  })

  ### Random effects ----
  ranef_metafor_har <- metafor::ranef(fit_metafor_har)[[1L]][["intrcpt"]]
  ranef_metafor_ar  <- metafor::ranef(fit_metafor_ar)[[1L]][["intrcpt"]]
  scenario_plot("ranef-agreement", {
    ranef_brma_har <- as.data.frame(ranef(fit_brma_har, component = "study"))[["Mean"]]
    ranef_brma_ar  <- as.data.frame(ranef(fit_brma_ar, component = "study"))[["Mean"]]

    par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(ranef_metafor_har, ranef_brma_har, "HAR: RoBMA vs metafor")
    scenario_agreement_plot(ranef_metafor_ar,  ranef_brma_ar,  "AR: RoBMA vs metafor")
    scenario_agreement_plot(ranef_brma_ar, ranef_brma_har, "RoBMA: HAR vs AR", reference_label = "RoBMA AR", estimate_label = "RoBMA HAR")
  })

  ### random-effect comparisons ----
  scenario_plot("ranef-comparison", {
    ranef_metafor <- metafor::ranef(ranef_metafor_har)
    ranef_brma    <- ranef(fit_brma_har)

    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(ranef_metafor[["study_id"]][["intrcpt"]], as.data.frame(ranef_brma$study)[["Mean"]], "Study effects")
    scenario_agreement_plot(ranef_metafor[["obs"]][["intrcpt"]], as.data.frame(ranef_brma$observation)[["Mean"]], "Study effects")
  })

  ### diagnostics ----
  scenario_plot("marginal_diagnostics", plot_marginal_diagnostics(ranef_metafor_har, fit_brma_har))

  ### diagnostic plots ----
  set.seed(1)
  scenario_plot("funnel", funnel(fit_brma_har))
  scenario_plot("qqnorm", qqnorm(fit_brma_har))
  scenario_plot("zplot",  zplot(fit_brma_har))
})
