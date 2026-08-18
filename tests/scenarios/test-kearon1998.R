if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("kearon1998")
# testthat::test_file("tests/scenarios/test-kearon1998.R")

### Description
# compare a bivariate unstructured random-effects model with metafor and with
# the corresponding model that omits the random-effects correlation

testthat::test_that("Kearon bivariate diagnostic-accuracy model", {

  set.seed(1)
  data("dat.kearon1998", package = "metadat")

  dat <- suppressWarnings(metafor::to.long(measure = "OR", ai = tp, n1i = np, ci = tn, n2i = nn, data = dat.kearon1998, subset = patients == "asymptomatic", append = FALSE))
  dat$group <- factor(dat$group, levels = c(1, 2), labels = c("sensitivity", "specificity"))
  dat <- metafor::escalc(measure = "PLO", xi = out1, mi = out2, data = dat, add = 1 / 2, to = "all")

  unit_information_sd <- estimate_unit_information_sd(sei = sqrt(dat$vi), ni = dat$out1 + dat$out2)
  group_prior         <- prior_factor(distribution = "normal", parameters = list(mean = 0, sd = unit_information_sd), contrast = "independent")
  heterogeneity_prior <- prior("normal", parameters = list(mean = 0, sd = unit_information_sd * 0.25), truncation = list(0, Inf))
  random_prior_us     <- prior_random(study = random_block(contrasts = c(group = "independent")))
  random_prior_diag   <- prior_random(study = random_block(contrasts = c(group = "independent")))
  random_prior_hcs0   <- prior_random(study = random_block(covariance = random_covariance(cor = prior("spike", list(0)))))

  ### Model fits ----
  fit_metafor_us   <- metafor::rma.mv(yi, vi, mods = ~ 0 + group, random = ~ group | study, struct = "UN", data = dat)
  fit_metafor_diag <- metafor::rma.mv(yi, vi, mods = ~ 0 + group, random = ~ group | study, struct = "DIAG", data = dat)

  fit_brma_us <- scenario_fit("fit_brma_us", {
    temp_fit <- brma.mv(yi = yi, V = vi, ni = out1 + out2, mods = ~ 0 + group, random = ~ us(0 + group | study),
                        set_contrast_factor_predictors = "independent", prior_heterogeneity = random_prior_us,
                        data = dat, measure = "GEN", seed = 1)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_hcs <- scenario_fit("fit_brma_hcs", {
    temp_fit <- brma.mv(yi = yi, V = vi, ni = out1 + out2, mods = ~ 0 + group, random = ~ hcs(group | study),
                        set_contrast_factor_predictors = "independent",
                        data = dat, measure = "GEN", seed = 1)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_diag <- scenario_fit("fit_brma_diag", {
    temp_fit <- brma.mv(yi = yi, V = vi, ni = out1 + out2, mods = ~ 0 + group, random = ~ diag(0 + group | study),
                        set_contrast_factor_predictors = "independent", prior_heterogeneity = random_prior_diag,
                        data = dat, measure = "GEN", seed = 1)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })
  fit_brma_hcs0 <- scenario_fit("fit_brma_hcs0", {
    temp_fit <- brma.mv(yi = yi, V = vi, ni = out1 + out2, mods = ~ 0 + group, random = ~ hcs(group | study),
                        set_contrast_factor_predictors = "independent", prior_heterogeneity = random_prior_hcs0,
                        data = dat, measure = "GEN", seed = 1)
    temp_fit <- add_loo(temp_fit)
    temp_fit <- add_marglik(temp_fit)
    return(temp_fit)
  })

  ### model summaries ----
  fit_metafor_us
  scenario_text("summary-fit_brma_us",  summary(fit_brma_us))
  scenario_text("summary-fit_brma_hcs", summary(fit_brma_hcs))
  fit_metafor_diag
  scenario_text("summary-fit_brma_diag", summary(fit_brma_diag))
  scenario_text("summary-fit_brma_hcs0", summary(fit_brma_hcs0))

  summary_heterogeneity(fit_brma_us)
  summary(marginal_means(fit_brma_us))
  summary(marginal_means(fit_brma_hcs))
  summary(marginal_means(fit_brma_diag))
  summary(marginal_means(fit_brma_hcs0))

  ### hypothesis ----
  set.seed(1)
  qcmde_control <- list(samples = 1000L)
  bf_hyp_us   <- hypothesis(fit_brma_us,  "study: cor(group[sensitivity],group[specificity]) = 0")
  bf_hyp_hcs  <- hypothesis(fit_brma_hcs, "study: cor = 0")
  bf_hyp_us_qCMDE  <- hypothesis(fit_brma_us,  "study: cor(group[sensitivity],group[specificity]) = 0", density_method = "qCMDE", density_control = qcmde_control)
  bf_hyp_hcs_qCMDE <- hypothesis(fit_brma_hcs, "study: cor = 0", density_method = "qCMDE", density_control = qcmde_control)
  bf_hcs_hcs0  <- bf(fit_brma_hcs, fit_brma_hcs0)
  bf_us_diag   <- bf(fit_brma_us, fit_brma_diag)

  anova(fit_metafor_us, fit_metafor_diag)
  scenario_text("bridge_density_eq", data.frame(
    bridge = c(bf_us_diag$bf, bf_hcs_hcs0$bf),
    KDE    = c(bf_hyp_us$BF,  bf_hyp_hcs$BF),
    qCMDE  = c(bf_hyp_us_qCMDE$BF,  bf_hyp_hcs_qCMDE$BF)))
  scenario_text("bf_eq_1", bf(fit_brma_hcs0, fit_brma_diag))  # these models are equivalent
  scenario_text("bf_eq_0", bf(fit_brma_hcs, fit_brma_us))     # these models are equivalent

  ### simple plots ----
  scenario_plot("mu", {
    plot(fit_brma_us, "group", prior = TRUE)
    lines(fit_brma_us, "group", density_method = "qCMDE", lty = 2)

    lines(fit_brma_hcs, "group", col = c("green", "blue"))
    lines(fit_brma_hcs, "group", col = c("green", "blue"), lty = 2, density_method = "qCMDE")
  })

  scenario_plot("sd_common", {
    plot(fit_brma_us, "sd_common", prior = TRUE)
    lines(fit_brma_us, "sd_common", density_method = "qCMDE", lty = 2)

    lines(fit_brma_hcs, "sd_common", col = "blue")
    lines(fit_brma_hcs, "sd_common", col = "blue", lty = 2, density_method = "qCMDE")
  })

  scenario_plot("cor", {
    plot(fit_brma_us, "study: cor(group[sensitivity],group[specificity])", prior = TRUE)
    lines(fit_brma_us, "study: cor(group[sensitivity],group[specificity])", density_method = "qCMDE", lty = 2)

    lines(fit_brma_hcs, "study: cor(group[sensitivity],group[specificity])", col = "blue")
    lines(fit_brma_hcs, "study: cor(group[sensitivity],group[specificity])", col = "blue", lty = 2, density_method = "qCMDE")
  })

  # other random parameters
  summary_heterogeneity(fit_brma_us)
  summary_heterogeneity(fit_brma_hcs)

  scenario_plot("random_us", {
    par(mfrow = c(3, 2)) # HERE extend

    plot(fit_brma_us, "sd_common", prior = TRUE)
    plot(fit_brma_us, "var_common", prior = TRUE)

    plot(fit_brma_us, "study: sd(group[sensitivity])", prior = TRUE)  #TODO: these should allow analytical prior visualization
    plot(fit_brma_us, "study: sd(group[specificity])", prior = TRUE)  #TODO: these should allow analytical prior visualization

    plot(fit_brma_us, "var_ratio(group[sensitivity])", prior = TRUE) #TODO: these should allow analytical prior visualization
    plot(fit_brma_us, "var_ratio(group[specificity])", prior = TRUE) #TODO: these should allow analytical prior visualization
  })

  scenario_plot("random_hcs", {
    par(mfrow = c(3, 2)) # HERE extend

    plot(fit_brma_hcs, "sd_common", prior = TRUE)
    plot(fit_brma_hcs, "var_common", prior = TRUE)

    plot(fit_brma_hcs, "study: sd(group[sensitivity])", prior = TRUE)  #TODO: these should allow analytical prior visualization
    plot(fit_brma_hcs, "study: sd(group[specificity])", prior = TRUE)  #TODO: these should allow analytical prior visualization

    plot(fit_brma_hcs, "var_ratio(group[sensitivity])", prior = TRUE) #TODO: these should allow analytical prior visualization
    plot(fit_brma_hcs, "var_ratio(group[specificity])", prior = TRUE) #TODO: these should allow analytical prior visualization
  })

  ### random-effect comparisons ----
  ranef_metafor_us   <- metafor::ranef(fit_metafor_us)[[1L]]
  ranef_metafor_diag <- metafor::ranef(fit_metafor_diag)[[1L]]
  ranef_brma_us      <- ranef(fit_brma_us)
  ranef_brma_hcs     <- ranef(fit_brma_hcs)
  ranef_brma_diag    <- ranef(fit_brma_diag)
  ranef_brma_hcs0    <- ranef(fit_brma_hcs0)

  scenario_plot("ranef-structure-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    # one of the groups is shrunk towards zero becuase the priors are too narrow for the data (systematic errors in first row)
    scenario_agreement_plot(ranef_metafor_us[["intrcpt"]],   as.data.frame(ranef_brma_us)[["Mean"]], "UN: RoBMA vs metafor")
    scenario_agreement_plot(ranef_metafor_diag[["intrcpt"]], as.data.frame(ranef_brma_diag)[["Mean"]], "DIAG: RoBMA vs metafor")
    scenario_agreement_plot(as.data.frame(ranef_brma_us)[["Mean"]],  as.data.frame(ranef_brma_hcs)[["Mean"]], "RoBMA: US vs HCS")
    scenario_agreement_plot(as.data.frame(ranef_brma_diag)[["Mean"]], as.data.frame(ranef_brma_hcs0)[["Mean"]], "RoBMA: DIAG vs HCS0")
  })


  ### diagnostics ----
  plot_marginal_diagnostics <- function(fit_metafor, fit_robma) {
    metafor_values <- list(
      "Residuals"      = as.numeric(stats::residuals(fit_metafor)),
      "Rstandard"      = stats::rstandard(fit_metafor)[["z"]],
      "Hat values"     = as.numeric(stats::hatvalues(fit_metafor)),
      "Cooks distance" = cooks.distance(fit_metafor),
      "DFBETAS"        = unlist(dfbetas(fit_metafor))
    )

    robma_values <- list(
      "Residuals"      = residuals(fit_robma, type = "outcome", conditioning_depth = "marginal"),
      "Rstandard"      = suppressWarnings(rstandard(fit_robma, conditioning_depth = "marginal"))[["z"]],
      "Hat values"     = suppressWarnings(hatvalues(fit_robma)),
      "Cooks distance" = cooks.distance(fit_robma),
      "DFBETAS"        = unlist(dfbetas(fit_robma))
    )

    par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
    for (diagnostic in names(metafor_values)) {
      scenario_agreement_plot(
        metafor_values[[diagnostic]], robma_values[[diagnostic]], diagnostic
      )
    }

    return(invisible(NULL))
  }

  scenario_plot("marginal_diagnostics_us",   plot_marginal_diagnostics(fit_metafor_us, fit_brma_us))
  scenario_plot("marginal_diagnostics_diag", plot_marginal_diagnostics(fit_metafor_diag, fit_brma_diag))
  # the main differences bellow seem to be due to unstable loo
  scenario_plot("robma_diagnostics_1",       plot_marginal_diagnostics(fit_brma_us, fit_brma_hcs))
  scenario_plot("robma_diagnostics_2",       plot_marginal_diagnostics(fit_brma_diag, fit_brma_hcs0))

  ### diagnostic plots ----
  set.seed(1)
  scenario_plot("funnel",  {
    par(mfrow = c(1, 2))
    funnel(fit_brma_us, main = "US")
    funnel(fit_brma_hcs, main = "HCS")
  })
  scenario_plot("qqnorm",  {
    par(mfrow = c(2, 2))
    qqnorm(fit_brma_us,     main = "US")
    qqnorm(fit_brma_hcs,    main = "HCS")
    qqnorm(fit_brma_diag, main = "DIAG")
    qqnorm(fit_brma_hcs0,   main = "HCS0")
  })

  scenario_plot("zplot",  zplot(fit_brma_hcs, to = 10, by.hist = 1))

})
