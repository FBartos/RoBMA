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
  fit_metafor_no_effect <- metafor::rma.mv(yi, V_assink, random = ~ 1 | study, data = dat.assink2016)
  fit_metafor_fixed     <- metafor::rma.mv(yi, V_assink, data = dat.assink2016)
  fit_metafor_diag_no_study   <- metafor::rma.mv(yi, V_assink_diagonal, random = ~ 1 | id, data = dat.assink2016)
  fit_metafor_diag_no_effect  <- metafor::rma.mv(yi, V_assink_diagonal, random = ~ 1 | study, data = dat.assink2016)
  fit_metafor_diag_fixed      <- metafor::rma.mv(yi, V_assink_diagonal, data = dat.assink2016)
  fit_metafor_reg      <- metafor::rma.mv(yi, V_assink,          mods = ~ deltype, random = ~ 1 | study / esid, data = dat.assink2016)
  fit_metafor_diag_reg <- metafor::rma.mv(yi, V_assink_diagonal, mods = ~ deltype, random = ~ 1 | study / esid, data = dat.assink2016)

  fit_brma.mv <- scenario_fit("fit_brma.mv", {
    tmp <- brma.mv(yi = yi, V = V_assink, measure = "SMD", random = ~ 1 | study / esid, data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_brma.mv_reg <- scenario_fit("fit_brma.mv_reg", {
    tmp <- brma.mv(yi = yi, V = V_assink, mods = ~ deltype, measure = "SMD", random = ~ 1 | study / esid, data = dat.assink2016, seed = 1)
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
  fit_brma.mv_fixed <- scenario_fit("fit_brma.mv_fixed", {
    tmp <- brma.mv(yi = yi, V = V_assink, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  # diagonal for comparison with brma.uni
  fit_brma.mv_diag <- scenario_fit("fit_brma.mv_diag", {
    tmp <- brma.mv(yi = yi, V = V_assink_diagonal, measure = "SMD", random = ~ 1 | study / esid, data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_brma.mv_diag_reg <- scenario_fit("fit_brma.mv_diag_reg", {
    tmp <- brma.mv(yi = yi, V = V_assink_diagonal, mods = ~ deltype, measure = "SMD", random = ~ 1 | study / esid, data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_brma.mv_diag_no_study <- scenario_fit("fit_brma.mv_diag_no_study", {
    tmp <- brma.mv(yi = yi, V = V_assink_diagonal, measure = "SMD", random = ~ 1 | study:esid, data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_brma.mv_diag_no_effect <- scenario_fit("fit_brma.mv_diag_no_effect", {
    tmp <- brma.mv(yi = yi, V = V_assink_diagonal, measure = "SMD", random = ~ 1 | study , data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_brma.mv_diag_fixed <- scenario_fit("fit_brma.mv_diag_fixed", {
    tmp <- brma.mv(yi = yi, V = V_assink_diagonal, measure = "SMD", prior_heterogeneity = NULL , data = dat.assink2016, seed = 1)
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
  fit_brma_cluster_reg <- scenario_fit("fit_brma_cluster_reg", {
    tmp <- brma(yi = yi, vi = vi, mods = ~ deltype, cluster = study, measure = "SMD", data = dat.assink2016, seed = 1)
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
  fit_brma_fixed <- scenario_fit("fit_brma_fixed", {
    tmp <- brma(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })

  ### model summary ----
  # assess equal models
  fit_metafor
  fit_metafor.cs
  scenario_text("summary-fit_brma.mv", summary(fit_brma.mv))

  fit_metafor_no_study
  scenario_text("summary-fit_brma.mv_no_study", summary(fit_brma.mv_no_study))

  fit_metafor_no_effect
  scenario_text("summary-fit_brma.mv_no_effect", summary(fit_brma.mv_no_effect))

  fit_metafor_fixed
  scenario_text("summary-fit_brma.mv_fixed", summary(fit_brma.mv_fixed))

  # univariate
  fit_metafor_diag
  fit_metafor_diag.cs
  scenario_text("summary-fit_brma.mv_diag", summary(fit_brma.mv_diag))
  scenario_text("summary-fit_brma_cluster", summary(fit_brma_cluster))

  fit_metafor_diag_no_study
  scenario_text("summary-fit_brma.mv_diag_no_study", summary(fit_brma.mv_diag_no_study))
  scenario_text("summary-fit_brma",                  summary(fit_brma))

  fit_metafor_diag_no_effect
  scenario_text("summary-fit_brma.mv_diag_no_effect", summary(fit_brma.mv_diag_no_effect))
  # cannot fit brma with study only random effects

  fit_metafor_diag_fixed
  scenario_text("summary-fit_brma.mv_diag_fixed", summary(fit_brma.mv_diag_fixed))
  scenario_text("summary-fit_brma_fixed",         summary(fit_brma_fixed))

  # meta-regression
  fit_metafor_reg # the estimates are shrunk -- the original ones are quite large
  scenario_text("summary-fit_brma.mv_reg", summary(fit_brma.mv_reg))

  fit_metafor_diag_reg # the estimates are shrunk -- the original ones are quite large
  scenario_text("summary-fit_brma.mv_diag_reg", summary(fit_brma.mv_diag_reg))
  scenario_text("summary-fit_brma_cluster_reg", summary(fit_brma_cluster_reg))

  # model comparison equivalence
  scenario_text("model-fit-equivalent", cbind.data.frame(
    "logml.mv"      = c(logml(fit_brma.mv),      logml(fit_brma.mv_no_study),      logml(fit_brma.mv_no_effect),      logml(fit_brma.mv_fixed)), # different from 2 and 3
    "logml.mv_diag" = c(logml(fit_brma.mv_diag), logml(fit_brma.mv_diag_no_study), logml(fit_brma.mv_diag_no_effect), logml(fit_brma.mv_diag_fixed)), # equal to 3
    "logml.brma"    = c(logml(fit_brma_cluster), logml(fit_brma),                  NA,                                logml((fit_brma_fixed)))
  ))
  getloo <- function(fit) loo(fit)[["estimates"]]["looic",1]
  scenario_text("model-fit-loo", cbind.data.frame(
    "getloo.mv"      = c(getloo(fit_brma.mv),      getloo(fit_brma.mv_no_study),      getloo(fit_brma.mv_no_effect),      getloo(fit_brma.mv_fixed)), # different from 2 and 3
    "getloo.mv_diag" = c(getloo(fit_brma.mv_diag), getloo(fit_brma.mv_diag_no_study), getloo(fit_brma.mv_diag_no_effect), getloo(fit_brma.mv_diag_fixed)), # equal to 3
    "getloo.brma"    = c(getloo(fit_brma_cluster), getloo(fit_brma),                  NA,                                 getloo((fit_brma_fixed)))
  ))

  # some additional checks
  scenario_text("comapre-loo",   loo_model_weights(fit_brma.mv, fit_brma.mv_no_study, fit_brma.mv_no_effect, fit_brma.mv_fixed))
  scenario_text("comapre-logml", t(t(round(post_prob(fit_brma.mv, fit_brma.mv_no_study, fit_brma.mv_no_effect, fit_brma.mv_fixed), 3))))

  scenario_text("comapre-diag-loo",   loo_model_weights(fit_brma.mv_diag, fit_brma.mv_diag_no_study, fit_brma.mv_diag_no_effect, fit_brma.mv_diag_fixed))
  scenario_text("comapre-diag-logml", t(t(round(post_prob(fit_brma.mv_diag, fit_brma.mv_diag_no_study, fit_brma.mv_diag_no_effect, fit_brma.mv_diag_fixed), 3))))

  scenario_text("comapre-uni-loo",   loo_model_weights(fit_brma_cluster, fit_brma, fit_brma_fixed))
  scenario_text("comapre-uni-logml", t(t(round(post_prob(fit_brma_cluster, fit_brma, fit_brma_fixed), 3))))

  # metareg comparion
  scenario_text("model-fit-reg-equivalent", cbind.data.frame(
    "logml.mv"      = logml(fit_brma.mv_reg), # different from 2 and 3
    "logml.mv_diag" = logml(fit_brma.mv_diag_reg), # equal to 3
    "logml.brma"    = logml(fit_brma_cluster_reg)
  ))


  ### basic fit plots ----
  set.seed(1)
  scenario_plot("fit.uni_posterior_rho", {
    plot(fit_brma_cluster, "rho", prior = TRUE)
    lines(fit_brma_cluster, "rho", density_method = "IWMDE", lty = 2)

    lines(fit_brma.mv_diag, "var_frac(random_total: study)", col = "blue")
    lines(fit_brma.mv_diag, "var_frac(random_total: study)", density_method = "qCMDE", col = "blue", lty = 2)
  })

  set.seed(1)
  scenario_plot("fit.mv_posterior_mu", {
    plot(fit_brma.mv, "mu", prior = TRUE, xlim = c(-1, 1))
    lines(fit_brma.mv, "mu", density_method = "IWMDE", lty = 2)

    # diagonal needs to be wider
    lines(fit_brma.mv_diag, "mu", col = "blue")
    lines(fit_brma.mv_diag, "mu", density_method = "qCMDE", col = "blue", lty = 2, density_control = list(samples = 2000))
  })

  set.seed(1)
  scenario_plot("fit.mv_posterior_mod", {
    plot(fit_brma.mv_reg, "deltype", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 3))
    lines(fit_brma.mv_reg, "deltype", density_method = "qCMDE", lty = 2, density_control = list(samples = 2000))
  })

  set.seed(1)
  scenario_plot("fit.mv_posterior_random", {
    par(mfrow = c(2, 3))

    plot(fit_brma.mv, "sd_total(random_total)", prior = TRUE)
    lines(fit_brma.mv, "sd_total(random_total)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "sd(intercept | study)", prior = TRUE)
    # FUTURE: would be nice to have but not essential
    # lines(fit_brma.mv, "sd(intercept | study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "sd(intercept | esid:study)", prior = TRUE)
    # FUTURE: would be nice to have but not essential
    # lines(fit_brma.mv, "sd(intercept | esid:study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "var_frac(random_total: esid_study)", prior = TRUE)
    lines(fit_brma.mv, "var_frac(random_total: esid_study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "var_frac(random_total: study)", prior = TRUE)
    lines(fit_brma.mv, "var_frac(random_total: study)", density_method = "qCMDE", lty = 2)
  })

  set.seed(1)
  scenario_plot("fit.mv_diag_posterior_random", {
    par(mfrow = c(2, 3))

    plot(fit_brma.mv_diag, "sd_total(random_total)", prior = TRUE)
    lines(fit_brma.mv_diag, "sd_total(random_total)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv_diag, "sd(intercept | study)", prior = TRUE)
    # FUTURE: would be nice to have but not essential
    # lines(fit_brma.mv, "sd(intercept | study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv_diag, "sd(intercept | esid:study)", prior = TRUE)
    # FUTURE: would be nice to have but not essential
    # lines(fit_brma.mv, "sd(intercept | esid:study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv_diag, "var_frac(random_total: esid_study)", prior = TRUE)
    lines(fit_brma.mv_diag, "var_frac(random_total: esid_study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv_diag, "var_frac(random_total: study)", prior = TRUE)
    lines(fit_brma.mv_diag, "var_frac(random_total: study)", density_method = "qCMDE", lty = 2)
  })


  ### hypothesis ----
  set.seed(1)
  BF_brma_rho    <- hypothesis(fit_brma_cluster, c("rho != 0 vs rho = 0", "rho != 1 vs rho = 1"), density_method = "qCMDE", density_control = list(samples = 2000))
  BF_mv_diag_rho <- hypothesis(fit_brma.mv_diag, c("var_frac(random_total: study) != 0 vs var_frac(random_total: study) = 0", "var_frac(random_total: study) != 1 vs var_frac(random_total: study) = 1"), density_method = "qCMDE", density_control = list(samples = 2000))
  scenario_text("fit_rho_bayes_factor_comparison", data.frame(
    rho                = c(0, 1),
    density_brma_BF    = BF_brma_rho[["BF"]],
    density_mv_diag_BF = BF_mv_diag_rho[["BF"]],
    marglik_brma_BF    = c(bf(fit_brma_cluster, fit_brma)$bf, NA),
    marglik_mv_diag_BF = c(bf(fit_brma.mv_diag, fit_brma.mv_diag_no_study)$bf, bf(fit_brma.mv_diag, fit_brma.mv_diag_no_effect)$bf)
  ))

  set.seed(1)
  BF_random      <- hypothesis(fit_brma.mv, c(
    "var_frac(random_total: study) != 0 vs var_frac(random_total: study) = 0", "var_frac(random_total: study) != 1 vs var_frac(random_total: study) = 1",
    "sd_total(random_total) = 0"
    ),density_method = "qCMDE", density_control = list(samples = 2000))
  scenario_text("fit_random_bayes_factor_comparison", data.frame(
    hypothesis = c("rho != 0", "rho != 1", "sd != 0"),
    density_BF = BF_random[["BF"]],
    marglik_BF = c(bf(fit_brma.mv, fit_brma.mv_no_study)$bf, bf(fit_brma.mv, fit_brma.mv_no_effect)$bf,
                   bf(fit_brma.mv, fit_brma.mv_fixed)$bf)
  ))

  set.seed(1)
  BF_mods <- hypothesis(fit_brma.mv_reg, c("deltype[general] = 0 vs deltype[general] != 0", "deltype[general] = 0 vs deltype[general] > 0", "deltype[general] > 0 vs deltype[general] < 0"),density_method = "qCMDE", density_control = list(samples = 2000))
  scenario_text("fit_mods", BF_mods)

  ### pooled effects ----
  compare_preds <- function(fit_metafor, fit_RoBMA, fit_RoBMA2 = NULL) {
    cbind.data.frame(
      "metafor"  = t(data.frame(predict(fit_metafor))[,-2]),
      "brma.mv"  = t(unname(data.frame(pooled_effect(fit_RoBMA))[,-2])),
      "brma.uni" = if (!is.null(fit_RoBMA2)) t(unname(data.frame(pooled_effect(fit_RoBMA2))[,-2])) else rep(NA, 5))
  }

  set.seed(1)
  scenario_text("pooled-effect-1",  compare_preds(fit_metafor,           fit_brma.mv))
  scenario_text("pooled-effect-2",  compare_preds(fit_metafor_no_effect, fit_brma.mv_no_effect))
  scenario_text("pooled-effect-3",  compare_preds(fit_metafor_no_study,  fit_brma.mv_no_study))
  scenario_text("pooled-effect-4",  compare_preds(fit_metafor_fixed,     fit_brma.mv_fixed))

  scenario_text("pooled-effect-1-diag",  compare_preds(fit_metafor_diag,           fit_brma.mv_diag,           fit_brma_cluster))
  scenario_text("pooled-effect-2-diag",  compare_preds(fit_metafor_diag_no_effect, fit_brma.mv_diag_no_effect))
  scenario_text("pooled-effect-3-diag",  compare_preds(fit_metafor_diag_no_study,  fit_brma.mv_diag_no_study,  fit_brma))
  scenario_text("pooled-effect-4-diag",  compare_preds(fit_metafor_diag_fixed,     fit_brma.mv_diag_fixed,     fit_brma_fixed))

  ### predictions ----
  set.seed(1)
  compare_preds_reg <- function(fit_metafor, fit_RoBMA, fit_RoBMA2 = NULL) {
    cbind.data.frame(
      "metafor"  = unlist(data.frame(predict(fit_metafor))[c(1, 50, 80),-c(2, 5, 6)]),
      "brma.mv"  = unlist(unname(data.frame(predict(fit_RoBMA))[c(1, 50, 80),-2])),
      "brma.uni" = if (!is.null(fit_RoBMA2)) unlist(unname(data.frame(predict(fit_RoBMA2))[c(1, 50, 80),-2])) else rep(NA, 9))
  }
  compare_preds_reg_pi <- function(fit_metafor, fit_RoBMA, fit_RoBMA2 = NULL) {
    cbind.data.frame(
      "metafor"  = unlist(data.frame(predict(fit_metafor))[c(1, 50, 80),-c(2, 3, 4)]),
      "brma.mv"  = unlist(unname(data.frame(predict(fit_RoBMA, type = "estimate"))[c(1, 50, 80),-2])),
      "brma.uni" = if (!is.null(fit_RoBMA2)) unlist(unname(data.frame(predict(fit_RoBMA2, type = "estimate"))[c(1, 50, 80),-2])) else rep(NA, 9))
  }

  scenario_text("pooled-effect-reg",          compare_preds_reg(fit_metafor_reg,      fit_brma.mv_reg))
  scenario_text("pooled-effect-diag-reg",     compare_preds_reg(fit_metafor_diag_reg, fit_brma.mv_diag_reg, fit_brma_cluster_reg))
  scenario_text("pooled-effect-reg_pi",       compare_preds_reg_pi(fit_metafor_reg,      fit_brma.mv_reg)) # TODO examine
  scenario_text("pooled-effect-diag-reg_pi",  compare_preds_reg_pi(fit_metafor_diag_reg, fit_brma.mv_diag_reg, fit_brma_cluster_reg))

  ### marginal means ----
  set.seed(1)
  scenario_text("marginal_means", marginal_means(fit_brma.mv_reg))
  scenario_plot("marginal_means_plot", plot(marginal_means(fit_brma.mv_reg), "deltype", prior = TRUE, xlim = c(-2, 2)))

  ### summary heterogeneity ----
  set.seed(1)
  scenario_text("summary_heterogeneity-1a", summary_heterogeneity(fit_brma.mv_diag))
  scenario_text("summary_heterogeneity-1b", summary_heterogeneity(fit_brma_cluster))

  scenario_text("summary_heterogeneity-2", summary_heterogeneity(fit_brma.mv))
  scenario_text("summary_heterogeneity-3", summary_heterogeneity(fit_brma.mv_no_effect))
  scenario_text("summary_heterogeneity-4", summary_heterogeneity(fit_brma.mv_fixed))

  ### random effects ----
  ranef_metafor <- metafor::ranef(fit_metafor)
  ranef_brma.mv <- ranef(fit_brma.mv)

  ranef_metafor_diag <- metafor::ranef(fit_metafor_diag)
  ranef_brma.mv_diag <- ranef(fit_brma.mv_diag)
  ranef_brma         <- ranef(fit_brma_cluster)

  scenario_plot("ranef_mv", {
    par(mfrow = c(1, 2))
    scenario_agreement_plot(ranef_metafor$study[["intrcpt"]], as.data.frame(ranef_brma.mv$study)[["Mean"]], main = "study")
    scenario_agreement_plot(ranef_metafor$`study/esid`[["intrcpt"]], as.data.frame(ranef_brma.mv$esid_study)[["Mean"]], main = "esid_study")
  })
  scenario_plot("ranef_mv_diag", {
    par(mfrow = c(2, 2))
    scenario_agreement_plot(ranef_metafor_diag$study[["intrcpt"]], as.data.frame(ranef_brma.mv_diag$study)[["Mean"]], main = "study")
    scenario_agreement_plot(ranef_metafor_diag$`study/esid`[["intrcpt"]], as.data.frame(ranef_brma.mv_diag$esid_study)[["Mean"]], main = "esid_study")

    scenario_agreement_plot(ranef_metafor_diag$study[["intrcpt"]], as.data.frame(ranef_brma$cluster)[["Mean"]], main = "study")
    scenario_agreement_plot(ranef_metafor_diag$`study/esid`[["intrcpt"]], as.data.frame(ranef_brma$estimate)[["Mean"]], main = "esid_study")
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

  scenario_plot("fit_mv_marginal_diagnostics",      {plot_marginal_diagnostics(fit_metafor, fit_brma.mv)})
  scenario_plot("fit_mv_diag_marginal_diagnostics", {plot_marginal_diagnostics(fit_metafor_diag, fit_brma.mv_diag)})
  scenario_plot("fit_mv_reg_marginal_diagnostics",  {plot_marginal_diagnostics(fit_metafor_reg, fit_brma.mv_reg)})

  ### diagnostic plots ----
  set.seed(1)
  scenario_plot("funnel_brma_mv",  {
    par(mfrow = c(1, 2))
    funnel(fit_brma.mv, main = "funnel")
    bfunnel(fit_brma.mv, main = "bfunnel")
  })
  scenario_plot("qqnorm_brma_mv", qqnorm(fit_brma.mv, main = "RoBMA"))
  scenario_plot("zplot_brma_mv",  zplot(fit_brma.mv, to = 10))

  set.seed(1)
  scenario_plot("funnel_brma_mv_reg", funnel(fit_brma.mv_reg, main = "funnel"))
  scenario_plot("qqnorm_brma_mv_reg", qqnorm(fit_brma.mv_reg))
  scenario_plot("zplot_brma_mv_reg",  zplot(fit_brma.mv_reg, to = 10))

})
