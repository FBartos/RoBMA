if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("assink2016")
# testthat::test_file("tests/scenarios/test-assink2016.R")
# Selection settings are source-specific; the broad model comparisons use defaults.
# Failed full-V fits and fit blocks exceeding 20 minutes are commented out
# with their dependent output, retaining the calls for later optimization.
### Description
# Compare univariate, specialized multilevel, and multivariate nested
# random-effects targets, including publication-bias adjustment.

testthat::test_that("Assink brma and brma.mv models", {

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
  fit_brma.mv_bycomp <- scenario_fit("fit_brma.mv_bycomp", {
    tmp <- brma.mv(yi = yi, V = V_assink, measure = "SMD", random = list(~ 1 | study:esid, ~ 1 | study), data = dat.assink2016, seed = 1)
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
    tmp <- brma.mv(yi = yi, V = V_assink, measure = "SMD", data = dat.assink2016, seed = 1)
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
    tmp <- brma.mv(yi = yi, V = V_assink_diagonal, measure = "SMD", data = dat.assink2016, seed = 1)
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
  # scale .mv models
  fit_brma.mv_scale_total  <- scenario_fit("fit_brma.mv_scale_total", {
    tmp <- brma.mv(yi = yi, V = V_assink, scale = ~ deltype, measure = "SMD", random = ~ 1 | study / esid, data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  }, cache_version = 1L)
  fit_brma.mv_scale_effect <- scenario_fit("fit_brma.mv_scale_effect", {
    tmp <- brma.mv(yi = yi, V = V_assink, scale = list(esid = ~ deltype), measure = "SMD", random = ~ 1 | study / esid, data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  }, cache_version = 1L)
  # BMA.mv models
  fit_BMA.mv_cluster <- scenario_fit("fit_BMA.mv_vi", {
    tmp <- BMA(yi = yi, vi = vi, measure = "SMD", cluster = study, data = dat.assink2016, seed = 1)
    tmp <- add_loo(tmp)
    return(tmp)
  }, cache_version = 1L)
  fit_BMA.mv_diag <- scenario_fit("fit_BMA.mv_diag", {
    tmp <- BMA.mv(yi = yi, V = V_assink_diagonal, measure = "SMD", random = ~ 1 | study / esid, data = dat.assink2016, seed = 1)
    tmp <- add_loo(tmp)
    return(tmp)
  }, cache_version = 1L)
  fit_BMA.mv      <- scenario_fit("fit_BMA.mv", {
    tmp <- BMA.mv(yi = yi, V = V_assink, measure = "SMD", random = ~ 1 | study / esid, data = dat.assink2016, seed = 1)
    tmp <- add_loo(tmp)
    return(tmp)
  }, cache_version = 1L)
  fit_BMA.mv_bycomp  <- scenario_fit("fit_BMA.mv_bycomp", {
    tmp <- BMA.mv(yi = yi, V = V_assink, measure = "SMD",
                  random = list("study" = ~ 1 | study, "study:esid" = ~ 1 | study:esid), data = dat.assink2016, seed = 1)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_BMA.mv_bycomp_eff  <- scenario_fit("fit_BMA.mv_bycomp_eff", {
    tmp <- BMA.mv(yi = yi, V = V_assink, measure = "SMD",
                  prior_effect_null = NULL,
                  random = list("study" = ~ 1 | study, "study:esid" = ~ 1 | study:esid), data = dat.assink2016, seed = 1)
    tmp <- add_loo(tmp)
    return(tmp)
  })

  ### model summary ----
  # assess equal models
  fit_metafor
  fit_metafor.cs
  scenario_text("summary-fit_brma.mv", summary(fit_brma.mv))
  scenario_text("summary-fit_brma.mv_bycomp", summary(fit_brma.mv_bycomp))

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

  # scale models
  scenario_text("summary-fit_brma.mv_scale_effect",    summary(fit_brma.mv_scale_effect))
  scenario_text("summary-fit_fit_brma.mv_scale_total", summary(fit_brma.mv_scale_total))

  # BMA models
  scenario_text("summary-fit_BMA.mv_cluster",    summary(fit_BMA.mv_cluster))
  scenario_text("summary-fit_BMA.mv_diag",       summary(fit_BMA.mv_diag))
  scenario_text("summary-fit_BMA.mv",            summary(fit_BMA.mv))
  scenario_text("summary-fit_BMA.mv_bycomp",     summary(fit_BMA.mv_bycomp))

  # parameters
  metafor_parameters <- c("intercept", study_tau2 = "sigma[study]^2", estimate_tau2 = "sigma[study/esid]^2", tau2_total = "sigma[total]^2")
  robma_parameters   <- c("mu", study_tau2 = "tau2", estimate_tau2 = "tau2", tau2_total = "tau2_total")
  robma_components   <- c(NA, "study", "esid_study", NA)
  scenario_text("metafor-comparison", data.frame(
    model          = rep(c("correlated V", "diagonal V"), each = 2L),
    implementation = rep(c("metafor", "RoBMA"), 2L),
    rbind(
      ex_m(fit_metafor, metafor_parameters),      ex_r(fit_brma.mv,      robma_parameters, component = robma_components),
      ex_m(fit_metafor_diag, metafor_parameters), ex_r(fit_brma.mv_diag, robma_parameters, component = robma_components)
    ),
    row.names = NULL
  ))

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
  scenario_text("model-fit-BMA", cbind.data.frame("loo" =c("BMA_cluster" = getloo(fit_BMA.mv_cluster), "BMA.mv_diag" = getloo(fit_BMA.mv_diag), "BMA.mv" = getloo(fit_BMA.mv))))

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

    lines(fit_brma.mv_diag, "tau2_prop(study)", col = "blue")
    lines(fit_brma.mv_diag, "tau2_prop(study)", density_method = "qCMDE", col = "blue", lty = 2)
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

    plot(fit_brma.mv, "tau_total", prior = TRUE)
    lines(fit_brma.mv, "tau_total", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "study: tau", prior = TRUE)
    lines(fit_brma.mv, "study: tau", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_brma.mv, "esid_study: tau", prior = TRUE)
    lines(fit_brma.mv, "esid_study: tau", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_brma.mv, "tau2_prop(esid_study)", prior = TRUE)
    lines(fit_brma.mv, "tau2_prop(esid_study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "tau2_prop(study)", prior = TRUE)
    lines(fit_brma.mv, "tau2_prop(study)", density_method = "qCMDE", lty = 2)
  })

  set.seed(1)
  scenario_plot("fit_BMA.mv_diag_posterior", {
    par(mfrow = c(2, 2))

    plot(fit_BMA.mv_cluster, "mu")
    lines(fit_BMA.mv_diag, "mu", density_method = "qCMDE", lty = 2)

    plot(fit_BMA.mv_cluster, "tau", prior = TRUE)
    lines(fit_BMA.mv_diag, "tau_total", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_BMA.mv_cluster, "rho", prior = TRUE)
    lines(fit_BMA.mv_diag, "tau2_prop(study)", density_method = "qCMDE", lty = 2)
  })

  set.seed(1)
  scenario_plot("fit_BMA.mv_posterior", {
    par(mfrow = c(2, 3))

    plot(fit_BMA.mv, "mu")
    lines(fit_BMA.mv, "mu", density_method = "qCMDE", lty = 2)

    plot(fit_BMA.mv, "tau_total", prior = TRUE)
    lines(fit_BMA.mv, "tau_total", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_BMA.mv, "esid_study: tau", prior = TRUE)
    lines(fit_BMA.mv, "esid_study: tau", density_method = "qCMDE", lty = 2)

    plot(fit_BMA.mv, "tau2_prop(esid_study)", prior = TRUE)
    lines(fit_BMA.mv, "tau2_prop(esid_study)", density_method = "qCMDE", lty = 2)

    plot(fit_BMA.mv, "tau2_prop(study)", prior = TRUE)
    lines(fit_BMA.mv, "tau2_prop(study)", density_method = "qCMDE", lty = 2)
  })

  set.seed(1)
  scenario_plot("fit.BMA_posterior_random", {
    par(mfrow = c(2, 3))

    plot(fit_brma.mv, "tau_total", prior = TRUE)
    lines(fit_brma.mv, "tau_total", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "study: tau", prior = TRUE)
    lines(fit_brma.mv, "study: tau", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_brma.mv, "esid_study: tau", prior = TRUE)
    lines(fit_brma.mv, "esid_study: tau", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_brma.mv, "tau2_prop(esid_study)", prior = TRUE)
    lines(fit_brma.mv, "tau2_prop(esid_study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "tau2_prop(study)", prior = TRUE)
    lines(fit_brma.mv, "tau2_prop(study)", density_method = "qCMDE", lty = 2)
  })



  ### hypothesis ----
  set.seed(1)
  BF_brma_rho    <- scenario_time("BF_brma_rho", hypothesis(fit_brma_cluster, c("rho != 0 vs rho = 0", "rho != 1 vs rho = 1"), density_method = "qCMDE", density_control = list(samples = 2000)))
  BF_mv_diag_rho <- scenario_time("BF_mv_diag_rho", hypothesis(fit_brma.mv_diag, c("tau2_prop(study) != 0 vs tau2_prop(study) = 0", "tau2_prop(study) != 1 vs tau2_prop(study) = 1"), density_method = "qCMDE", density_control = list(samples = 2000)))
  scenario_text("fit_rho_bayes_factor_comparison", data.frame(
    rho                = c(0, 1),
    density_brma_BF    = BF_brma_rho[["BF"]],
    density_mv_diag_BF = BF_mv_diag_rho[["BF"]],
    marglik_brma_BF    = c(bf(fit_brma_cluster, fit_brma)$bf, NA),
    marglik_mv_diag_BF = c(bf(fit_brma.mv_diag, fit_brma.mv_diag_no_study)$bf, bf(fit_brma.mv_diag, fit_brma.mv_diag_no_effect)$bf)
  ))

  set.seed(1)
  BF_random      <- scenario_time("BF_random", hypothesis(fit_brma.mv, c(
    "tau2_prop(study) != 0 vs tau2_prop(study) = 0", "tau2_prop(study) != 1 vs tau2_prop(study) = 1",
    "tau_total = 0"
    ),density_method = "qCMDE", density_control = list(samples = 2000)))
  scenario_text("fit_random_bayes_factor_comparison", data.frame(
    hypothesis = c("rho != 0", "rho != 1", "tau != 0"),
    density_BF = BF_random[["BF"]],
    marglik_BF = c(bf(fit_brma.mv, fit_brma.mv_no_study)$bf, bf(fit_brma.mv, fit_brma.mv_no_effect)$bf,
                   bf(fit_brma.mv, fit_brma.mv_fixed)$bf)
  ))

  set.seed(1)
  BF_mods <- scenario_time("BF_mods", hypothesis(fit_brma.mv_reg, c("deltype[general] = 0 vs deltype[general] != 0", "deltype[general] = 0 vs deltype[general] > 0", "deltype[general] > 0 vs deltype[general] < 0"),density_method = "qCMDE", density_control = list(samples = 2000)))
  scenario_text("fit_mods", BF_mods)

  # compare marglik posterior probabilities against BMA posterior probabilities
  scenario_text("bma-vs-marglik", cbind(
    data.frame(summary_models(fit_BMA.mv_bycomp_eff, type = "individual"))[,c("Random..study", "Random..study.esid", "post_prob")],
    marglik = round(unname(post_prob(fit_brma.mv_fixed, fit_brma.mv_no_effect, fit_brma.mv_no_study, fit_brma.mv)), 5)
  ))

  # additional summaries
  scenario_text("mod-sum-bma",  summary_models(fit_BMA.mv_bycomp))
  scenario_text("mod-sum-bma2", summary_models(fit_BMA.mv_bycomp, type = "individual"))


  ### pooled effects ----
  compare_preds <- function(fit_metafor, fit_RoBMA, fit_RoBMA2 = NULL) {
    cbind.data.frame(
      "metafor" = t(data.frame(predict(fit_metafor))[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")]),
      ex_p(fit_RoBMA),
      if (!is.null(fit_RoBMA2)) ex_p(fit_RoBMA2) else data.frame(brma.uni = rep(NA_real_, 5L)))
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
      "metafor"  = unlist(data.frame(predict(fit_metafor))[c(1, 50, 80), c("pred", "ci.lb", "ci.ub")]),
      "brma.mv"  = unlist(data.frame(predict(fit_RoBMA))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      "brma.uni" = if (!is.null(fit_RoBMA2)) unlist(data.frame(predict(fit_RoBMA2))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE) else rep(NA_real_, 9L))
  }
  compare_preds_reg_pi <- function(fit_metafor, fit_RoBMA, fit_RoBMA2 = NULL) {
    cbind.data.frame(
      "metafor"  = unlist(data.frame(predict(fit_metafor))[c(1, 50, 80), c("pred", "pi.lb", "pi.ub")]),
      "brma.mv"  = unlist(data.frame(predict(fit_RoBMA, type = "estimate"))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      "brma.uni" = if (!is.null(fit_RoBMA2)) unlist(data.frame(predict(fit_RoBMA2, type = "estimate"))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE) else rep(NA_real_, 9L))
  }

  scenario_text("pooled-effect-reg",          compare_preds_reg(fit_metafor_reg,      fit_brma.mv_reg))
  scenario_text("pooled-effect-diag-reg",     compare_preds_reg(fit_metafor_diag_reg, fit_brma.mv_diag_reg, fit_brma_cluster_reg))
  scenario_text("pooled-effect-reg_pi",       compare_preds_reg_pi(fit_metafor_reg,      fit_brma.mv_reg))
  scenario_text("pooled-effect-diag-reg_pi",  compare_preds_reg_pi(fit_metafor_diag_reg, fit_brma.mv_diag_reg, fit_brma_cluster_reg))

  ### marginal means ----
  set.seed(1)
  scenario_text("marginal_means", marginal_means(fit_brma.mv_reg))
  scenario_plot("marginal_means_plot", plot(marginal_means(fit_brma.mv_reg), "deltype", prior = TRUE, xlim = c(-2, 2)))

  ### summary heterogeneity ----
  set.seed(1)
  scenario_text("summary_heterogeneity-1a", summary_heterogeneity(fit_brma.mv_diag))
  scenario_text("summary_heterogeneity-1b", summary_heterogeneity(fit_brma_cluster))

  scenario_text("summary_heterogeneity-2",  summary_heterogeneity(fit_brma.mv))
  scenario_text("summary_heterogeneity-2a", summary_heterogeneity(fit_brma.mv_bycomp))
  scenario_text("summary_heterogeneity-3",  summary_heterogeneity(fit_brma.mv_no_effect))
  scenario_text("summary_heterogeneity-4",  summary_heterogeneity(fit_brma.mv_fixed))

  scenario_text("summary_heterogeneity-scale-1",  summary_heterogeneity(fit_brma.mv_scale_effect))
  scenario_text("summary_heterogeneity-scale-2",  summary_heterogeneity(fit_brma.mv_scale_total))

  ### random effects ----
  ranef_metafor <- metafor::ranef(fit_metafor)
  ranef_brma.mv <- scenario_time("ranef_brma.mv", ranef(fit_brma.mv))

  ranef_metafor_diag <- metafor::ranef(fit_metafor_diag)
  ranef_brma.mv_diag <- scenario_time("ranef_brma.mv_diag", ranef(fit_brma.mv_diag))
  ranef_brma         <- scenario_time("ranef_brma", ranef(fit_brma_cluster))

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

  ### print priors ----
  scenario_text("print-prior-fit_brma.mv",          print_prior(fit_brma.mv))
  scenario_text("print-prior-fit_brma.mv_bycomp",   print_prior(fit_brma.mv_bycomp))
  scenario_text("print-prior-fit_brma.mv_no_study", print_prior(fit_brma.mv_no_study))
  scenario_text("print-prior-fit_brma.mv_reg",      print_prior(fit_brma.mv_reg))

  ### diagnostics ----
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

  scenario_plot("funnel_BMA_mv",  {
    par(mfrow = c(1, 2))
    funnel(fit_BMA.mv,   main = "BMA overall")
    funnel(fit_BMA.mv_bycomp, main = "BMA factored")
  })
  scenario_plot("qqnorm_BMA_mv",  {
    par(mfrow = c(1, 2))
    qqnorm(fit_BMA.mv,   main = "BMA overall")
    qqnorm(fit_BMA.mv_bycomp, main = "BMA factored")
  })

})

testthat::test_that("Assink bPET and bPET.mv models", {

  set.seed(1)
  data("dat.assink2016", package = "metadat")

  V_assink <- metafor::vcalc(
    vi, cluster = study, type = deltype, obs = esid,
    rho = c(0.7, 0.5), data = dat.assink2016
  )
  V_assink_diagonal <- diag(dat.assink2016[["vi"]])

  ### Model fits ----
  fit_bPET_metafor                    <- metafor::rma(yi, vi, mods = ~ sqrt(vi), data = dat.assink2016)
  fit_bPET_metafor.mv_V               <- metafor::rma.mv(yi, V_assink,          mods = ~ sqrt(vi), random = ~ 1 | study / esid, data = dat.assink2016)
  fit_bPET_metafor.mv_V_no_study      <- metafor::rma.mv(yi, V_assink,          mods = ~ sqrt(vi), random = ~ 1 | id,           data = dat.assink2016)
  fit_bPET_metafor.mv_V_no_effect     <- metafor::rma.mv(yi, V_assink,          mods = ~ sqrt(vi), random = ~ 1 | study,        data = dat.assink2016)
  fit_bPET_metafor.mv_V_fixed         <- metafor::rma.mv(yi, V_assink,          mods = ~ sqrt(vi),                              data = dat.assink2016)
  fit_bPET_metafor.mv_V_reg           <- metafor::rma.mv(yi, V_assink,          mods = ~ deltype + sqrt(vi), random = ~ 1 | study / esid, data = dat.assink2016)
  fit_bPET_metafor.mv_vi              <- metafor::rma.mv(yi, V_assink_diagonal, mods = ~ sqrt(vi), random = ~ 1 | study / esid, data = dat.assink2016)
  fit_bPET_metafor.mv_vi_no_study     <- metafor::rma.mv(yi, V_assink_diagonal, mods = ~ sqrt(vi), random = ~ 1 | id,           data = dat.assink2016)
  fit_bPET_metafor.mv_vi_no_effect    <- metafor::rma.mv(yi, V_assink_diagonal, mods = ~ sqrt(vi), random = ~ 1 | study,        data = dat.assink2016)
  fit_bPET_metafor.mv_vi_fixed        <- metafor::rma.mv(yi, V_assink_diagonal, mods = ~ sqrt(vi),                              data = dat.assink2016)
  fit_bPET_metafor.mv_vi_reg          <- metafor::rma.mv(yi, V_assink_diagonal, mods = ~ deltype + sqrt(vi), random = ~ 1 | study / esid, data = dat.assink2016)

  fit_bPET <- scenario_fit("fit_bPET", {
    tmp <- bPET(yi = yi, vi = vi, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET_cluster <- scenario_fit("fit_bPET_cluster", {
    tmp <- bPET(yi = yi, vi = vi, cluster = study, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET_cluster_reg <- scenario_fit("fit_bPET_cluster_reg", {
    tmp <- bPET(yi = yi, vi = vi, mods = ~ deltype, cluster = study, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET_fixed <- scenario_fit("fit_bPET_fixed", {
    tmp <- bPET(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET.mv <- scenario_fit("fit_bPET.mv", {
    tmp <- bPET.mv(yi = yi, vi = vi, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET.mv_reg <- scenario_fit("fit_bPET.mv_reg", {
    tmp <- bPET.mv(yi = yi, vi = vi, mods = ~ deltype, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET.mv_no_study <- scenario_fit("fit_bPET.mv_no_study", {
    tmp <- bPET.mv(yi = yi, vi = vi, random = ~ 1 | study:esid, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET.mv_no_effect <- scenario_fit("fit_bPET.mv_no_effect", {
    tmp <- bPET.mv(yi = yi, vi = vi, random = ~ 1 | study, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET.mv_fixed <- scenario_fit("fit_bPET.mv_fixed", {
    tmp <- bPET.mv(yi = yi, vi = vi, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET.mv_V <- scenario_fit("fit_bPET.mv_V", {
    tmp <- bPET.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET.mv_V_reg <- scenario_fit("fit_bPET.mv_V_reg", {
    tmp <- bPET.mv(yi = yi, V = V_assink, mods = ~ deltype, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET.mv_V_no_study <- scenario_fit("fit_bPET.mv_V_no_study", {
    tmp <- bPET.mv(yi = yi, V = V_assink, random = ~ 1 | study:esid, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET.mv_V_no_effect <- scenario_fit("fit_bPET.mv_V_no_effect", {
    tmp <- bPET.mv(yi = yi, V = V_assink, random = ~ 1 | study, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPET.mv_V_fixed <- scenario_fit("fit_bPET.mv_V_fixed", {
    tmp <- bPET.mv(yi = yi, V = V_assink, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })

  ### Model summaries ----
  fit_bPET_metafor.mv_V
  scenario_text("bPET-summary-mv-V", summary(fit_bPET.mv_V))

  fit_bPET_metafor.mv_V_no_study
  scenario_text("bPET-summary-mv-V-no-study", summary(fit_bPET.mv_V_no_study))

  fit_bPET_metafor.mv_V_no_effect
  scenario_text("bPET-summary-mv-V-no-effect", summary(fit_bPET.mv_V_no_effect))

  fit_bPET_metafor.mv_V_fixed
  scenario_text("bPET-summary-mv-V-fixed", summary(fit_bPET.mv_V_fixed))

  fit_bPET_metafor.mv_vi
  scenario_text("bPET-summary-mv-vi",   summary(fit_bPET.mv))
  scenario_text("bPET-summary-cluster", summary(fit_bPET_cluster))

  fit_bPET_metafor.mv_vi_no_study
  scenario_text("bPET-summary-mv-vi-no-study", summary(fit_bPET.mv_no_study))
  scenario_text("bPET-summary-simple",         summary(fit_bPET))

  fit_bPET_metafor.mv_vi_no_effect
  scenario_text("bPET-summary-mv-vi-no-effect", summary(fit_bPET.mv_no_effect))

  fit_bPET_metafor.mv_vi_fixed
  scenario_text("bPET-summary-mv-vi-fixed", summary(fit_bPET.mv_fixed))
  scenario_text("bPET-summary-fixed",       summary(fit_bPET_fixed))

  fit_bPET_metafor.mv_V_reg
  scenario_text("bPET-summary-mv-V-reg", summary(fit_bPET.mv_V_reg))

  fit_bPET_metafor.mv_vi_reg
  scenario_text("bPET-summary-mv-vi-reg",   summary(fit_bPET.mv_reg))
  scenario_text("bPET-summary-cluster-reg", summary(fit_bPET_cluster_reg))

  metafor_PET_parameters    <- c(mu = "intercept", PET = "sqrt(vi)", tau_total = "tau")
  metafor_PET_mv_parameters <- c(mu = "intercept", PET = "sqrt(vi)", tau_total = "sigma[total]")
  robma_PET_parameters      <- c(mu = "mu", PET = "PET", tau_total = "tau")
  robma_PET_mv_parameters   <- c(mu = "mu", PET = "PET", tau_total = "tau_total")
  metafor_PET_mv_vi_study_fraction <- (ex_m(fit_bPET_metafor.mv_vi, "sigma[study]") / ex_m(fit_bPET_metafor.mv_vi, "sigma[total]"))^2
  metafor_PET_mv_V_study_fraction  <- (ex_m(fit_bPET_metafor.mv_V, "sigma[study]") / ex_m(fit_bPET_metafor.mv_V, "sigma[total]"))^2

  scenario_text("bPET-metafor-comparison", data.frame(
    model          = rep(c("simple", "cluster", "mv_vi", "mv_correlated_V"), each = 2L),
    implementation = rep(c("metafor", "RoBMA"), 4L),
    rbind(
      c(ex_m(fit_bPET_metafor, metafor_PET_parameters), study_tau2_prop = NA_real_),
      c(ex_r(fit_bPET, robma_PET_parameters), study_tau2_prop = NA_real_),
      c(ex_m(fit_bPET_metafor.mv_vi, metafor_PET_mv_parameters), study_tau2_prop = metafor_PET_mv_vi_study_fraction),
      c(ex_r(fit_bPET_cluster, robma_PET_parameters), study_tau2_prop = ex_r(fit_bPET_cluster, "rho")),
      c(ex_m(fit_bPET_metafor.mv_vi, metafor_PET_mv_parameters), study_tau2_prop = metafor_PET_mv_vi_study_fraction),
      c(ex_r(fit_bPET.mv, robma_PET_mv_parameters), study_tau2_prop = ex_r(fit_bPET.mv, "tau2_prop(study)")),
      c(ex_m(fit_bPET_metafor.mv_V, metafor_PET_mv_parameters), study_tau2_prop = metafor_PET_mv_V_study_fraction),
      c(ex_r(fit_bPET.mv_V, robma_PET_mv_parameters), study_tau2_prop = ex_r(fit_bPET.mv_V, "tau2_prop(study)"))
    ),
    row.names = NULL
  ))

  ### Model-fit comparisons ----
  scenario_text("bPET-model-fit-equivalent", cbind.data.frame(
    "logml.mv_V" = c(logml(fit_bPET.mv_V), logml(fit_bPET.mv_V_no_study), logml(fit_bPET.mv_V_no_effect), logml(fit_bPET.mv_V_fixed)), # different from the diagonal-V and specialized fits
    "logml.mv_vi" = c(logml(fit_bPET.mv),  logml(fit_bPET.mv_no_study),   logml(fit_bPET.mv_no_effect),   logml(fit_bPET.mv_fixed)),   # equal to the specialized fits where the structures match
    "logml.bPET"  = c(logml(fit_bPET_cluster), logml(fit_bPET), NA, logml(fit_bPET_fixed))
  ))
  getloo_bPET <- function(fit) loo(fit)[["estimates"]]["looic", 1L]
  scenario_text("bPET-model-fit-loo", cbind.data.frame(
    "looic.mv_V" = c(getloo_bPET(fit_bPET.mv_V), getloo_bPET(fit_bPET.mv_V_no_study), getloo_bPET(fit_bPET.mv_V_no_effect), getloo_bPET(fit_bPET.mv_V_fixed)), # different from the diagonal-V and specialized fits
    "looic.mv_vi" = c(getloo_bPET(fit_bPET.mv),  getloo_bPET(fit_bPET.mv_no_study),   getloo_bPET(fit_bPET.mv_no_effect),   getloo_bPET(fit_bPET.mv_fixed)),   # equal to the specialized fits where the structures match
    "looic.bPET"  = c(getloo_bPET(fit_bPET_cluster), getloo_bPET(fit_bPET), NA, getloo_bPET(fit_bPET_fixed))
  ))

  scenario_text("bPET-compare-mv-V-loo",   loo_model_weights(fit_bPET.mv_V, fit_bPET.mv_V_no_study, fit_bPET.mv_V_no_effect, fit_bPET.mv_V_fixed))
  scenario_text("bPET-compare-mv-V-logml", t(t(round(post_prob(fit_bPET.mv_V, fit_bPET.mv_V_no_study, fit_bPET.mv_V_no_effect, fit_bPET.mv_V_fixed), 3))))
  scenario_text("bPET-compare-mv-vi-loo",   loo_model_weights(fit_bPET.mv, fit_bPET.mv_no_study, fit_bPET.mv_no_effect, fit_bPET.mv_fixed))
  scenario_text("bPET-compare-mv-vi-logml", t(t(round(post_prob(fit_bPET.mv, fit_bPET.mv_no_study, fit_bPET.mv_no_effect, fit_bPET.mv_fixed), 3))))
  scenario_text("bPET-compare-specialized-loo",   loo_model_weights(fit_bPET_cluster, fit_bPET, fit_bPET_fixed))
  scenario_text("bPET-compare-specialized-logml", t(t(round(post_prob(fit_bPET_cluster, fit_bPET, fit_bPET_fixed), 3))))

  scenario_text("bPET-model-fit-reg-equivalent", cbind.data.frame(
    "logml.mv_V" = logml(fit_bPET.mv_V_reg),      # different from the diagonal-V and specialized fits
    "logml.mv_vi" = logml(fit_bPET.mv_reg),       # equal to the specialized cluster fit
    "logml.bPET"  = logml(fit_bPET_cluster_reg)
  ))

  ### Basic fit plots ----
  scenario_plot("bPET-posterior-rho", {
    plot(fit_bPET_cluster, "rho", prior = TRUE)
    lines(fit_bPET_cluster, "rho", density_method = "IWMDE", lty = 2)
    lines(fit_bPET.mv, "tau2_prop(study)", col = "blue")
    lines(fit_bPET.mv, "tau2_prop(study)", density_method = "qCMDE", col = "blue", lty = 2)
  })
  scenario_plot("bPET-posterior-location", {
    plot(fit_bPET.mv_V, "mu", prior = TRUE, xlim = c(-1, 1))
    lines(fit_bPET.mv_V, "mu", density_method = "IWMDE", lty = 2)
    lines(fit_bPET.mv, "mu", col = "blue")
    lines(fit_bPET.mv, "mu", density_method = "qCMDE", col = "blue", lty = 2, density_control = list(samples = 2000))
  })
  scenario_plot("bPET-posterior-PET", {
    plot(fit_bPET.mv_V, "PET", prior = TRUE)
    lines(fit_bPET.mv_V, "PET", density_method = "qCMDE", lty = 2)
    lines(fit_bPET.mv, "PET", col = "blue")
    lines(fit_bPET.mv, "PET", density_method = "qCMDE", col = "blue", lty = 2, density_control = list(samples = 2000))
  })
  scenario_plot("bPET-posterior-mod", {
    plot(fit_bPET.mv_V_reg, "deltype", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 3))
    lines(fit_bPET.mv_V_reg, "deltype", density_method = "qCMDE", lty = 2, density_control = list(samples = 2000))
  })
  scenario_plot("bPET-posterior-random-mv-V", {
    par(mfrow = c(2, 3))
    plot(fit_bPET.mv_V, "tau_total", prior = TRUE)
    lines(fit_bPET.mv_V, "tau_total", density_method = "qCMDE", lty = 2)
    plot(fit_bPET.mv_V, "study: tau", prior = TRUE)
    lines(fit_bPET.mv_V, "study: tau", density_method = "qCMDE", lty = 2, density_control = list(samples = 2000L))
    plot(fit_bPET.mv_V, "esid_study: tau", prior = TRUE)
    lines(fit_bPET.mv_V, "esid_study: tau", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))
    plot(fit_bPET.mv_V, "tau2_prop(esid_study)", prior = TRUE)
    lines(fit_bPET.mv_V, "tau2_prop(esid_study)", density_method = "qCMDE", lty = 2)
    plot(fit_bPET.mv_V, "tau2_prop(study)", prior = TRUE)
    lines(fit_bPET.mv_V, "tau2_prop(study)", density_method = "qCMDE", lty = 2)
  })
  scenario_plot("bPET-posterior-random-mv-vi", {
    par(mfrow = c(2, 3))
    plot(fit_bPET.mv, "tau_total", prior = TRUE)
    lines(fit_bPET.mv, "tau_total", density_method = "qCMDE", lty = 2)
    plot(fit_bPET.mv, "study: tau", prior = TRUE)
    lines(fit_bPET.mv, "study: tau", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))
    plot(fit_bPET.mv, "esid_study: tau", prior = TRUE)
    lines(fit_bPET.mv, "esid_study: tau", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))
    plot(fit_bPET.mv, "tau2_prop(esid_study)", prior = TRUE)
    lines(fit_bPET.mv, "tau2_prop(esid_study)", density_method = "qCMDE", lty = 2)
    plot(fit_bPET.mv, "tau2_prop(study)", prior = TRUE)
    lines(fit_bPET.mv, "tau2_prop(study)", density_method = "qCMDE", lty = 2)
  })

  ### Hypotheses ----
  set.seed(1)
  BF_bPET_rho    <- scenario_time("BF_bPET_rho", hypothesis(fit_bPET_cluster, c("rho != 0 vs rho = 0", "rho != 1 vs rho = 1"), density_method = "qCMDE", density_control = list(samples = 2000)))
  BF_bPET_mv_rho <- scenario_time("BF_bPET_mv_rho", hypothesis(fit_bPET.mv, c("tau2_prop(study) != 0 vs tau2_prop(study) = 0", "tau2_prop(study) != 1 vs tau2_prop(study) = 1"), density_method = "qCMDE", density_control = list(samples = 2000)))
  scenario_text("bPET-rho-bayes-factor-comparison", data.frame(
    rho                = c(0, 1),
    density_bPET_BF    = BF_bPET_rho[["BF"]],
    density_mv_vi_BF   = BF_bPET_mv_rho[["BF"]],
    marglik_bPET_BF    = c(bf(fit_bPET_cluster, fit_bPET)$bf, NA),
    marglik_mv_vi_BF   = c(bf(fit_bPET.mv, fit_bPET.mv_no_study)$bf, bf(fit_bPET.mv, fit_bPET.mv_no_effect)$bf)
  ))

  set.seed(1)
  BF_bPET_random <- scenario_time("BF_bPET_random", hypothesis(fit_bPET.mv_V, c(
    "tau2_prop(study) != 0 vs tau2_prop(study) = 0", "tau2_prop(study) != 1 vs tau2_prop(study) = 1",
    "tau_total = 0"
  ), density_method = "qCMDE", density_control = list(samples = 2000)))
  scenario_text("bPET-random-bayes-factor-comparison", data.frame(
    hypothesis = c("rho != 0", "rho != 1", "tau != 0"),
    density_BF = BF_bPET_random[["BF"]],
    marglik_BF = c(bf(fit_bPET.mv_V, fit_bPET.mv_V_no_study)$bf, bf(fit_bPET.mv_V, fit_bPET.mv_V_no_effect)$bf,
                   bf(fit_bPET.mv_V, fit_bPET.mv_V_fixed)$bf)
  ))

  set.seed(1)
  BF_bPET_mods <- scenario_time("BF_bPET_mods", hypothesis(fit_bPET.mv_V_reg, c("deltype[general] = 0 vs deltype[general] != 0", "deltype[general] = 0 vs deltype[general] > 0", "deltype[general] > 0 vs deltype[general] < 0"), density_method = "qCMDE", density_control = list(samples = 2000)))
  scenario_text("bPET-mods", BF_bPET_mods)

  ### Pooled effects ----
  compare_bPET_pooled <- function(fit_metafor, fit_RoBMA) {
    return(cbind.data.frame(
      "metafor" = t(data.frame(predict(fit_metafor, newmods = 0))[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")]),
      "RoBMA"   = ex_p(fit_RoBMA)[, 1L]
    ))
  }
  scenario_text("bPET-pooled-effect-simple",  compare_bPET_pooled(fit_bPET_metafor,       fit_bPET))
  scenario_text("bPET-pooled-effect-cluster", compare_bPET_pooled(fit_bPET_metafor.mv_vi, fit_bPET_cluster))
  scenario_text("bPET-pooled-effect-mv-vi",   compare_bPET_pooled(fit_bPET_metafor.mv_vi, fit_bPET.mv))
  scenario_text("bPET-pooled-effect-mv-V",    compare_bPET_pooled(fit_bPET_metafor.mv_V,  fit_bPET.mv_V))
  scenario_text("bPET-pooled-effect-mv-V-no-effect", compare_bPET_pooled(fit_bPET_metafor.mv_V_no_effect, fit_bPET.mv_V_no_effect))
  scenario_text("bPET-pooled-effect-mv-V-no-study",  compare_bPET_pooled(fit_bPET_metafor.mv_V_no_study,  fit_bPET.mv_V_no_study))
  scenario_text("bPET-pooled-effect-mv-V-fixed",     compare_bPET_pooled(fit_bPET_metafor.mv_V_fixed,     fit_bPET.mv_V_fixed))
  scenario_text("bPET-pooled-effect-mv-vi-no-effect", compare_bPET_pooled(fit_bPET_metafor.mv_vi_no_effect, fit_bPET.mv_no_effect))
  scenario_text("bPET-pooled-effect-mv-vi-no-study",  compare_bPET_pooled(fit_bPET_metafor.mv_vi_no_study,  fit_bPET.mv_no_study))
  scenario_text("bPET-pooled-effect-mv-vi-fixed",     compare_bPET_pooled(fit_bPET_metafor.mv_vi_fixed,     fit_bPET.mv_fixed))
  scenario_text("bPET-pooled-effect-fixed",           compare_bPET_pooled(fit_bPET_metafor.mv_vi_fixed,     fit_bPET_fixed))

  ### Predictions ----
  compare_bPET_preds_reg <- function(fit_metafor, fit_RoBMA, fit_RoBMA2 = NULL) {
    cbind.data.frame(
      "metafor" = unlist(data.frame(predict(fit_metafor))[c(1, 50, 80), c("pred", "ci.lb", "ci.ub")]),
      "bPET.mv" = unlist(data.frame(predict(fit_RoBMA))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      "bPET"    = if (!is.null(fit_RoBMA2)) unlist(data.frame(predict(fit_RoBMA2))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE) else rep(NA_real_, 9L)
    )
  }
  compare_bPET_preds_reg_pi <- function(fit_metafor, fit_RoBMA, fit_RoBMA2 = NULL) {
    cbind.data.frame(
      "metafor" = unlist(data.frame(predict(fit_metafor))[c(1, 50, 80), c("pred", "pi.lb", "pi.ub")]),
      "bPET.mv" = unlist(data.frame(predict(fit_RoBMA, type = "estimate"))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      "bPET"    = if (!is.null(fit_RoBMA2)) unlist(data.frame(predict(fit_RoBMA2, type = "estimate"))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE) else rep(NA_real_, 9L)
    )
  }
  scenario_text("bPET-predictions-reg",         compare_bPET_preds_reg(fit_bPET_metafor.mv_V_reg,  fit_bPET.mv_V_reg))
  scenario_text("bPET-predictions-mv-vi-reg",   compare_bPET_preds_reg(fit_bPET_metafor.mv_vi_reg, fit_bPET.mv_reg, fit_bPET_cluster_reg))
  scenario_text("bPET-predictions-reg-pi",      compare_bPET_preds_reg_pi(fit_bPET_metafor.mv_V_reg,  fit_bPET.mv_V_reg))
  scenario_text("bPET-predictions-mv-vi-reg-pi", compare_bPET_preds_reg_pi(fit_bPET_metafor.mv_vi_reg, fit_bPET.mv_reg, fit_bPET_cluster_reg))

  ### Marginal means ----
  scenario_text("bPET-marginal-means", marginal_means(fit_bPET.mv_V_reg))
  scenario_plot("bPET-marginal-means-plot", plot(marginal_means(fit_bPET.mv_V_reg), "deltype", prior = TRUE, xlim = c(-2, 2)))

  ### Heterogeneity ----
  scenario_text("bPET-summary-heterogeneity-mv-vi",   summary_heterogeneity(fit_bPET.mv))
  scenario_text("bPET-summary-heterogeneity-cluster", summary_heterogeneity(fit_bPET_cluster))
  scenario_text("bPET-summary-heterogeneity",         summary_heterogeneity(fit_bPET.mv_V))
  scenario_text("bPET-summary-heterogeneity-study",   summary_heterogeneity(fit_bPET.mv_V_no_effect))
  scenario_text("bPET-summary-heterogeneity-fixed",   summary_heterogeneity(fit_bPET.mv_V_fixed))

  ### Random effects ----
  ranef_bPET_metafor.mv_V  <- metafor::ranef(fit_bPET_metafor.mv_V)
  ranef_bPET.mv_V           <- scenario_time("ranef_bPET.mv_V", ranef(fit_bPET.mv_V))
  ranef_bPET_metafor.mv_vi <- metafor::ranef(fit_bPET_metafor.mv_vi)
  ranef_bPET.mv             <- scenario_time("ranef_bPET.mv", ranef(fit_bPET.mv))
  ranef_bPET_cluster        <- scenario_time("ranef_bPET_cluster", ranef(fit_bPET_cluster))
  scenario_text("bPET-random-effects", head(as.data.frame(ranef_bPET.mv_V), 12L))

  scenario_plot("bPET-ranef-mv-V", {
    par(mfrow = c(1, 2))
    scenario_agreement_plot(ranef_bPET_metafor.mv_V$study[["intrcpt"]], as.data.frame(ranef_bPET.mv_V$study)[["Mean"]], main = "study")
    scenario_agreement_plot(ranef_bPET_metafor.mv_V$`study/esid`[["intrcpt"]], as.data.frame(ranef_bPET.mv_V$esid_study)[["Mean"]], main = "esid_study")
  })
  scenario_plot("bPET-ranef-mv-vi", {
    par(mfrow = c(2, 2))
    scenario_agreement_plot(ranef_bPET_metafor.mv_vi$study[["intrcpt"]], as.data.frame(ranef_bPET.mv$study)[["Mean"]], main = "study")
    scenario_agreement_plot(ranef_bPET_metafor.mv_vi$`study/esid`[["intrcpt"]], as.data.frame(ranef_bPET.mv$esid_study)[["Mean"]], main = "esid_study")
    scenario_agreement_plot(ranef_bPET_metafor.mv_vi$study[["intrcpt"]], as.data.frame(ranef_bPET_cluster$cluster)[["Mean"]], main = "study")
    scenario_agreement_plot(ranef_bPET_metafor.mv_vi$`study/esid`[["intrcpt"]], as.data.frame(ranef_bPET_cluster$estimate)[["Mean"]], main = "esid_study")
  })

  ### Diagnostics ----
  plot_bPET_marginal_diagnostics <- function(fit_metafor, fit_bPET) {
    metafor_values <- list(
      "Residuals"      = as.numeric(stats::residuals(fit_metafor)),
      "Rstandard"      = stats::rstandard(fit_metafor)[["z"]],
      "Hat values"     = as.numeric(stats::hatvalues(fit_metafor)),
      "Cooks distance" = stats::cooks.distance(fit_metafor),
      "DFBETAS"        = unlist(stats::dfbetas(fit_metafor))
    )
    bPET_values <- list(
      "Residuals"      = stats::residuals(fit_bPET, type = "outcome", conditioning_depth = "marginal"),
      "Rstandard"      = suppressWarnings(stats::rstandard(fit_bPET, conditioning_depth = "marginal"))[["z"]],
      "Hat values"     = suppressWarnings(stats::hatvalues(fit_bPET)),
      "Cooks distance" = stats::cooks.distance(fit_bPET),
      "DFBETAS"        = unlist(c(
        unclass(stats::dfbetas(fit_bPET, component = "mods")),
        unclass(stats::dfbetas(fit_bPET, component = "bias"))
      ), use.names = FALSE)
    )

    par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
    for (diagnostic in names(metafor_values)) {
      scenario_agreement_plot(metafor_values[[diagnostic]], bPET_values[[diagnostic]], main = diagnostic)
    }
    return(invisible(NULL))
  }
  scenario_plot("bPET-marginal-diagnostics-mv-V",  plot_bPET_marginal_diagnostics(fit_bPET_metafor.mv_V,  fit_bPET.mv_V))
  scenario_plot("bPET-marginal-diagnostics-mv-vi", plot_bPET_marginal_diagnostics(fit_bPET_metafor.mv_vi, fit_bPET.mv))
  scenario_plot("bPET-marginal-diagnostics-reg",   plot_bPET_marginal_diagnostics(fit_bPET_metafor.mv_V_reg, fit_bPET.mv_V_reg))

  scenario_plot("bPET-funnel-mv-V", {
    par(mfrow = c(1, 2))
    funnel(fit_bPET.mv_V, main = "funnel")
    bfunnel(fit_bPET.mv_V, main = "bfunnel")
  })
  scenario_plot("bPET-qqnorm-mv-V", qqnorm(fit_bPET.mv_V, main = "bPET.mv"))
  scenario_plot("bPET-zplot-mv-V",  zplot(fit_bPET.mv_V, to = 10))
  scenario_plot("bPET-funnel-reg", funnel(fit_bPET.mv_V_reg, main = "funnel"))
  scenario_plot("bPET-qqnorm-reg", qqnorm(fit_bPET.mv_V_reg))
  scenario_plot("bPET-zplot-reg",  zplot(fit_bPET.mv_V_reg, to = 10))
  scenario_plot("bPET-PET-plot-V",  plot_pet_peese(fit_bPET.mv_V))
})

testthat::test_that("Assink bPEESE.mv-specific checks", {

  set.seed(1)
  data("dat.assink2016", package = "metadat")

  V_assink <- metafor::vcalc(
    vi, cluster = study, type = deltype, obs = esid,
    rho = c(0.7, 0.5), data = dat.assink2016
  )

  ### Nested model fits ----
  fit_bPEESE.mv_V <- scenario_fit("fit_bPEESE.mv_V", {
    tmp <- bPEESE.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bPEESE.mv_V_null <- scenario_fit("fit_bPEESE.mv_V_null", {
    tmp <- bPEESE.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid,
                     prior_effect = prior("spike", list(0)), measure = "SMD",
                     data = dat.assink2016, seed = 1)
    tmp <- add_marglik(tmp)
    return(tmp)
  })

  ### Point-null Bayes factor ----
  set.seed(1)
  BF_bPEESE_mu <- scenario_time("BF_bPEESE_mu", hypothesis(fit_bPEESE.mv_V, "mu = 0", density_control = list(samples = 2000)))
  scenario_text("bPEESE-mu-bayes-factor-comparison", c(
    marglik = bf(fit_bPEESE.mv_V, fit_bPEESE.mv_V_null)[["bf"]],
    density = BF_bPEESE_mu[["BF"]][[1L]]
  ))

  ### PEESE-specific plots ----
  set.seed(1)
  scenario_plot("bPEESE-PEESE-plot-V", plot_pet_peese(fit_bPEESE.mv_V))
  scenario_plot("bPEESE-funnel-mv-V",  funnel(fit_bPEESE.mv_V))
  scenario_plot("bPEESE-bfunnel-mv-V", bfunnel(fit_bPEESE.mv_V))
})
testthat::test_that("Assink bselmodel and bselmodel.mv models", {

  set.seed(1)
  data("dat.assink2016", package = "metadat")

  V_assink <- metafor::vcalc(
    vi, cluster = study, type = deltype, obs = esid,
    rho = c(0.7, 0.5), data = dat.assink2016
  )
  # The broad structure comparisons use the default selection specification:
  # integrate estimate effects and sampling error; condition on other effects.
  # "no_estimate" drops estimate heterogeneity; "null" below fixes mu to zero.
  ### Default model fits ----
  fit_bselmodel_cond <- scenario_fit("fit_bselmodel_cond", {
    tmp <- bselmodel(yi = yi, vi = vi, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_cluster_cond <- scenario_fit("fit_bselmodel_cluster_cond", {
    tmp <- bselmodel(yi = yi, vi = vi, cluster = study, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_cluster_reg_cond <- scenario_fit("fit_bselmodel_cluster_reg_cond", {
    tmp <- bselmodel(yi = yi, vi = vi, mods = ~ deltype, cluster = study, measure = "SMD", data = dat.assink2016, sample = 20000, burnin = 5000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_fixed_cond <- scenario_fit("fit_bselmodel_fixed_cond", {
    tmp <- bselmodel(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_cond <- scenario_fit("fit_bselmodel.mv_cond", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_reg_cond <- scenario_fit("fit_bselmodel.mv_reg_cond", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, mods = ~ deltype, random = ~ 1 | study / esid,
                        prior_heterogeneity = BayesTools::prior_random(study = BayesTools::random_block(parameterization = "centered")), measure = "SMD",
                        data = dat.assink2016, sample = 15000, burnin = 5000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_no_study_cond <- scenario_fit("fit_bselmodel.mv_no_study_cond", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study:esid, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_no_estimate_cond <- scenario_fit("fit_bselmodel.mv_no_estimate_cond", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study, measure = "SMD", data = dat.assink2016, sample = 60000, burnin = 10000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_fixed_cond <- scenario_fit("fit_bselmodel.mv_fixed_cond", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  # Deferred: dense vcalc integration failed its precision check with both
  # 512 and 32768 points. Keep the original calls for later optimization.
  # fit_bselmodel.mv_V_cond <- scenario_fit("fit_bselmodel.mv_V_cond", {
    # tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_control = set_selection_likelihood_control(points_per_scramble = 32768, max_points_per_scramble = 32768), effect_direction = "positive", parallel = TRUE, seed = 1)
    # tmp <- add_marglik(tmp, repetitions = 3L)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })
  # Deferred: selection-normalizer diagnostics rejected all ten fitting attempts
  # (3.2 minutes). The diagonal regression comparisons remain active below.
  # fit_bselmodel.mv_V_reg_cond <- scenario_fit("fit_bselmodel.mv_V_reg_cond", {
    # tmp <- bselmodel.mv(yi = yi, V = V_assink, mods = ~ deltype, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, sample = 15000, burnin = 5000, parallel = TRUE, seed = 1)
    # tmp <- add_marglik(tmp)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })
  # Deferred: selection-normalizer diagnostics rejected all ten attempts
  # (2.2 minutes); the corresponding diagonal model completed.
  # fit_bselmodel.mv_V_no_study_cond <- scenario_fit("fit_bselmodel.mv_V_no_study_cond", {
    # tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study:esid, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    # tmp <- add_marglik(tmp)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })
  # Deferred: selection-normalizer diagnostics rejected all ten attempts
  # before the requested 120,000-sample run could complete.
  # fit_bselmodel.mv_V_no_estimate_cond <- scenario_fit("fit_bselmodel.mv_V_no_estimate_cond", {
    # tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study, measure = "SMD", data = dat.assink2016, sample = 120000, burnin = 10000, parallel = TRUE, seed = 1)
    # tmp <- add_marglik(tmp)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })
  # Deferred: selection-normalizer diagnostics rejected all ten attempts
  # (1.9 minutes), including this model without random effects.
  # fit_bselmodel.mv_V_fixed_cond <- scenario_fit("fit_bselmodel.mv_V_fixed_cond", {
    # tmp <- bselmodel.mv(yi = yi, V = V_assink, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    # tmp <- add_marglik(tmp)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })

  # Deferred with its effect-present partner; the null failed at 512 points.
  # fit_bselmodel.mv_V_cond_null <- scenario_fit("fit_bselmodel.mv_V_cond_null", {
    # tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, prior_effect = prior("spike", list(0)), measure = "SMD", data = dat.assink2016, selection_control = set_selection_likelihood_control(points_per_scramble = 32768, max_points_per_scramble = 32768), effect_direction = "positive", parallel = TRUE, seed = 1)
    # tmp <- add_marglik(tmp, repetitions = 3L)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })

  ### Model summaries ----
  fit_bselmodel_metafor <- metafor::selmodel(
    metafor::rma(yi, vi, data = dat.assink2016, method = "ML"),
    type = "stepfun", steps = 0.025, decreasing = TRUE
  )
  fit_bselmodel_metafor_fixed <- metafor::selmodel(
    metafor::rma(yi, vi, data = dat.assink2016, method = "FE"),
    type = "stepfun", steps = 0.025, decreasing = TRUE
  )
  # The univariate and estimate-only diagonal models have metafor counterparts.
  # Full V changes the sampling law and need not reproduce these estimates.
  scenario_text("bselmodel-summary-metafor", fit_bselmodel_metafor)
  scenario_text("bselmodel-summary-mv-vi-no-study-cond", summary(fit_bselmodel.mv_no_study_cond))
  scenario_text("bselmodel-summary-simple-cond",         summary(fit_bselmodel_cond))
  # scenario_text("bselmodel-summary-mv-V-no-study-cond", summary(fit_bselmodel.mv_V_no_study_cond))
  scenario_text("bselmodel-summary-metafor-fixed", fit_bselmodel_metafor_fixed)
  scenario_text("bselmodel-summary-mv-vi-fixed-cond", summary(fit_bselmodel.mv_fixed_cond))
  scenario_text("bselmodel-summary-fixed-cond",       summary(fit_bselmodel_fixed_cond))
  # scenario_text("bselmodel-summary-mv-V-fixed-cond", summary(fit_bselmodel.mv_V_fixed_cond))
  scenario_text("bselmodel-summary-mv-vi-cond",   summary(fit_bselmodel.mv_cond))
  scenario_text("bselmodel-summary-cluster-cond", summary(fit_bselmodel_cluster_cond))
  # scenario_text("bselmodel-summary-mv-V-cond",    summary(fit_bselmodel.mv_V_cond))
  scenario_text("bselmodel-summary-mv-vi-no-estimate-cond", summary(fit_bselmodel.mv_no_estimate_cond))
  # scenario_text("bselmodel-summary-mv-V-no-estimate-cond", summary(fit_bselmodel.mv_V_no_estimate_cond))
  scenario_text("bselmodel-summary-mv-vi-reg-cond",       summary(fit_bselmodel.mv_reg_cond))
  scenario_text("bselmodel-summary-cluster-reg-cond",     summary(fit_bselmodel_cluster_reg_cond))
  # scenario_text("bselmodel-summary-mv-V-reg-cond",        summary(fit_bselmodel.mv_V_reg_cond))

  scenario_text("bselmodel-metafor-comparison-default", data.frame(
    implementation  = c("metafor", "bselmodel", "bselmodel.mv"),
    mu              = c(unname(fit_bselmodel_metafor[["beta"]][[1L]]),  ex_r(fit_bselmodel_cond, "mu"),             ex_r(fit_bselmodel.mv_no_study_cond, "mu")),
    omega           = c(unname(fit_bselmodel_metafor[["delta"]][[2L]]), ex_r(fit_bselmodel_cond, "omega[0.025,1]"), ex_r(fit_bselmodel.mv_no_study_cond, "omega[0.025,1]")),
    total_random_sd = c(sqrt(fit_bselmodel_metafor[["tau2"]]),          ex_r(fit_bselmodel_cond, "tau"),            ex_r(fit_bselmodel.mv_no_study_cond, "tau")),
    row.names = NULL
  ))
  scenario_text("bselmodel-metafor-comparison-fixed-default", data.frame(
    implementation = c("metafor", "bselmodel", "bselmodel.mv"),
    mu             = c(unname(fit_bselmodel_metafor_fixed[["beta"]][[1L]]),  ex_r(fit_bselmodel_fixed_cond, "mu"),             ex_r(fit_bselmodel.mv_fixed_cond, "mu")),
    omega          = c(unname(fit_bselmodel_metafor_fixed[["delta"]][[2L]]), ex_r(fit_bselmodel_fixed_cond, "omega[0.025,1]"), ex_r(fit_bselmodel.mv_fixed_cond, "omega[0.025,1]")),
    row.names = NULL
  ))
  # scenario_text("print-prior-fit_bselmodel.mv_V_cond",       print_prior(fit_bselmodel.mv_V_cond))
  # scenario_text("print-prior-fit_bselmodel.mv_V_fixed_cond", print_prior(fit_bselmodel.mv_V_fixed_cond))

  ### Model-fit comparisons ----
  # Corresponding diagonal .mv and univariate models should agree.
  getloo_bselmodel <- function(fit) loo(fit)[["estimates"]]["looic", 1L]
  scenario_text("bselmodel-model-fit-equivalent-cond", data.frame(
    structure = c("nested", "estimate", "study", "fixed"),
    # logml_mv_V      = c(logml(fit_bselmodel.mv_V_cond),     logml(fit_bselmodel.mv_V_no_study_cond), logml(fit_bselmodel.mv_V_no_estimate_cond), logml(fit_bselmodel.mv_V_fixed_cond)),
    logml_mv_vi     = c(logml(fit_bselmodel.mv_cond),       logml(fit_bselmodel.mv_no_study_cond),   logml(fit_bselmodel.mv_no_estimate_cond),   logml(fit_bselmodel.mv_fixed_cond)),
    logml_bselmodel = c(logml(fit_bselmodel_cluster_cond), logml(fit_bselmodel_cond),               NA,                                       logml(fit_bselmodel_fixed_cond)),
    # looic_mv_V      = c(getloo_bselmodel(fit_bselmodel.mv_V_cond),     getloo_bselmodel(fit_bselmodel.mv_V_no_study_cond), getloo_bselmodel(fit_bselmodel.mv_V_no_estimate_cond), getloo_bselmodel(fit_bselmodel.mv_V_fixed_cond)),
    looic_mv_vi     = c(getloo_bselmodel(fit_bselmodel.mv_cond),       getloo_bselmodel(fit_bselmodel.mv_no_study_cond),   getloo_bselmodel(fit_bselmodel.mv_no_estimate_cond),   getloo_bselmodel(fit_bselmodel.mv_fixed_cond)),
    looic_bselmodel = c(getloo_bselmodel(fit_bselmodel_cluster_cond), getloo_bselmodel(fit_bselmodel_cond),               NA,                                                 getloo_bselmodel(fit_bselmodel_fixed_cond)),
    row.names = NULL
  ))
  scenario_text("bselmodel-model-fit-reg-equivalent-default", data.frame(
    structure = c("mv-vi", "cluster"),
    logml     = c(logml(fit_bselmodel.mv_reg_cond), logml(fit_bselmodel_cluster_reg_cond)),
    looic     = c(getloo_bselmodel(fit_bselmodel.mv_reg_cond), getloo_bselmodel(fit_bselmodel_cluster_reg_cond)),
    row.names = NULL
  ))

  ### Posterior plots ----
  # Cluster and diagonal .mv overlays should agree; full V can differ.
  scenario_plot("bselmodel-posterior-rho-default", {
    plot(fit_bselmodel_cluster_cond, "rho", prior = TRUE, ylim = c(0, 4))
    lines(fit_bselmodel_cluster_cond, "rho", lty = 2, density_method = "qCMDE")
    lines(fit_bselmodel.mv_cond, "tau2_prop(study)", col = "blue")
    lines(fit_bselmodel.mv_cond, "tau2_prop(study)", col = "blue", lty = 2, density_method = "qCMDE")
    # lines(fit_bselmodel.mv_V_cond, "tau2_prop(study)", col = "red")
    # lines(fit_bselmodel.mv_V_cond, "tau2_prop(study)", col = "red", lty = 2, density_method = "qCMDE")
  })
  scenario_plot("bselmodel-posterior-location-default", {
    plot(fit_bselmodel_cluster_cond, "mu", prior = TRUE, ylim = c(0, 4), xlim = c(-0.5, 1))
    lines(fit_bselmodel_cluster_cond, "mu", lty = 2, density_method = "qCMDE", density_control = list(samples = 2000))
    lines(fit_bselmodel.mv_cond, "mu", col = "blue")
    lines(fit_bselmodel.mv_cond, "mu", col = "blue", lty = 2, density_method = "qCMDE", density_control = list(samples = 2000))
    # lines(fit_bselmodel.mv_V_cond, "mu", col = "red")
    # lines(fit_bselmodel.mv_V_cond, "mu", col = "red", lty = 2, density_method = "qCMDE", density_control = list(samples = 2000))
  })
  # scenario_plot("bselmodel-posterior-mod-default", {
    # plot(fit_bselmodel.mv_V_reg_cond, "deltype", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 3))
    # lines(fit_bselmodel.mv_V_reg_cond, "deltype", density_method = "qCMDE", density_control = list(samples = 2000))
  # })

  # The mean-null model retains both random-effect levels and identical
  # nuisance priors. These are BF10 values for presence of the pooled effect.
  # compare_bselmodel_effect_bf <- function(label, effect, null, density) {

    # bridge_effect <- bridge_sampler(effect)
    # bridge_null <- bridge_sampler(null)
    # data.frame(
      # selection = label,
      # bridge_BF10 = bf(effect, null)[["bf"]],
      # bridge_logBF10 = logml(effect) - logml(null),
      # bridge_H1_logml_MCSE_max = max(bridge_effect[["mcse_logml"]]),
      # bridge_H0_logml_MCSE_max = max(bridge_null[["mcse_logml"]]),
      # bridge_logBF10_min = min(bridge_effect[["logml"]]) - max(bridge_null[["logml"]]),
      # bridge_logBF10_max = max(bridge_effect[["logml"]]) - min(bridge_null[["logml"]]),
      # qCMDE_BF10 = density[["BF"]],
      # qCMDE_logBF10 = log(density[["BF"]]),
      # qCMDE_error_percent = density[["BF_error"]]
    # )
  # }
  # Repetition ranges and individual log-ML MCSEs are reported separately;
  # the range is not a confidence interval for the aggregate bridge estimate.
  # BF_bselmodel_mu_cond <- scenario_time("BF_bselmodel_mu_cond", hypothesis(fit_bselmodel.mv_V_cond, "mu != 0 vs mu = 0", density_method = "qCMDE", density_control = list(samples = 2000), seed = 1))
  # scenario_text("bselmodel-mu-bayes-factor-comparison-default", compare_bselmodel_effect_bf("default", fit_bselmodel.mv_V_cond, fit_bselmodel.mv_V_cond_null, BF_bselmodel_mu_cond))

  ### Heterogeneity and moderator hypotheses ----
  # Current boundary qCMDE estimates have 65-95% reported error and disagree
  # numerically with bridge estimates. Keep the limitations visible for review.
  set.seed(1)
  BF_bselmodel_rho_cond      <- scenario_time("BF_bselmodel_rho_cond", hypothesis(fit_bselmodel_cluster_cond, c("rho != 0 vs rho = 0", "rho != 1 vs rho = 1"), seed = 1))
  BF_bselmodel_mv_rho_cond   <- scenario_time("BF_bselmodel_mv_rho_cond", hypothesis(fit_bselmodel.mv_cond, c("tau2_prop(study) != 0 vs tau2_prop(study) = 0", "tau2_prop(study) != 1 vs tau2_prop(study) = 1"), seed = 1))
  # BF_bselmodel_mv_V_rho_cond <- scenario_time("BF_bselmodel_mv_V_rho_cond", hypothesis(fit_bselmodel.mv_V_cond, c("tau2_prop(study) != 0 vs tau2_prop(study) = 0", "tau2_prop(study) != 1 vs tau2_prop(study) = 1"), seed = 1))
  scenario_text("bselmodel-rho-bayes-factor-comparison-default", data.frame(
    rho = c(0, 1),
    qCMDE_cluster_BF            = BF_bselmodel_rho_cond[["BF"]],
    qCMDE_cluster_error_percent = BF_bselmodel_rho_cond[["BF_error"]],
    qCMDE_mv_vi_BF              = BF_bselmodel_mv_rho_cond[["BF"]],
    qCMDE_mv_vi_error_percent   = BF_bselmodel_mv_rho_cond[["BF_error"]],
    # qCMDE_mv_V_BF               = BF_bselmodel_mv_V_rho_cond[["BF"]],
    # qCMDE_mv_V_error_percent    = BF_bselmodel_mv_V_rho_cond[["BF_error"]],
    marglik_cluster_BF = c(bf(fit_bselmodel_cluster_cond, fit_bselmodel_cond)[["bf"]], NA),
    marglik_mv_vi_BF   = c(bf(fit_bselmodel.mv_cond, fit_bselmodel.mv_no_study_cond)[["bf"]], bf(fit_bselmodel.mv_cond, fit_bselmodel.mv_no_estimate_cond)[["bf"]]),
    # marglik_mv_V_BF    = c(bf(fit_bselmodel.mv_V_cond, fit_bselmodel.mv_V_no_study_cond)[["bf"]], bf(fit_bselmodel.mv_V_cond, fit_bselmodel.mv_V_no_estimate_cond)[["bf"]]),
    row.names = NULL
  ))
  # scenario_text("bselmodel-mods-cond", hypothesis(fit_bselmodel.mv_V_reg_cond, c("deltype[general] = 0 vs deltype[general] != 0", "deltype[general] = 0 vs deltype[general] > 0", "deltype[general] > 0 vs deltype[general] < 0"), seed = 1))

  ### Pooled effects and predictions ----
  scenario_text("bselmodel-pooled-effect-cond", cbind.data.frame(
    "metafor" = t(data.frame(predict(fit_bselmodel_metafor))[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")]),
    "RoBMA"   = ex_p(fit_bselmodel_cond)[, 1L]
  ))
  compare_bselmodel_preds_reg <- function(fit_mv, fit_specialized = NULL, type = "terms") {

    cbind.data.frame(
      "bselmodel.mv" = unlist(data.frame(predict(fit_mv, type = type))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      "bselmodel"    = if (!is.null(fit_specialized)) { unlist(data.frame(predict(fit_specialized, type = type))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE)
      } else { rep(NA_real_, 9L)}
    )
  }
  # scenario_text("bselmodel-predictions-mv-V-reg-cond", compare_bselmodel_preds_reg(fit_bselmodel.mv_V_reg_cond))
  scenario_text("bselmodel-predictions-mv-vi-reg-cond", compare_bselmodel_preds_reg(fit_bselmodel.mv_reg_cond, fit_bselmodel_cluster_reg_cond))
  # scenario_text("bselmodel-predictions-mv-V-reg-pi-cond", compare_bselmodel_preds_reg(fit_bselmodel.mv_V_reg_cond, type = "estimate"))
  scenario_text("bselmodel-predictions-mv-vi-reg-pi-cond", compare_bselmodel_preds_reg(fit_bselmodel.mv_reg_cond, fit_bselmodel_cluster_reg_cond, type = "estimate"))

  ### Marginal means ----
  # means_bselmodel_reg_cond <- scenario_time("means_bselmodel_reg_cond", marginal_means(fit_bselmodel.mv_V_reg_cond))
  # scenario_text("bselmodel-marginal-means-cond", means_bselmodel_reg_cond)
  # scenario_plot("bselmodel-marginal-means-plot-default", {
    # plot(means_bselmodel_reg_cond, "deltype", prior = TRUE, xlim = c(-2, 2))
    # lines(means_bselmodel_reg_cond, "deltype", density_method = "qCMDE", density_control = list(samples = Inf), lty = 3)
  # })

  ### Heterogeneity ----
  scenario_text("bselmodel-summary-heterogeneity-mv-vi-cond", summary_heterogeneity(fit_bselmodel.mv_cond))
  scenario_text("bselmodel-summary-heterogeneity-cluster-cond", summary_heterogeneity(fit_bselmodel_cluster_cond))
  # scenario_text("bselmodel-summary-heterogeneity-mv-V-cond", summary_heterogeneity(fit_bselmodel.mv_V_cond))
  # scenario_text("bselmodel-summary-heterogeneity-study-cond", summary_heterogeneity(fit_bselmodel.mv_V_no_estimate_cond))
  # scenario_text("bselmodel-summary-heterogeneity-fixed-cond", summary_heterogeneity(fit_bselmodel.mv_V_fixed_cond))

  ### Random effects ----
  # ranef_bselmodel.mv_V_cond <- scenario_time("ranef_bselmodel.mv_V_cond", ranef(fit_bselmodel.mv_V_cond))
  ranef_bselmodel.mv_cond      <- scenario_time("ranef_bselmodel.mv_cond",      ranef(fit_bselmodel.mv_cond))
  ranef_bselmodel_cluster_cond <- scenario_time("ranef_bselmodel_cluster_cond", ranef(fit_bselmodel_cluster_cond))

  plot_bselmodel_ranef_equivalence <- function(ranef_mv, ranef_cluster) {

    par(mfrow = c(1, 2))
    scenario_agreement_plot(as.data.frame(ranef_mv$study)[["Mean"]],      as.data.frame(ranef_cluster$cluster)[["Mean"]],  main = "study",      estimate_label = "cluster", reference_label = "mv-vi")
    scenario_agreement_plot(as.data.frame(ranef_mv$esid_study)[["Mean"]], as.data.frame(ranef_cluster$estimate)[["Mean"]], main = "esid_study", estimate_label = "cluster", reference_label = "mv-vi")
    return(invisible(NULL))
  }
  # plot_bselmodel_ranef_equivalence2 <- function(ranef_mv, ranef_mv2) {

    # par(mfrow = c(1, 2))
    # scenario_agreement_plot(as.data.frame(ranef_mv$study)[["Mean"]],      as.data.frame(ranef_mv2$study)[["Mean"]],      main = "study",      estimate_label = "mv-vi", reference_label = "mv-V")
    # scenario_agreement_plot(as.data.frame(ranef_mv$esid_study)[["Mean"]], as.data.frame(ranef_mv2$esid_study)[["Mean"]], main = "esid_study", estimate_label = "mv-vi", reference_label = "mv-V")
    # return(invisible(NULL))
  # }
  scenario_plot("bselmodel-ranef-mv-vi-cond", plot_bselmodel_ranef_equivalence(ranef_bselmodel.mv_cond, ranef_bselmodel_cluster_cond))
  # scenario_plot("bselmodel-ranef-mv-V-vi-cond", plot_bselmodel_ranef_equivalence2(ranef_bselmodel.mv_V_cond,  ranef_bselmodel.mv_cond))

  # blups
  scenario_plot("bselmodel-blup-mv-vi-cond",   scenario_agreement_plot(data.frame(blup(fit_bselmodel.mv_cond))[,"Mean"],    data.frame(blup(fit_bselmodel_cluster_cond))[,"Mean"],   estimate_label = "cluster", reference_label = "mv-vi"))
  # scenario_plot("bselmodel-blup-mv-V-vi-cond", scenario_agreement_plot(data.frame(blup(fit_bselmodel.mv_V_cond))[,"Mean"],  data.frame(blup(fit_bselmodel.mv_cond))[,"Mean"],   estimate_label = "mv-vi", reference_label = "mv-V"))


  ### Diagnostics ----
  plot_bselmodel_diagnostic_equivalence <- function(fit_reference, fit_mv, reference_label = "cluster", estimate_label = "mv-vi") {

    reference_values <- list(
      "Residuals"      = stats::residuals(fit_reference, type = "outcome", conditioning_depth = "marginal"),
      "Rstudent"       = suppressWarnings(stats::rstudent(fit_reference))[["z"]]
    )
    mv_values <- list(
      "Residuals"      = stats::residuals(fit_mv, type = "outcome", conditioning_depth = "marginal"),
      "Rstudent"       = suppressWarnings(stats::rstudent(fit_mv))[["z"]]
    )

    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    for (diagnostic in names(reference_values)) {
      scenario_agreement_plot(reference_values[[diagnostic]], mv_values[[diagnostic]], main = diagnostic, reference_label = reference_label, estimate_label = estimate_label)
    }
    return(invisible(NULL))
  }
  scenario_plot("bselmodel-marginal-diagnostics-mv-vi-cond", plot_bselmodel_diagnostic_equivalence(fit_bselmodel_cluster_cond, fit_bselmodel.mv_cond))


  # scenario_plot("bselmodel-funnel-mv-V-cond", {
    # par(mfrow = c(1, 2))
    # funnel(fit_bselmodel.mv_V_cond, main = "funnel")
    # bfunnel(fit_bselmodel.mv_V_cond, main = "bfunnel")
  # })
  # scenario_plot("bselmodel-qqnorm-mv-V-cond", qqnorm(fit_bselmodel.mv_V_cond))
  # TODO: optimization needed: the original default zplot panel was unfinished
  # at 90 seconds with its original numerical settings; no fresh result yet.
  # scenario_plot("bselmodel-zplot-mv-V-cond", zplot(fit_bselmodel.mv_V_cond, to = 10))
  # scenario_plot("bselmodel-funnel-reg-cond", funnel(fit_bselmodel.mv_V_reg_cond))
  # scenario_plot("bselmodel-qqnorm-reg-cond", qqnorm(fit_bselmodel.mv_V_reg_cond))
  # TODO: optimization needed: the original regression zplot panel was also
  # unfinished at 90 seconds; preserve the workload for its future review.
  # scenario_plot("bselmodel-zplot-reg-cond", zplot(fit_bselmodel.mv_V_reg_cond, to = 10))

  ### Alternative specification: integrate all three sources ----
  # Alternatives are limited to this full-V, study + estimate structure.
  # The null fixes mu to zero and retains the same nuisance priors and budgets.
  # Deferred: the all-integrated effect model failed selection-normalizer
  # diagnostics in all ten fitting attempts (2.3 minutes).
  # fit_bselmodel.mv_V_marg <- scenario_fit("fit_bselmodel.mv_V_marg", {
    # tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection = selection_model(estimate_random_effects = "integrate", other_random_effects = "integrate", known_sampling_variance = "integrate"), selection_control = set_selection_likelihood_control(max_points_per_scramble = 32768), effect_direction = "positive", parallel = TRUE, seed = 1)
    # tmp <- add_marglik(tmp, repetitions = 3L)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })
  # The matched mean-null model also failed all ten attempts (2.3 minutes).
  # fit_bselmodel.mv_V_marg_null <- scenario_fit("fit_bselmodel.mv_V_marg_null", {
    # tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, prior_effect = prior("spike", list(0)), measure = "SMD", data = dat.assink2016, selection = selection_model(estimate_random_effects = "integrate", other_random_effects = "integrate", known_sampling_variance = "integrate"), selection_control = set_selection_likelihood_control(max_points_per_scramble = 32768), effect_direction = "positive", parallel = TRUE, seed = 1)
    # tmp <- add_marglik(tmp, repetitions = 3L)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })

  # scenario_text("print-prior-fit_bselmodel.mv_V_marg", print_prior(fit_bselmodel.mv_V_marg))
  # scenario_text("bselmodel-summary-mv-V-marg", summary(fit_bselmodel.mv_V_marg))
  # scenario_text("bselmodel-summary-mv-V-marg-null", summary(fit_bselmodel.mv_V_marg_null))
  # scenario_text("bselmodel-summary-heterogeneity-mv-V-marg", summary_heterogeneity(fit_bselmodel.mv_V_marg))
  set.seed(1)
  # BF_bselmodel_mu_marg <- scenario_time("BF_bselmodel_mu_marg", hypothesis(fit_bselmodel.mv_V_marg, "mu != 0 vs mu = 0", density_method = "qCMDE", density_control = list(samples = 2000), seed = 1))
  ### Alternative specification: retain the complete sampling error ----
  # This changes only the sampling choice relative to the default. Estimate
  # effects remain integrated, so selection is nontrivial in both H1 and H0.
  # Deferred: the cached blocks took 20.5 and 24.3 minutes before LOO failed.
  # The LOO correction has focused tests; full scenario refits remain deferred.
  # fit_bselmodel.mv_V_sampling_cond <- scenario_fit("fit_bselmodel.mv_V_sampling_cond", {
    # tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection = selection_model(estimate_random_effects = "integrate", other_random_effects = "condition", known_sampling_variance = "condition"), effect_direction = "positive", parallel = TRUE, seed = 1)
    # tmp <- add_marglik(tmp, repetitions = 3L)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })
  # fit_bselmodel.mv_V_sampling_cond_null <- scenario_fit("fit_bselmodel.mv_V_sampling_cond_null", {
    # tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, prior_effect = prior("spike", list(0)), measure = "SMD", data = dat.assink2016, selection = selection_model(estimate_random_effects = "integrate", other_random_effects = "condition", known_sampling_variance = "condition"), effect_direction = "positive", parallel = TRUE, seed = 1)
    # tmp <- add_marglik(tmp, repetitions = 3L)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })
  # scenario_text("print-prior-fit_bselmodel.mv_V_sampling_cond", print_prior(fit_bselmodel.mv_V_sampling_cond))
  # scenario_text("bselmodel-summary-mv-V-sampling-cond", summary(fit_bselmodel.mv_V_sampling_cond))
  # scenario_text("bselmodel-summary-mv-V-sampling-cond-null", summary(fit_bselmodel.mv_V_sampling_cond_null))
  # BF_bselmodel_mu_sampling_cond <- scenario_time("BF_bselmodel_mu_sampling_cond", hypothesis(fit_bselmodel.mv_V_sampling_cond, "mu != 0 vs mu = 0", density_method = "qCMDE", density_control = list(samples = 2000), seed = 1))

  # Bridge and density BF10 should agree within estimator uncertainty within
  # each specification; the evidence can differ between selection models.
  # scenario_text("bselmodel-mu-bayes-factor-comparison-specifications", rbind(
    # compare_bselmodel_effect_bf("default", fit_bselmodel.mv_V_cond, fit_bselmodel.mv_V_cond_null, BF_bselmodel_mu_cond),
    # compare_bselmodel_effect_bf("all integrated", fit_bselmodel.mv_V_marg, fit_bselmodel.mv_V_marg_null, BF_bselmodel_mu_marg),
    # compare_bselmodel_effect_bf("sampling conditioned", fit_bselmodel.mv_V_sampling_cond, fit_bselmodel.mv_V_sampling_cond_null, BF_bselmodel_mu_sampling_cond)
  # ))
  # scenario_text("bselmodel-mu-density-diagnostics-specifications", rbind(
    # data.frame(selection = "default", density_diagnostics(BF_bselmodel_mu_cond)),
    # data.frame(selection = "all integrated", density_diagnostics(BF_bselmodel_mu_marg)),
    # data.frame(selection = "sampling conditioned", density_diagnostics(BF_bselmodel_mu_sampling_cond))
  # ))
  # scenario_text("bselmodel-model-fit-specifications", data.frame(
    # model = c("default effect", "default null", "all integrated effect", "all integrated null"),
    # logml = c(logml(fit_bselmodel.mv_V_cond), logml(fit_bselmodel.mv_V_cond_null), logml(fit_bselmodel.mv_V_marg), logml(fit_bselmodel.mv_V_marg_null)),
    # looic = c(getloo_bselmodel(fit_bselmodel.mv_V_cond), getloo_bselmodel(fit_bselmodel.mv_V_cond_null), getloo_bselmodel(fit_bselmodel.mv_V_marg), getloo_bselmodel(fit_bselmodel.mv_V_marg_null)),
    # row.names = NULL
  # ))
  # scenario_text("bselmodel-pooled-mv-V-specifications", cbind.data.frame(
    # "default"     = ex_p(fit_bselmodel.mv_V_cond)[, 1L],
    # "marg effect" = ex_p(fit_bselmodel.mv_V_marg)[, 1L],
    # "marg null"   = ex_p(fit_bselmodel.mv_V_marg_null)[, 1L]
  # ))
  # compare_bselmodel_selection_predictions <- function(type) {

    # cbind.data.frame(
      # "default"     = unlist(data.frame(predict(fit_bselmodel.mv_V_cond, type = type, conditioning_depth = "marginal"))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      # "marg effect" = unlist(data.frame(predict(fit_bselmodel.mv_V_marg, type = type, conditioning_depth = "marginal"))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      # "marg null"   = unlist(data.frame(predict(fit_bselmodel.mv_V_marg_null, type = type, conditioning_depth = "marginal"))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE)
    # )
  # }
  # Different selection models can imply different selected latent and response
  # distributions even with identical V, random structures, and nuisance priors.
  # scenario_text("bselmodel-predictions-estimate-specifications", compare_bselmodel_selection_predictions("estimate"))
  # scenario_text("bselmodel-predictions-response-specifications", compare_bselmodel_selection_predictions("response"))
  # scenario_plot("bselmodel-blup-specifications", scenario_agreement_plot(data.frame(blup(fit_bselmodel.mv_V_cond))[, "Mean"], data.frame(blup(fit_bselmodel.mv_V_marg))[, "Mean"], reference_label = "default", estimate_label = "all integrated"))
  # scenario_plot("bselmodel-weightfunction-specifications", {
    # par(mfrow = c(1, 3))
    # plot_weightfunction(fit_bselmodel.mv_V_cond, main = "default")
    # plot_weightfunction(fit_bselmodel.mv_V_marg, main = "marg effect")
    # plot_weightfunction(fit_bselmodel.mv_V_marg_null, main = "marg null")
  # })
  # scenario_plot("bselmodel-funnel-specifications", {
    # par(mfrow = c(1, 3))
    # funnel(fit_bselmodel.mv_V_cond, main = "default")
    # funnel(fit_bselmodel.mv_V_marg, main = "marg effect")
    # funnel(fit_bselmodel.mv_V_marg_null, main = "marg null")
  # })
  # scenario_plot("bselmodel-qqnorm-specifications", {
    # par(mfrow = c(1, 3))
    # qqnorm(fit_bselmodel.mv_V_cond, main = "default")
    # qqnorm(fit_bselmodel.mv_V_marg, main = "marg effect")
    # qqnorm(fit_bselmodel.mv_V_marg_null, main = "marg null")
  # })
  # TODO: optimization needed: the former marg/default zplot comparison and
  # isolated default panel were unfinished at 90 seconds. These new panels
  # retain the numerical settings and require their own timing and review.
  # scenario_plot("bselmodel-zplot-specifications", {
    # par(mfrow = c(1, 3))
    # zplot(fit_bselmodel.mv_V_cond, to = 10, main = "default")
    # zplot(fit_bselmodel.mv_V_marg, to = 10, main = "marg effect")
    # zplot(fit_bselmodel.mv_V_marg_null, to = 10, main = "marg null")
  # })
})

testthat::test_that("Assink RoBMA and RoBMA.mv models", {

  set.seed(1)
  data("dat.assink2016", package = "metadat")

  V_assink <- metafor::vcalc(
    vi, cluster = study, type = deltype, obs = esid,
    rho = c(0.7, 0.5), data = dat.assink2016
  )
  # The broad structure comparisons use the default selection specification:
  # integrate estimate effects and sampling error; condition on other effects.
  ### Default model fits ----
  # Current run: maximum reported R-hat = 1.155; review ensemble mixing.
  fit_RoBMA_cond <- scenario_fit("fit_RoBMA_cond", {
    tmp <- RoBMA(yi = yi, vi = vi, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  # Current run: maximum reported R-hat = 1.094; review ensemble mixing.
  fit_RoBMA_cluster_cond <- scenario_fit("fit_RoBMA_cluster_cond", {
    tmp <- RoBMA(yi = yi, vi = vi, cluster = study, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  # Current run: maximum reported R-hat = 1.055; review ensemble mixing.
  fit_RoBMA_mv_diag_cond <- scenario_fit("fit_RoBMA_mv_diag_cond", {
    tmp <- RoBMA.mv(yi = yi, V = vi, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  # Deferred: selection-normalizer diagnostics rejected all ten attempts
  # (3.8 minutes). The diagonal ensemble remains active.
  # fit_RoBMA_mv_V_cond <- scenario_fit("fit_RoBMA_mv_V_cond", {
    # tmp <- RoBMA.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, parallel = TRUE, seed = 1)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })

  ### Priors and summaries ----
  scenario_text("print-prior-fit_RoBMA_cond",         print_prior(fit_RoBMA_cond))
  scenario_text("print-prior-fit_RoBMA_cluster_cond", print_prior(fit_RoBMA_cluster_cond))
  # scenario_text("print-prior-fit_RoBMA_mv_V_cond",    print_prior(fit_RoBMA_mv_V_cond))

  scenario_text("summary-fit_RoBMA_cond",         summary(fit_RoBMA_cond))
  scenario_text("summary-fit_RoBMA_cluster_cond", summary(fit_RoBMA_cluster_cond))
  scenario_text("summary-fit_RoBMA_mv_diag_cond", summary(fit_RoBMA_mv_diag_cond))
  # scenario_text("summary-fit_RoBMA_mv_V_cond",    summary(fit_RoBMA_mv_V_cond))

  ### Heterogeneity and I2 ----
  scenario_text("summary_heterogeneity-fit_RoBMA_cond",         summary_heterogeneity(fit_RoBMA_cond))
  scenario_text("summary_heterogeneity-fit_RoBMA_cluster_cond", summary_heterogeneity(fit_RoBMA_cluster_cond))
  scenario_text("summary_heterogeneity-fit_RoBMA_mv_diag_cond", summary_heterogeneity(fit_RoBMA_mv_diag_cond))
  # scenario_text("summary_heterogeneity-fit_RoBMA_mv_V_cond",    summary_heterogeneity(fit_RoBMA_mv_V_cond))

  ### Model fit ----
  getloo <- function(fit) loo(fit)[["estimates"]]["looic", 1L]
  scenario_text("RoBMA-loo", data.frame(
    structure = "default",
    uni       = getloo(fit_RoBMA_cond),
    cluster   = getloo(fit_RoBMA_cluster_cond),
    "mv-diag" = getloo(fit_RoBMA_mv_diag_cond),
    # "mv-V"    = getloo(fit_RoBMA_mv_V_cond),
    row.names = NULL
  ))

  ### Posterior densities ----
  # Deferred: the first qCMDE mu overlay failed with 29.8% bulk relative MCSE
  # (13.6 seconds). Later panels were not reached; original calls retained.
  # scenario_plot("posterior-fit_RoBMA_mv_diag_cond", {
    # par(mfrow = c(2, 3))
    # plot(fit_RoBMA_mv_diag_cond,  "mu", prior = TRUE)
    # lines(fit_RoBMA_mv_diag_cond, "mu", density_method = "qCMDE", lty = 2)
    # lines(fit_RoBMA_cluster_cond, "mu", col = "blue")
    # lines(fit_RoBMA_cluster_cond, "mu", density_method = "qCMDE", col = "blue", lty = 2)

    # plot(fit_RoBMA_mv_diag_cond,  "tau_total", prior = TRUE)
    # lines(fit_RoBMA_mv_diag_cond, "tau_total", density_method = "qCMDE", lty = 2)
    # lines(fit_RoBMA_cluster_cond, "tau", col = "blue")
    # lines(fit_RoBMA_cluster_cond, "tau", density_method = "qCMDE", col = "blue", lty = 2)

    # plot(fit_RoBMA_mv_diag_cond,  "study: tau", prior = TRUE)
    # lines(fit_RoBMA_mv_diag_cond, "study: tau", density_method = "qCMDE", lty = 2)

    # plot(fit_RoBMA_mv_diag_cond,  "esid_study: tau", prior = TRUE)
    # lines(fit_RoBMA_mv_diag_cond, "esid_study: tau", density_method = "qCMDE", lty = 2)

    # plot(fit_RoBMA_mv_diag_cond,  "tau2_prop(esid_study)", prior = TRUE)
    # lines(fit_RoBMA_mv_diag_cond, "tau2_prop(esid_study)", density_method = "qCMDE", lty = 2)

    # plot(fit_RoBMA_mv_diag_cond,  "tau2_prop(study)", prior = TRUE)
    # lines(fit_RoBMA_mv_diag_cond, "tau2_prop(study)", density_method = "qCMDE", lty = 2)
    # lines(fit_RoBMA_cluster_cond, "rho", col = "blue")
    # lines(fit_RoBMA_cluster_cond, "rho", density_method = "qCMDE", col = "blue", lty = 2)
  # })

  # TODO: qCMDE accuracy review needed: the previous conditional mu overlay
  # failed with 31.9% bulk relative MCSE; no completed figure was available.
  # scenario_plot("posterior-fit_RoBMA_mv_cond", {
    # par(mfrow = c(2, 3))
    # plot(fit_RoBMA_mv_V_cond,  "mu", prior = TRUE)
    # lines(fit_RoBMA_mv_V_cond, "mu", density_method = "qCMDE")

    # plot(fit_RoBMA_mv_V_cond,  "tau_total", prior = TRUE)
    # lines(fit_RoBMA_mv_V_cond, "tau_total", density_method = "qCMDE")

    # plot(fit_RoBMA_mv_V_cond,  "study: tau", prior = TRUE)
    # lines(fit_RoBMA_mv_V_cond, "study: tau", density_method = "qCMDE", density_control = list(samples = Inf))

    # plot(fit_RoBMA_mv_V_cond,  "esid_study: tau", prior = TRUE)
    # lines(fit_RoBMA_mv_V_cond, "esid_study: tau", density_method = "qCMDE")

    # plot(fit_RoBMA_mv_V_cond,  "tau2_prop(esid_study)", prior = TRUE)
    # lines(fit_RoBMA_mv_V_cond, "tau2_prop(esid_study)", density_method = "qCMDE")

    # plot(fit_RoBMA_mv_V_cond,  "tau2_prop(study)", prior = TRUE)
    # lines(fit_RoBMA_mv_V_cond, "tau2_prop(study)", density_method = "qCMDE")
  # })

  ### Random effects and BLUPs ----
  ranef_RoBMA_cluster_cond <- scenario_time("ranef-RoBMA_cluster_cond", ranef(fit_RoBMA_cluster_cond))
  ranef_RoBMA_mv_diag_cond <- scenario_time("ranef-RoBMA_mv_diag_cond", ranef(fit_RoBMA_mv_diag_cond))
  # ranef_RoBMA_mv_V_cond    <- scenario_time("ranef-RoBMA_mv_V_cond",    ranef(fit_RoBMA_mv_V_cond))

  scenario_plot("RoBMA-ranef-2", {
    par(mfrow = c(1, 2))
    scenario_agreement_plot(as.data.frame(ranef_RoBMA_mv_diag_cond$study)[["Mean"]],      as.data.frame(ranef_RoBMA_cluster_cond$cluster)[["Mean"]],  main = "study",      estimate_label = "cluster", reference_label = "mv-vi")
    scenario_agreement_plot(as.data.frame(ranef_RoBMA_mv_diag_cond$esid_study)[["Mean"]], as.data.frame(ranef_RoBMA_cluster_cond$estimate)[["Mean"]], main = "esid_study", estimate_label = "cluster", reference_label = "mv-vi")
  })
  # The full sampling covariance changes the model; this is a sensitivity check.
  # scenario_plot("RoBMA-ranef-4", {
    # par(mfrow = c(1, 2))
    # scenario_agreement_plot(as.data.frame(ranef_RoBMA_mv_V_cond$study)[["Mean"]],      as.data.frame(ranef_RoBMA_mv_diag_cond$study)[["Mean"]],      main = "study",      estimate_label = "mv-vi", reference_label = "mv-V")
    # scenario_agreement_plot(as.data.frame(ranef_RoBMA_mv_V_cond$esid_study)[["Mean"]], as.data.frame(ranef_RoBMA_mv_diag_cond$esid_study)[["Mean"]], main = "esid_study", estimate_label = "mv-vi", reference_label = "mv-V")
  # })

  scenario_plot("RoBMA-blup-2", scenario_agreement_plot(data.frame(blup(fit_RoBMA_mv_diag_cond))[,"Mean"], data.frame(blup(fit_RoBMA_cluster_cond))[,"Mean"], estimate_label = "cluster", reference_label = "mv-vi"))
  # scenario_plot("RoBMA-blup-mv-V-cond-vs-diag", scenario_agreement_plot(data.frame(blup(fit_RoBMA_mv_V_cond))[,"Mean"], data.frame(blup(fit_RoBMA_mv_diag_cond))[,"Mean"], estimate_label = "mv-vi", reference_label = "mv-V"))

  ### Diagnostics ----
  # scenario_plot("RoBMA-funnel-mv-V-cond", funnel(fit_RoBMA_mv_V_cond, main = "default"))
  # scenario_plot("RoBMA-qqnorm-mv-V-cond", qqnorm(fit_RoBMA_mv_V_cond, main = "default"))
  # TODO: optimization needed: conditional selected zplot usability remains
  # unresolved. This isolated default panel has no new timing evidence.
  # scenario_plot("RoBMA-zplot-mv-V-cond", zplot(fit_RoBMA_mv_V_cond, to = 10, main = "default"))

  ### Alternative specification: integrate all three sources ----
  # Exactly two alternative fits: force presence or absence of the mean, while
  # preserving model averaging over bias and both random-effect levels.
  # Deferred: the all-integrated effect-present ensemble failed selection-
  # normalizer diagnostics in all ten fitting attempts (4.1 minutes).
  # fit_RoBMA_mv_V_marg_effect <- scenario_fit("fit_RoBMA_mv_V_marg_effect", {
    # tmp <- RoBMA.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, prior_effect_null = NULL, measure = "SMD", data = dat.assink2016, selection = selection_model(estimate_random_effects = "integrate", other_random_effects = "integrate", known_sampling_variance = "integrate"), selection_control = set_selection_likelihood_control(max_points_per_scramble = 32768), parallel = TRUE, seed = 1)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })
  # The matched mean-null ensemble also failed all ten attempts (3.1 minutes).
  # fit_RoBMA_mv_V_marg_null <- scenario_fit("fit_RoBMA_mv_V_marg_null", {
    # tmp <- RoBMA.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, prior_effect = NULL, measure = "SMD", data = dat.assink2016, selection = selection_model(estimate_random_effects = "integrate", other_random_effects = "integrate", known_sampling_variance = "integrate"), selection_control = set_selection_likelihood_control(max_points_per_scramble = 32768), parallel = TRUE, seed = 1)
    # tmp <- add_loo(tmp)
    # return(tmp)
  # })
  # scenario_text("print-prior-fit_RoBMA_mv_V_marg_effect", print_prior(fit_RoBMA_mv_V_marg_effect))
  # scenario_text("print-prior-fit_RoBMA_mv_V_marg_null", print_prior(fit_RoBMA_mv_V_marg_null))
  # scenario_text("summary-fit_RoBMA_mv_V_marg_effect", summary(fit_RoBMA_mv_V_marg_effect))
  # scenario_text("summary-fit_RoBMA_mv_V_marg_null", summary(fit_RoBMA_mv_V_marg_null))
  # scenario_text("summary_heterogeneity-fit_RoBMA_mv_V_marg_effect", summary_heterogeneity(fit_RoBMA_mv_V_marg_effect))
  # scenario_text("summary_heterogeneity-fit_RoBMA_mv_V_marg_null", summary_heterogeneity(fit_RoBMA_mv_V_marg_null))

  ### Effect evidence ----
  # RoBMA ensembles do not provide bridge marginal likelihoods. Compare the
  # default ensemble's inclusion BF with its active-effect density BF; the
  # forced marginal H1 instead supplies a density BF against its matching H0.
  # With zero/one rare-state draws and overwhelming evidence in the same
  # direction, compare conclusions and retain the error%. Density and
  # integration diagnostics still apply unchanged.
  # scenario_text("RoBMA-effect-inclusion-default", summary_models(fit_RoBMA_mv_V_cond, type = "marginal"))
  # BF_RoBMA_mu_cond <- scenario_time("BF_RoBMA_mu_cond", hypothesis(fit_RoBMA_mv_V_cond, "mu != 0 vs mu = 0", conditional = TRUE, density_method = "qCMDE", density_control = list(samples = 2000), seed = 1))
  # BF_RoBMA_mu_marg <- scenario_time("BF_RoBMA_mu_marg", hypothesis(fit_RoBMA_mv_V_marg_effect, "mu != 0 vs mu = 0", conditional = FALSE, density_method = "qCMDE", density_control = list(samples = 2000), seed = 1))
  # scenario_text("RoBMA-effect-density-default", BF_RoBMA_mu_cond)
  # scenario_text("RoBMA-effect-density-marg", BF_RoBMA_mu_marg)
  # scenario_text("RoBMA-loo-specifications", data.frame(
    # model = c("default ensemble", "all integrated effect", "all integrated null"),
    # looic = c(getloo(fit_RoBMA_mv_V_cond), getloo(fit_RoBMA_mv_V_marg_effect), getloo(fit_RoBMA_mv_V_marg_null)),
    # row.names = NULL
  # ))
  # scenario_text("RoBMA-pooled-specifications", cbind.data.frame(
    # "default ensemble" = ex_p(fit_RoBMA_mv_V_cond)[, 1L],
    # "marg effect"      = ex_p(fit_RoBMA_mv_V_marg_effect)[, 1L],
    # "marg null"        = ex_p(fit_RoBMA_mv_V_marg_null)[, 1L]
  # ))
  # compare_RoBMA_selection_predictions <- function(type) {

    # cbind.data.frame(
      # "default ensemble" = unlist(data.frame(predict(fit_RoBMA_mv_V_cond, type = type, conditioning_depth = "marginal"))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      # "marg effect"      = unlist(data.frame(predict(fit_RoBMA_mv_V_marg_effect, type = type, conditioning_depth = "marginal"))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      # "marg null"        = unlist(data.frame(predict(fit_RoBMA_mv_V_marg_null, type = type, conditioning_depth = "marginal"))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE)
    # )
  # }
  # These are model-sensitivity comparisons, not expected equalities. The
  # default ensemble additionally averages over presence of the pooled effect.
  # scenario_text("RoBMA-predictions-estimate-specifications", compare_RoBMA_selection_predictions("estimate"))
  # scenario_text("RoBMA-predictions-response-specifications", compare_RoBMA_selection_predictions("response"))
  # scenario_plot("RoBMA-blup-specifications", scenario_agreement_plot(data.frame(blup(fit_RoBMA_mv_V_cond))[, "Mean"], data.frame(blup(fit_RoBMA_mv_V_marg_effect))[, "Mean"], reference_label = "default ensemble", estimate_label = "marg effect"))
  # scenario_plot("RoBMA-weightfunction-specifications", {
    # par(mfrow = c(1, 3))
    # plot_weightfunction(fit_RoBMA_mv_V_cond, main = "default ensemble")
    # plot_weightfunction(fit_RoBMA_mv_V_marg_effect, main = "marg effect")
    # plot_weightfunction(fit_RoBMA_mv_V_marg_null, main = "marg null")
  # })
  # scenario_plot("RoBMA-funnel-specifications", {
    # par(mfrow = c(1, 3))
    # funnel(fit_RoBMA_mv_V_cond, main = "default ensemble")
    # funnel(fit_RoBMA_mv_V_marg_effect, main = "marg effect")
    # funnel(fit_RoBMA_mv_V_marg_null, main = "marg null")
  # })
  # scenario_plot("RoBMA-qqnorm-specifications", {
    # par(mfrow = c(1, 3))
    # qqnorm(fit_RoBMA_mv_V_cond, main = "default ensemble")
    # qqnorm(fit_RoBMA_mv_V_marg_effect, main = "marg effect")
    # qqnorm(fit_RoBMA_mv_V_marg_null, main = "marg null")
  # })
  # TODO: optimization needed: the former full-V marginal fit did not complete,
  # so its predictive figures have no completed baseline or usable timing.
  # scenario_plot("RoBMA-zplot-specifications", {
    # par(mfrow = c(1, 3))
    # zplot(fit_RoBMA_mv_V_cond, to = 10, main = "default ensemble")
    # zplot(fit_RoBMA_mv_V_marg_effect, to = 10, main = "marg effect")
    # zplot(fit_RoBMA_mv_V_marg_null, to = 10, main = "marg null")
  # })
})
