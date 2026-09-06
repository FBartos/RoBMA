if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("assink2016")
# testthat::test_file("tests/scenarios/test-assink2016.R")
# in progress -- do not change anything here
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

  metafor_parameters <- c("intercept", study_variance = "sigma[study]^2", estimate_variance = "sigma[study/esid]^2", total_random_variance = "sigma[total]^2")
  robma_parameters   <- c("mu", study_variance = "var", estimate_variance = "var", total_random_variance = "var_total")
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

    lines(fit_brma.mv_diag, "var_prop(study)", col = "blue")
    lines(fit_brma.mv_diag, "var_prop(study)", density_method = "qCMDE", col = "blue", lty = 2)
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

    plot(fit_brma.mv, "sd_total", prior = TRUE)
    lines(fit_brma.mv, "sd_total", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "study: sd", prior = TRUE)
    lines(fit_brma.mv, "study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_brma.mv, "esid_study: sd", prior = TRUE)
    lines(fit_brma.mv, "esid_study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_brma.mv, "var_prop(esid_study)", prior = TRUE)
    lines(fit_brma.mv, "var_prop(esid_study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "var_prop(study)", prior = TRUE)
    lines(fit_brma.mv, "var_prop(study)", density_method = "qCMDE", lty = 2)
  })

  set.seed(1)
  scenario_plot("fit_BMA.mv_posterior", {
    par(mfrow = c(2, 3))

    plot(fit_BMA.mv, "mu")
    lines(fit_BMA.mv, "mu", density_method = "qCMDE", lty = 2)

    plot(fit_BMA.mv, "sd_total", prior = TRUE)
    lines(fit_BMA.mv, "sd_total", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_BMA.mv, "esid_study: sd", prior = TRUE)
    lines(fit_BMA.mv, "esid_study: sd", density_method = "qCMDE", lty = 2)

    plot(fit_BMA.mv, "var_prop(esid_study)", prior = TRUE)
    lines(fit_BMA.mv, "var_prop(esid_study)", density_method = "qCMDE", lty = 2)

    plot(fit_BMA.mv, "var_prop(study)", prior = TRUE)
    lines(fit_BMA.mv, "var_prop(study)", density_method = "qCMDE", lty = 2)
  })

  set.seed(1)
  scenario_plot("fit.BMA_posterior_random", {
    par(mfrow = c(2, 3))

    plot(fit_brma.mv, "sd_total", prior = TRUE)
    lines(fit_brma.mv, "sd_total", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "study: sd", prior = TRUE)
    lines(fit_brma.mv, "study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_brma.mv, "esid_study: sd", prior = TRUE)
    lines(fit_brma.mv, "esid_study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_brma.mv, "var_prop(esid_study)", prior = TRUE)
    lines(fit_brma.mv, "var_prop(esid_study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv, "var_prop(study)", prior = TRUE)
    lines(fit_brma.mv, "var_prop(study)", density_method = "qCMDE", lty = 2)
  })



  ### hypothesis ----
  set.seed(1)
  BF_brma_rho    <- scenario_time("BF_brma_rho", hypothesis(fit_brma_cluster, c("rho != 0 vs rho = 0", "rho != 1 vs rho = 1"), density_method = "qCMDE", density_control = list(samples = 2000)))
  BF_mv_diag_rho <- scenario_time("BF_mv_diag_rho", hypothesis(fit_brma.mv_diag, c("var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1"), density_method = "qCMDE", density_control = list(samples = 2000)))
  scenario_text("fit_rho_bayes_factor_comparison", data.frame(
    rho                = c(0, 1),
    density_brma_BF    = BF_brma_rho[["BF"]],
    density_mv_diag_BF = BF_mv_diag_rho[["BF"]],
    marglik_brma_BF    = c(bf(fit_brma_cluster, fit_brma)$bf, NA),
    marglik_mv_diag_BF = c(bf(fit_brma.mv_diag, fit_brma.mv_diag_no_study)$bf, bf(fit_brma.mv_diag, fit_brma.mv_diag_no_effect)$bf)
  ))

  set.seed(1)
  BF_random      <- scenario_time("BF_random", hypothesis(fit_brma.mv, c(
    "var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1",
    "sd_total = 0"
    ),density_method = "qCMDE", density_control = list(samples = 2000)))
  scenario_text("fit_random_bayes_factor_comparison", data.frame(
    hypothesis = c("rho != 0", "rho != 1", "sd != 0"),
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
    tmp <- bPET.mv(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, seed = 1)
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
    tmp <- bPET.mv(yi = yi, V = V_assink, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, seed = 1)
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

  metafor_PET_parameters    <- c(mu = "intercept", PET = "sqrt(vi)", total_random_sd = "tau")
  metafor_PET_mv_parameters <- c(mu = "intercept", PET = "sqrt(vi)", total_random_sd = "sigma[total]")
  robma_PET_parameters      <- c(mu = "mu", PET = "PET", total_random_sd = "tau")
  robma_PET_mv_parameters   <- c(mu = "mu", PET = "PET", total_random_sd = "sd_total")
  metafor_PET_mv_vi_study_fraction <- (ex_m(fit_bPET_metafor.mv_vi, "sigma[study]") / ex_m(fit_bPET_metafor.mv_vi, "sigma[total]"))^2
  metafor_PET_mv_V_study_fraction  <- (ex_m(fit_bPET_metafor.mv_V, "sigma[study]") / ex_m(fit_bPET_metafor.mv_V, "sigma[total]"))^2

  scenario_text("bPET-metafor-comparison", data.frame(
    model          = rep(c("simple", "cluster", "mv_vi", "mv_correlated_V"), each = 2L),
    implementation = rep(c("metafor", "RoBMA"), 4L),
    rbind(
      c(ex_m(fit_bPET_metafor, metafor_PET_parameters), study_variance_fraction = NA_real_),
      c(ex_r(fit_bPET, robma_PET_parameters), study_variance_fraction = NA_real_),
      c(ex_m(fit_bPET_metafor.mv_vi, metafor_PET_mv_parameters), study_variance_fraction = metafor_PET_mv_vi_study_fraction),
      c(ex_r(fit_bPET_cluster, robma_PET_parameters), study_variance_fraction = ex_r(fit_bPET_cluster, "rho")),
      c(ex_m(fit_bPET_metafor.mv_vi, metafor_PET_mv_parameters), study_variance_fraction = metafor_PET_mv_vi_study_fraction),
      c(ex_r(fit_bPET.mv, robma_PET_mv_parameters), study_variance_fraction = ex_r(fit_bPET.mv, "var_prop(study)")),
      c(ex_m(fit_bPET_metafor.mv_V, metafor_PET_mv_parameters), study_variance_fraction = metafor_PET_mv_V_study_fraction),
      c(ex_r(fit_bPET.mv_V, robma_PET_mv_parameters), study_variance_fraction = ex_r(fit_bPET.mv_V, "var_prop(study)"))
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
    lines(fit_bPET.mv, "var_prop(study)", col = "blue")
    lines(fit_bPET.mv, "var_prop(study)", density_method = "qCMDE", col = "blue", lty = 2)
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
    plot(fit_bPET.mv_V, "sd_total", prior = TRUE)
    lines(fit_bPET.mv_V, "sd_total", density_method = "qCMDE", lty = 2)
    plot(fit_bPET.mv_V, "study: sd", prior = TRUE)
    lines(fit_bPET.mv_V, "study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 2000L))
    plot(fit_bPET.mv_V, "esid_study: sd", prior = TRUE)
    lines(fit_bPET.mv_V, "esid_study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))
    plot(fit_bPET.mv_V, "var_prop(esid_study)", prior = TRUE)
    lines(fit_bPET.mv_V, "var_prop(esid_study)", density_method = "qCMDE", lty = 2)
    plot(fit_bPET.mv_V, "var_prop(study)", prior = TRUE)
    lines(fit_bPET.mv_V, "var_prop(study)", density_method = "qCMDE", lty = 2)
  })
  scenario_plot("bPET-posterior-random-mv-vi", {
    par(mfrow = c(2, 3))
    plot(fit_bPET.mv, "sd_total", prior = TRUE)
    lines(fit_bPET.mv, "sd_total", density_method = "qCMDE", lty = 2)
    plot(fit_bPET.mv, "study: sd", prior = TRUE)
    lines(fit_bPET.mv, "study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))
    plot(fit_bPET.mv, "esid_study: sd", prior = TRUE)
    lines(fit_bPET.mv, "esid_study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))
    plot(fit_bPET.mv, "var_prop(esid_study)", prior = TRUE)
    lines(fit_bPET.mv, "var_prop(esid_study)", density_method = "qCMDE", lty = 2)
    plot(fit_bPET.mv, "var_prop(study)", prior = TRUE)
    lines(fit_bPET.mv, "var_prop(study)", density_method = "qCMDE", lty = 2)
  })

  ### Hypotheses ----
  set.seed(1)
  BF_bPET_rho    <- scenario_time("BF_bPET_rho", hypothesis(fit_bPET_cluster, c("rho != 0 vs rho = 0", "rho != 1 vs rho = 1"), density_method = "qCMDE", density_control = list(samples = 2000)))
  BF_bPET_mv_rho <- scenario_time("BF_bPET_mv_rho", hypothesis(fit_bPET.mv, c("var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1"), density_method = "qCMDE", density_control = list(samples = 2000)))
  scenario_text("bPET-rho-bayes-factor-comparison", data.frame(
    rho                = c(0, 1),
    density_bPET_BF    = BF_bPET_rho[["BF"]],
    density_mv_vi_BF   = BF_bPET_mv_rho[["BF"]],
    marglik_bPET_BF    = c(bf(fit_bPET_cluster, fit_bPET)$bf, NA),
    marglik_mv_vi_BF   = c(bf(fit_bPET.mv, fit_bPET.mv_no_study)$bf, bf(fit_bPET.mv, fit_bPET.mv_no_effect)$bf)
  ))

  set.seed(1)
  BF_bPET_random <- scenario_time("BF_bPET_random", hypothesis(fit_bPET.mv_V, c(
    "var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1",
    "sd_total = 0"
  ), density_method = "qCMDE", density_control = list(samples = 2000)))
  scenario_text("bPET-random-bayes-factor-comparison", data.frame(
    hypothesis = c("rho != 0", "rho != 1", "sd != 0"),
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

  V_assink <- vcalc2(
    vi, cluster = study, type = deltype, obs = esid,
    rho = c(0.7, 0.5), data = dat.assink2016
  )


  fit_uni_metafor <- metafor::rma(yi, vi, data = dat.assink2016, method = "REML")
  fit_mv_metafor  <- metafor::rma.mv(yi, V_assink, data = dat.assink2016, method = "REML")
  fit_bselmodel_metafor <- metafor::selmodel(
    metafor::rma(yi, vi, data = dat.assink2016, method = "ML"),
    type = "stepfun", steps = 0.025, decreasing = TRUE
  )
  fit_bselmodel_metafor_fixed <- metafor::selmodel(
    metafor::rma(yi, vi, data = dat.assink2016, method = "FE"),
    type = "stepfun", steps = 0.025, decreasing = TRUE
  )

  ### Exact model fits ----
  fit_bselmodel_exact <- scenario_fit("fit_bselmodel_exact", {
    tmp <- bselmodel(yi = yi, vi = vi, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_cluster_exact <- scenario_fit("fit_bselmodel_cluster_exact", {
    tmp <- bselmodel(yi = yi, vi = vi, cluster = study, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_cluster_reg_exact <- scenario_fit("fit_bselmodel_cluster_reg_exact", {
    tmp <- bselmodel(yi = yi, vi = vi, mods = ~ deltype, cluster = study, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_fixed_exact <- scenario_fit("fit_bselmodel_fixed_exact", {
    tmp <- bselmodel(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_exact <- scenario_fit("fit_bselmodel.mv_exact", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_reg_exact <- scenario_fit("fit_bselmodel.mv_reg_exact", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, mods = ~ deltype, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_no_study_exact <- scenario_fit("fit_bselmodel.mv_no_study_exact", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study:esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_no_effect_exact <- scenario_fit("fit_bselmodel.mv_no_effect_exact", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_fixed_exact <- scenario_fit("fit_bselmodel.mv_fixed_exact", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_exact <- scenario_fit("fit_bselmodel.mv_V_exact", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_reg_exact <- scenario_fit("fit_bselmodel.mv_V_reg_exact", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, mods = ~ deltype, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_no_study_exact <- scenario_fit("fit_bselmodel.mv_V_no_study_exact", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study:esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_no_effect_exact <- scenario_fit("fit_bselmodel.mv_V_no_effect_exact", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_fixed_exact <- scenario_fit("fit_bselmodel.mv_V_fixed_exact", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })

  ### Approximate model fits ----
  fit_bselmodel_approximate <- scenario_fit("fit_bselmodel_approximate", {
    tmp <- bselmodel(yi = yi, vi = vi, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_cluster_approximate <- scenario_fit("fit_bselmodel_cluster_approximate", {
    tmp <- bselmodel(yi = yi, vi = vi, cluster = study, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_cluster_reg_approximate <- scenario_fit("fit_bselmodel_cluster_reg_approximate", {
    tmp <- bselmodel(yi = yi, vi = vi, mods = ~ deltype, cluster = study, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_fixed_approximate <- scenario_fit("fit_bselmodel_fixed_approximate", {
    tmp <- bselmodel(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_approximate <- scenario_fit("fit_bselmodel.mv_approximate", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_reg_approximate <- scenario_fit("fit_bselmodel.mv_reg_approximate", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, mods = ~ deltype, random = ~ 1 | study / esid,
                        prior_heterogeneity = BayesTools::prior_random(study = BayesTools::random_block(parameterization = "centered")), measure = "SMD",
                        data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_no_study_approximate <- scenario_fit("fit_bselmodel.mv_no_study_approximate", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study:esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_no_effect_approximate <- scenario_fit("fit_bselmodel.mv_no_effect_approximate", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_fixed_approximate <- scenario_fit("fit_bselmodel.mv_fixed_approximate", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_approximate <- scenario_fit("fit_bselmodel.mv_V_approximate", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_reg_approximate <- scenario_fit("fit_bselmodel.mv_V_reg_approximate", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, mods = ~ deltype, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_no_study_approximate <- scenario_fit("fit_bselmodel.mv_V_no_study_approximate", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study:esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_no_effect_approximate <- scenario_fit("fit_bselmodel.mv_V_no_effect_approximate", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_fixed_approximate <- scenario_fit("fit_bselmodel.mv_V_fixed_approximate", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })

  ### Model summaries ----
  # the V adds much more weight to smaller estimates -- as such, we would expect seeing much lower pooled effect fo V vs vi models
  fit_mv_metafor
  fit_uni_metafor

  # these are the only one with corresponding metafor models
  scenario_text("bselmodel-summary-metafor", fit_bselmodel_metafor)
  scenario_text("bselmodel-summary-mv-vi-no-study-exact",       summary(fit_bselmodel.mv_no_study_exact))
  scenario_text("bselmodel-summary-mv-vi-no-study-approximate", summary(fit_bselmodel.mv_no_study_approximate))
  scenario_text("bselmodel-summary-simple-exact",               summary(fit_bselmodel_exact))
  scenario_text("bselmodel-summary-simple-approximate",         summary(fit_bselmodel_approximate))
  # with V
  scenario_text("bselmodel-summary-mv-V-no-study-exact",       summary(fit_bselmodel.mv_V_no_study_exact))
  scenario_text("bselmodel-summary-mv-V-no-study-approximate", summary(fit_bselmodel.mv_V_no_study_approximate))

  scenario_text("bselmodel-summary-metafor-fixed", fit_bselmodel_metafor_fixed)
  scenario_text("bselmodel-summary-mv-vi-fixed-exact",       summary(fit_bselmodel.mv_fixed_exact))
  scenario_text("bselmodel-summary-mv-vi-fixed-approximate", summary(fit_bselmodel.mv_fixed_approximate))
  scenario_text("bselmodel-summary-fixed-exact",             summary(fit_bselmodel_fixed_exact))
  scenario_text("bselmodel-summary-fixed-approximate",       summary(fit_bselmodel_fixed_approximate))
  # with V
  scenario_text("bselmodel-summary-mv-V-fixed-exact",       summary(fit_bselmodel.mv_V_fixed_exact))
  scenario_text("bselmodel-summary-mv-V-fixed-approximate", summary(fit_bselmodel.mv_V_fixed_approximate))

  # the rest has no direct comparison
  # 3lvl V vs vi
  scenario_text("bselmodel-summary-mv-vi-exact",         summary(fit_bselmodel.mv_exact))
  scenario_text("bselmodel-summary-mv-vi-approximate",   summary(fit_bselmodel.mv_approximate))
  scenario_text("bselmodel-summary-cluster-exact",       summary(fit_bselmodel_cluster_exact))
  scenario_text("bselmodel-summary-cluster-approximate", summary(fit_bselmodel_cluster_approximate))

  scenario_text("bselmodel-summary-mv-V-exact",          summary(fit_bselmodel.mv_V_exact))
  scenario_text("bselmodel-summary-mv-V-approximate",    summary(fit_bselmodel.mv_V_approximate))

  # study only
  scenario_text("bselmodel-summary-mv-vi-no-effect-exact",       summary(fit_bselmodel.mv_no_effect_exact))
  scenario_text("bselmodel-summary-mv-vi-no-effect-approximate", summary(fit_bselmodel.mv_no_effect_approximate))

  scenario_text("bselmodel-summary-mv-V-no-effect-exact",       summary(fit_bselmodel.mv_V_no_effect_exact))
  scenario_text("bselmodel-summary-mv-V-no-effect-approximate", summary(fit_bselmodel.mv_V_no_effect_approximate))

  # meta-reg
  scenario_text("bselmodel-summary-mv-vi-reg-exact",     summary(fit_bselmodel.mv_reg_exact))
  scenario_text("bselmodel-summary-cluster-reg-exact",   summary(fit_bselmodel_cluster_reg_exact))
  scenario_text("bselmodel-summary-mv-V-reg-exact",      summary(fit_bselmodel.mv_V_reg_exact))

  scenario_text("bselmodel-summary-mv-vi-reg-approximate",       summary(fit_bselmodel.mv_reg_approximate))
  scenario_text("bselmodel-summary-cluster-reg-approximate",     summary(fit_bselmodel_cluster_reg_approximate))
  scenario_text("bselmodel-summary-mv-V-reg-approximate",        summary(fit_bselmodel.mv_V_reg_approximate))
  # ABOVE checked, do not modify!

  scenario_text("bselmodel-metafor-comparison", data.frame(
    implementation        = c("metafor", "RoBMA exact", "RoBMA approximate"),
    mu                    = c(unname(fit_bselmodel_metafor[["beta"]][[1L]]) , ex_r(fit_bselmodel_exact, "mu"),             ex_r(fit_bselmodel_approximate, "mu") ),
    omega                 = c(unname(fit_bselmodel_metafor[["delta"]][[2L]]), ex_r(fit_bselmodel_exact, "omega[0.025,1]"), ex_r(fit_bselmodel_approximate, "omega[0.025,1]") ),
    total_random_sd       = c(sqrt(fit_bselmodel_metafor[["tau2"]]),          ex_r(fit_bselmodel_exact, "tau"),            ex_r(fit_bselmodel_approximate, "tau") ),
    row.names = NULL
  ))
  scenario_text("bselmodel-metafor-comparison-fixed", data.frame(
    implementation        = c("metafor", "RoBMA exact", "RoBMA approximate"),
    mu                    = c(unname(fit_bselmodel_metafor_fixed[["beta"]][[1L]]) , ex_r(fit_bselmodel.mv_fixed_exact, "mu"),             ex_r(fit_bselmodel_fixed_approximate, "mu")),
    omega                 = c(unname(fit_bselmodel_metafor_fixed[["delta"]][[2L]]), ex_r(fit_bselmodel.mv_fixed_exact, "omega[0.025,1]"), ex_r(fit_bselmodel_fixed_approximate, "omega[0.025,1]")),
    row.names = NULL
  ))

  ### Model-fit comparisons ----
  getloo_bselmodel <- function(fit) loo(fit)[["estimates"]]["looic", 1L]
  equivalent_bselmodel_fit <- function(
      fit_mv_V,    fit_mv_V_no_study,  fit_mv_V_no_effect,  fit_mv_V_fixed,
      fit_mv_vi,   fit_mv_vi_no_study, fit_mv_vi_no_effect, fit_mv_vi_fixed,
      fit_cluster, fit_simple, fit_fixed) {

    data.frame(
      structure = c("nested", "estimate", "study", "fixed"),
      logml_mv_V      = c(logml(fit_mv_V),    logml(fit_mv_V_no_study),  logml(fit_mv_V_no_effect),  logml(fit_mv_V_fixed)),
      logml_mv_vi     = c(logml(fit_mv_vi),   logml(fit_mv_vi_no_study), logml(fit_mv_vi_no_effect), logml(fit_mv_vi_fixed)),
      logml_bselmodel = c(logml(fit_cluster), logml(fit_simple),         NA,                         logml(fit_fixed)),
      looic_mv_V      = c(getloo_bselmodel(fit_mv_V),    getloo_bselmodel(fit_mv_V_no_study),  getloo_bselmodel(fit_mv_V_no_effect),  getloo_bselmodel(fit_mv_V_fixed)),
      looic_mv_vi     = c(getloo_bselmodel(fit_mv_vi),   getloo_bselmodel(fit_mv_vi_no_study), getloo_bselmodel(fit_mv_vi_no_effect), getloo_bselmodel(fit_mv_vi_fixed)),
      looic_bselmodel = c(getloo_bselmodel(fit_cluster), getloo_bselmodel(fit_simple),         NA,                                    getloo_bselmodel(fit_fixed)),
      row.names = NULL
    )
  }
  scenario_text("bselmodel-model-fit-equivalent-exact", equivalent_bselmodel_fit(
    fit_bselmodel.mv_V_exact, fit_bselmodel.mv_V_no_study_exact, fit_bselmodel.mv_V_no_effect_exact, fit_bselmodel.mv_V_fixed_exact,
    fit_bselmodel.mv_exact, fit_bselmodel.mv_no_study_exact, fit_bselmodel.mv_no_effect_exact, fit_bselmodel.mv_fixed_exact,
    fit_bselmodel_cluster_exact, fit_bselmodel_exact, fit_bselmodel_fixed_exact
  ))
  scenario_text("bselmodel-model-fit-equivalent-approximate", equivalent_bselmodel_fit(
    fit_bselmodel.mv_V_approximate, fit_bselmodel.mv_V_no_study_approximate, fit_bselmodel.mv_V_no_effect_approximate, fit_bselmodel.mv_V_fixed_approximate,
    fit_bselmodel.mv_approximate, fit_bselmodel.mv_no_study_approximate, fit_bselmodel.mv_no_effect_approximate, fit_bselmodel.mv_fixed_approximate,
    fit_bselmodel_cluster_approximate, fit_bselmodel_approximate, fit_bselmodel_fixed_approximate
  ))

  scenario_text("bselmodel-model-fit-reg-equivalent", data.frame(
    likelihood = c("exact", "approximate"),
    logml_mv_V  =     c(logml(fit_bselmodel.mv_V_reg_exact),    logml(fit_bselmodel.mv_V_reg_approximate)),
    logml_mv_vi =     c(logml(fit_bselmodel.mv_reg_exact),      logml(fit_bselmodel.mv_reg_approximate)),
    logml_bselmodel = c(logml(fit_bselmodel_cluster_reg_exact), logml(fit_bselmodel_cluster_reg_approximate)),
    row.names = NULL
  ))

  ### Basic fit plots ----
  scenario_plot("bselmodel-posterior-rho", {
    par(mfrow = c(1, 2))
    plot(fit_bselmodel_cluster_exact,   "rho", prior = TRUE, main = "exact", ylim = c(0, 4))
    lines(fit_bselmodel_cluster_exact,  "rho", lty = 2, density_method = "qCMDE")
    lines(fit_bselmodel.mv_exact, "var_prop(study)", col = "blue")
    lines(fit_bselmodel.mv_exact, "var_prop(study)", col = "blue", lty = 2, density_method = "qCMDE")
    lines(fit_bselmodel.mv_V_exact, "var_prop(study)", col = "red")
    lines(fit_bselmodel.mv_V_exact, "var_prop(study)", col = "red", lty = 2, density_method = "qCMDE")

    plot(fit_bselmodel_cluster_approximate,   "rho", prior = TRUE, main = "approximate", ylim = c(0, 4))
    lines(fit_bselmodel_cluster_approximate,  "rho", lty = 2, density_method = "qCMDE")
    lines(fit_bselmodel.mv_approximate, "var_prop(study)", col = "blue")
    lines(fit_bselmodel.mv_approximate, "var_prop(study)", col = "blue", lty = 2, density_method = "qCMDE")
    lines(fit_bselmodel.mv_V_approximate, "var_prop(study)", col = "red")
    lines(fit_bselmodel.mv_V_approximate, "var_prop(study)", col = "red", lty = 2, density_method = "qCMDE")
  })
  scenario_plot("bselmodel-posterior-location", {
    par(mfrow = c(1, 2))
    plot(fit_bselmodel_cluster_exact,   "mu", prior = TRUE, main = "exact", ylim = c(0, 4), xlim = c(-0.5, 1))
    lines(fit_bselmodel_cluster_exact,  "mu", lty = 2, density_method = "qCMDE")
    lines(fit_bselmodel.mv_exact, "mu", col = "blue")
    lines(fit_bselmodel.mv_exact, "mu", col = "blue", lty = 2, density_method = "qCMDE")
    lines(fit_bselmodel.mv_V_exact, "mu", col = "red")
    lines(fit_bselmodel.mv_V_exact, "mu", col = "red", lty = 2, density_method = "qCMDE")

    plot(fit_bselmodel_cluster_approximate,   "mu", prior = TRUE, main = "approximate", ylim = c(0, 4), xlim = c(-0.5, 1))
    lines(fit_bselmodel_cluster_approximate,  "mu", lty = 2, density_method = "qCMDE")
    lines(fit_bselmodel.mv_approximate, "mu", col = "blue")
    lines(fit_bselmodel.mv_approximate, "mu", col = "blue", lty = 2, density_method = "qCMDE")
    lines(fit_bselmodel.mv_V_approximate, "mu", col = "red")
    lines(fit_bselmodel.mv_V_approximate, "mu", col = "red", lty = 2, density_method = "qCMDE")
  })
  scenario_plot("bselmodel-posterior-mod", {
    par(mfrow = c(1, 2))
    plot(fit_bselmodel.mv_V_reg_exact,        "deltype", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 3), main = "exact")
    lines(fit_bselmodel.mv_V_reg_exact,       "deltype", density_method = "qCMDE", xlim = c(-1, 1), ylim = c(0, 3), main = "exact")

    plot(fit_bselmodel.mv_V_reg_approximate , "deltype", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 3), main = "approximate")
    lines(fit_bselmodel.mv_V_reg_approximate, "deltype", density_method = "qCMDE", xlim = c(-1, 1), ylim = c(0, 3), main = "approximate")
  })
  scenario_plot("bselmodel-posterior-random-mv-V", {
    par(mfrow = c(2, 3))
    plot(fit_bselmodel.mv_V_exact, "sd_total", prior = TRUE)
    lines(fit_bselmodel.mv_V_exact, "sd_total", density_method = "qCMDE")

    plot(fit_bselmodel.mv_V_exact, "study: sd", prior = TRUE)
    lines(fit_bselmodel.mv_V_exact, "study: sd", density_method = "qCMDE")

    plot(fit_bselmodel.mv_V_exact, "esid_study: sd", prior = TRUE)
    lines(fit_bselmodel.mv_V_exact, "esid_study: sd", density_method = "qCMDE")

    plot(fit_bselmodel.mv_V_exact, "var_prop(esid_study)", prior = TRUE)
    lines(fit_bselmodel.mv_V_exact, "var_prop(esid_study)", density_method = "qCMDE")

    plot(fit_bselmodel.mv_V_exact, "var_prop(study)", prior = TRUE)
    lines(fit_bselmodel.mv_V_exact, "var_prop(study)", density_method = "qCMDE")
  })
  scenario_plot("bselmodel-weightfunction", {
    par(mfrow = c(1, 2))
    plot_weightfunction(fit_bselmodel.mv_V_exact,       main = "exact")
    plot_weightfunction(fit_bselmodel.mv_V_approximate, main = "approximate")
  })

  ### Hypotheses ----
  set.seed(1)
  BF_bselmodel_rho_exact          <- scenario_time("BF_bselmodel_rho_exact",          hypothesis(fit_bselmodel_cluster_exact, c("rho != 0 vs rho = 0", "rho != 1 vs rho = 1"), seed = 1))
  BF_bselmodel_mv_rho_exact       <- scenario_time("BF_bselmodel_mv_rho_exact",       hypothesis(fit_bselmodel.mv_exact, c("var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1"), seed = 1))
  BF_bselmodel_rho_approximate    <- scenario_time("BF_bselmodel_rho_approximate",    hypothesis(fit_bselmodel_cluster_approximate, c("rho != 0 vs rho = 0", "rho != 1 vs rho = 1"), seed = 1))
  BF_bselmodel_mv_rho_approximate <- scenario_time("BF_bselmodel_mv_rho_approximate", hypothesis(fit_bselmodel.mv_approximate, c("var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1"), seed = 1))

  scenario_text("bselmodel-rho-bayes-factor-comparison", data.frame(
    likelihood = rep(c("exact", "approximate"), each = 2L),
    rho = rep(c(0, 1), 2L),
    qCMDE_bselmodel_BF            = c(BF_bselmodel_rho_exact[["BF"]],       BF_bselmodel_rho_approximate[["BF"]]),
    qCMDE_bselmodel_error_percent = c(BF_bselmodel_rho_exact[["BF_error"]], BF_bselmodel_rho_approximate[["BF_error"]]),
    qCMDE_mv_vi_BF             = c(BF_bselmodel_mv_rho_exact[["BF"]],       BF_bselmodel_mv_rho_approximate[["BF"]]),
    qCMDE_mv_vi_error_percent  = c(BF_bselmodel_mv_rho_exact[["BF_error"]], BF_bselmodel_mv_rho_approximate[["BF_error"]]),
    marglik_bselmodel_BF       = c(
      bf(fit_bselmodel_cluster_exact, fit_bselmodel_exact)[["bf"]], NA,
      bf(fit_bselmodel_cluster_approximate, fit_bselmodel_approximate)[["bf"]], NA),
    marglik_mv_vi_BF = c(
      bf(fit_bselmodel.mv_exact, fit_bselmodel.mv_no_study_exact)[["bf"]],
      bf(fit_bselmodel.mv_exact, fit_bselmodel.mv_no_effect_exact)[["bf"]],
      bf(fit_bselmodel.mv_approximate, fit_bselmodel.mv_no_study_approximate)[["bf"]],
      bf(fit_bselmodel.mv_approximate, fit_bselmodel.mv_no_effect_approximate)[["bf"]]
    ),
    row.names = NULL
  ))

  BF_bselmodel_mv_V_rho_exact       <- scenario_time("BF_bselmodel_mv_V_rho_exact",       hypothesis(fit_bselmodel.mv_V_exact, c("var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1"), seed = 1))
  BF_bselmodel_mv_V_rho_approximate <- scenario_time("BF_bselmodel_mv_V_rho_approximate", hypothesis(fit_bselmodel.mv_V_approximate, c("var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1"), seed = 1))

  scenario_text("bselmodel-rho-V-bayes-factor-comparison", data.frame(
    likelihood = rep(c("exact", "approximate"), each = 2L),
    rho = rep(c(0, 1), 2L),
    qCMDE_mv_V_BF             = c(BF_bselmodel_mv_V_rho_exact[["BF"]],             BF_bselmodel_mv_V_rho_exact[["BF"]]),
    qCMDE_mv_V_error_percent  = c(BF_bselmodel_mv_V_rho_approximate[["BF_error"]], BF_bselmodel_mv_V_rho_approximate[["BF_error"]]),
    marglik_mv_V_BF = c(
      bf(fit_bselmodel.mv_V_exact, fit_bselmodel.mv_V_no_study_exact)[["bf"]],
      bf(fit_bselmodel.mv_V_exact, fit_bselmodel.mv_V_no_effect_exact)[["bf"]],
      bf(fit_bselmodel.mv_V_approximate, fit_bselmodel.mv_V_no_study_approximate)[["bf"]],
      bf(fit_bselmodel.mv_V_approximate, fit_bselmodel.mv_V_no_effect_approximate)[["bf"]]
    ),
    row.names = NULL
  ))

  set.seed(1)
  scenario_text("bselmodel-mods-exact", scenario_time("BF_bselmodel_mods_exact", hypothesis(fit_bselmodel.mv_V_reg_exact, c("deltype[general] = 0 vs deltype[general] != 0", "deltype[general] = 0 vs deltype[general] > 0", "deltype[general] > 0 vs deltype[general] < 0"), seed = 1)))
  scenario_text("bselmodel-mods-approximate",  scenario_time("BF_bselmodel_mods_exact", hypothesis(fit_bselmodel.mv_V_reg_approximate, c("deltype[general] = 0 vs deltype[general] != 0", "deltype[general] = 0 vs deltype[general] > 0", "deltype[general] > 0 vs deltype[general] < 0"), seed = 1)))

  ### Pooled effects and predictions ----
  compare_bselmodel_pooled <- function(fit_RoBMA) {
    return(cbind.data.frame(
      "metafor" = t(data.frame(predict(fit_bselmodel_metafor))[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")]),
      "RoBMA"   = ex_p(fit_RoBMA)[, 1L]
    ))
  }
  scenario_text("bselmodel-pooled-effect-exact",       compare_bselmodel_pooled(fit_bselmodel_exact))
  scenario_text("bselmodel-pooled-effect-approximate", compare_bselmodel_pooled(fit_bselmodel_approximate))
  scenario_text("bselmodel-pooled-mv-V", cbind.data.frame(
    "exact"       = ex_p(fit_bselmodel.mv_V_exact)[, 1L],
    "approximate" = ex_p(fit_bselmodel.mv_V_approximate)[, 1L]
  ))

  compare_bselmodel_preds_reg <- function(fit_mv, fit_specialized = NULL, type = "terms") {
    cbind.data.frame(
      "bselmodel.mv" = unlist(data.frame(predict(fit_mv, type = type))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      "bselmodel"    = if (!is.null(fit_specialized)) { unlist(data.frame(predict(fit_specialized, type = type))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE)
      } else { rep(NA_real_, 9L)}
    )
  }
  scenario_text("bselmodel-predictions-mv-V-reg-exact", compare_bselmodel_preds_reg(fit_bselmodel.mv_V_reg_exact))
  scenario_text("bselmodel-predictions-mv-vi-reg-exact", compare_bselmodel_preds_reg(fit_bselmodel.mv_reg_exact, fit_bselmodel_cluster_reg_exact))
  scenario_text("bselmodel-predictions-mv-V-reg-pi-exact", compare_bselmodel_preds_reg(fit_bselmodel.mv_V_reg_exact, type = "estimate"))
  scenario_text("bselmodel-predictions-mv-vi-reg-pi-exact", compare_bselmodel_preds_reg(fit_bselmodel.mv_reg_exact, fit_bselmodel_cluster_reg_exact, type = "estimate"))
  scenario_text("bselmodel-predictions-mv-V-reg-approximate", compare_bselmodel_preds_reg(fit_bselmodel.mv_V_reg_approximate))
  scenario_text("bselmodel-predictions-mv-vi-reg-approximate", compare_bselmodel_preds_reg(fit_bselmodel.mv_reg_approximate, fit_bselmodel_cluster_reg_approximate))
  scenario_text("bselmodel-predictions-mv-V-reg-pi-approximate", compare_bselmodel_preds_reg(fit_bselmodel.mv_V_reg_approximate, type = "estimate"))
  scenario_text("bselmodel-predictions-mv-vi-reg-pi-approximate", compare_bselmodel_preds_reg(fit_bselmodel.mv_reg_approximate, fit_bselmodel_cluster_reg_approximate, type = "estimate"))

  ### Marginal means ----
  scenario_text("bselmodel-marginal-means-exact",       marginal_means(fit_bselmodel.mv_V_reg_exact))
  scenario_text("bselmodel-marginal-means-approximate", marginal_means(fit_bselmodel.mv_V_reg_approximate))
  scenario_plot("bselmodel-marginal-means-plot", {
    par(mfrow = c(1, 2))
    plot(marginal_means(fit_bselmodel.mv_V_reg_exact),       "deltype", prior = TRUE, xlim = c(-2, 2), main = "exact")
    plot(marginal_means(fit_bselmodel.mv_V_reg_exact),       "deltype", density_method = "qCMDE", lty = 3)

    plot(marginal_means(fit_bselmodel.mv_V_reg_approximate), "deltype", prior = TRUE, xlim = c(-2, 2), main = "approximate")
    plot(marginal_means(fit_bselmodel.mv_V_reg_approximate), "deltype", density_method = "qCMDE", lty = 3)
  })

  ### Heterogeneity ----
  scenario_text("bselmodel-summary-heterogeneity-mv-vi-exact",       summary_heterogeneity(fit_bselmodel.mv_exact))
  scenario_text("bselmodel-summary-heterogeneity-mv-vi-approximate", summary_heterogeneity(fit_bselmodel.mv_approximate))
  scenario_text("bselmodel-summary-heterogeneity-cluster-exact",       summary_heterogeneity(fit_bselmodel_cluster_exact))
  scenario_text("bselmodel-summary-heterogeneity-cluster-approximate", summary_heterogeneity(fit_bselmodel_cluster_approximate))
  scenario_text("bselmodel-summary-heterogeneity-mv-V-exact",       summary_heterogeneity(fit_bselmodel.mv_V_exact))
  scenario_text("bselmodel-summary-heterogeneity-mv-V-approximate", summary_heterogeneity(fit_bselmodel.mv_V_approximate))

  scenario_text("bselmodel-summary-heterogeneity-study-exact", summary_heterogeneity(fit_bselmodel.mv_V_no_effect_exact))
  scenario_text("bselmodel-summary-heterogeneity-study-approximate", summary_heterogeneity(fit_bselmodel.mv_V_no_effect_approximate))
  scenario_text("bselmodel-summary-heterogeneity-fixed-exact", summary_heterogeneity(fit_bselmodel.mv_V_fixed_exact))
  scenario_text("bselmodel-summary-heterogeneity-fixed-approximate", summary_heterogeneity(fit_bselmodel.mv_V_fixed_approximate))

  ### Random effects ----
  ranef_bselmodel.mv_V_exact       <- scenario_time("ranef_bselmodel.mv_V_exact",       ranef(fit_bselmodel.mv_V_exact))
  ranef_bselmodel.mv_V_approximate <- scenario_time("ranef_bselmodel.mv_V_approximate", ranef(fit_bselmodel.mv_V_approximate))
  ranef_bselmodel.mv_exact         <- scenario_time("ranef_bselmodel.mv_exact",         ranef(fit_bselmodel.mv_exact))
  ranef_bselmodel.mv_approximate      <- scenario_time("ranef_bselmodel.mv_approximate",      ranef(fit_bselmodel.mv_approximate))
  ranef_bselmodel_cluster_exact       <- scenario_time("ranef_bselmodel_cluster_exact",       ranef(fit_bselmodel_cluster_exact))
  ranef_bselmodel_cluster_approximate <- scenario_time("ranef_bselmodel_cluster_approximate", ranef(fit_bselmodel_cluster_approximate))

  plot_bselmodel_ranef_equivalence <- function(ranef_mv, ranef_cluster) {
    par(mfrow = c(1, 2))
    scenario_agreement_plot(as.data.frame(ranef_mv$study)[["Mean"]],      as.data.frame(ranef_cluster$cluster)[["Mean"]],  main = "study",      estimate_label = "mv", reference_label = "uni")
    scenario_agreement_plot(as.data.frame(ranef_mv$esid_study)[["Mean"]], as.data.frame(ranef_cluster$estimate)[["Mean"]], main = "esid_study", estimate_label = "mv", reference_label = "uni")
    return(invisible(NULL))
  }
  plot_bselmodel_ranef_equivalence2 <- function(ranef_mv, ranef_mv2) {
    par(mfrow = c(1, 2))
    scenario_agreement_plot(as.data.frame(ranef_mv$study)[["Mean"]],      as.data.frame(ranef_mv2$study)[["Mean"]],      main = "study",      estimate_label = "mv-V", reference_label = "mv-vi")
    scenario_agreement_plot(as.data.frame(ranef_mv$esid_study)[["Mean"]], as.data.frame(ranef_mv2$esid_study)[["Mean"]], main = "esid_study", estimate_label = "mv-V", reference_label = "mv-vi")
    return(invisible(NULL))
  }
  scenario_plot("bselmodel-ranef-mv-vi-exact",       plot_bselmodel_ranef_equivalence(ranef_bselmodel.mv_exact,       ranef_bselmodel_cluster_exact))
  scenario_plot("bselmodel-ranef-mv-vi-approximate", plot_bselmodel_ranef_equivalence(ranef_bselmodel.mv_approximate, ranef_bselmodel_cluster_approximate))
  scenario_plot("bselmodel-ranef-mv-V-vi",             plot_bselmodel_ranef_equivalence2(ranef_bselmodel.mv_V_exact,        ranef_bselmodel.mv_exact))
  scenario_plot("bselmodel-ranef-mv-V-vi-approximate", plot_bselmodel_ranef_equivalence2(ranef_bselmodel.mv_V_approximate,  ranef_bselmodel.mv_approximate))

  ### Diagnostics ----
  plot_bselmodel_diagnostic_equivalence <- function(fit_reference, fit_mv) {

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
      scenario_agreement_plot(reference_values[[diagnostic]], mv_values[[diagnostic]], main = diagnostic)
    }
    return(invisible(NULL))
  }
  scenario_plot("bselmodel-marginal-diagnostics-mv-vi-exact",       plot_bselmodel_diagnostic_equivalence(fit_bselmodel_cluster_exact,       fit_bselmodel.mv_exact))
  scenario_plot("bselmodel-marginal-diagnostics-mv-vi-approximate", plot_bselmodel_diagnostic_equivalence(fit_bselmodel_cluster_approximate, fit_bselmodel.mv_approximate))
  scenario_plot("bselmodel-marginal-diagnostics-mv-V-exact-vs-approx", plot_bselmodel_diagnostic_equivalence(fit_bselmodel.mv_V_exact,       fit_bselmodel.mv_V_approximate))

  ### Diagnostic plots ----
  scenario_plot("bselmodel-funnel-mv-V-exact", {
    par(mfrow = c(1, 2))
    funnel(fit_bselmodel.mv_V_exact, main = "funnel")
    bfunnel(fit_bselmodel.mv_V_exact, main = "bfunnel")
  })
  scenario_plot("bselmodel-funnel-mv-V-approximate", {
    par(mfrow = c(1, 2))
    funnel(fit_bselmodel.mv_V_approximate, main = "funnel")
    bfunnel(fit_bselmodel.mv_V_approximate, main = "bfunnel")
  })
  scenario_plot("bselmodel-qqnorm-mv-V", {
    par(mfrow = c(1, 2))
    qqnorm(fit_bselmodel.mv_V_exact, main = "exact")
    qqnorm(fit_bselmodel.mv_V_approximate, main = "approximate")
  })
  scenario_plot("bselmodel-zplot-mv-V", {
    par(mfrow = c(1, 2))
    zplot(fit_bselmodel.mv_V_exact, to = 10, main = "exact")
    zplot(fit_bselmodel.mv_V_approximate, to = 10, main = "approximate")
  })
  scenario_plot("bselmodel-funnel-reg", {
    par(mfrow = c(1, 2))
    funnel(fit_bselmodel.mv_V_reg_exact, main = "exact")
    funnel(fit_bselmodel.mv_V_reg_approximate, main = "approximate")
  })
  scenario_plot("bselmodel-qqnorm-reg", {
    par(mfrow = c(1, 2))
    qqnorm(fit_bselmodel.mv_V_reg_exact, main = "exact")
    qqnorm(fit_bselmodel.mv_V_reg_approximate, main = "approximate")
  })
  scenario_plot("bselmodel-zplot-reg", {
    par(mfrow = c(1, 2))
    zplot(fit_bselmodel.mv_V_reg_exact, to = 10, main = "exact")
    zplot(fit_bselmodel.mv_V_reg_approximate, to = 10, main = "approximate")
  })
})
