if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("assink2016")
# testthat::test_file("tests/scenarios/test-assink2016.R")

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
      ex_m(fit_metafor, metafor_parameters),      ex_r(fit_brma.mv, robma_parameters, component = robma_components),
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
  scenario_plot("fit.mv_diag_posterior_random", {
    par(mfrow = c(2, 3))

    plot(fit_brma.mv_diag, "sd_total", prior = TRUE)
    lines(fit_brma.mv_diag, "sd_total", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv_diag, "study: sd", prior = TRUE)
    lines(fit_brma.mv_diag, "study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_brma.mv_diag, "esid_study: sd", prior = TRUE)
    lines(fit_brma.mv_diag, "esid_study: sd", density_method = "qCMDE", lty = 2, density_control = list(samples = 1000L))

    plot(fit_brma.mv_diag, "var_prop(esid_study)", prior = TRUE)
    lines(fit_brma.mv_diag, "var_prop(esid_study)", density_method = "qCMDE", lty = 2)

    plot(fit_brma.mv_diag, "var_prop(study)", prior = TRUE)
    lines(fit_brma.mv_diag, "var_prop(study)", density_method = "qCMDE", lty = 2)
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

testthat::test_that("Assink bselmodel and bselmodel.mv models", {

  set.seed(1)
  data("dat.assink2016", package = "metadat")

  V_assink <- metafor::vcalc(
    vi, cluster = study, type = deltype, obs = esid,
    rho = c(0.7, 0.5), data = dat.assink2016
  )
  sampling_blocks <- split(
    seq_len(nrow(dat.assink2016)),
    dat.assink2016[["study"]]
  )
  sampling_factor_count <- sum(vapply(
    sampling_blocks,
    function(rows) length(unique(dat.assink2016[["deltype"]][rows])),
    integer(1L)
  ))
  sampling_loading <- matrix(
    0,
    nrow = nrow(dat.assink2016),
    ncol = sampling_factor_count
  )
  sampling_sei <- sqrt(dat.assink2016[["vi"]])
  column_start <- 1L
  for (rows in sampling_blocks) {
    types <- unique(dat.assink2016[["deltype"]][rows])
    type_covariance <- matrix(0.5, nrow = length(types), ncol = length(types))
    diag(type_covariance) <- 0.7
    type_loading <- t(chol(type_covariance))
    columns <- column_start:(column_start + length(types) - 1L)
    sampling_loading[rows, columns] <- sampling_sei[rows] *
      type_loading[
        match(dat.assink2016[["deltype"]][rows], types),
        , drop = FALSE
      ]
    column_start <- max(columns) + 1L
  }
  V_assink_factor <- known_v_factor(
    diagonal = 0.3 * dat.assink2016[["vi"]],
    loading  = sampling_loading
  )
  reconstructed_V_assink <-
    diag(V_assink_factor[["diagonal"]]) +
    tcrossprod(V_assink_factor[["loading"]])
  scenario_text("bselmodel-known-v-factor-reconstruction", data.frame(
    reconstructs_V = isTRUE(all.equal(
      unclass(V_assink), reconstructed_V_assink, tolerance = 1e-14
    )),
    row.names = NULL
  ))

  fit_bselmodel_metafor <- metafor::selmodel(
    metafor::rma(yi, vi, data = dat.assink2016, method = "ML"),
    type = "stepfun", steps = 0.025, decreasing = TRUE
  )

  ### Exact model fits ----
  fit_bselmodel_exact <- scenario_fit("fit_bselmodel_exact", {
    tmp <- bselmodel(yi = yi, vi = vi, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", seed = 1)
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
    tmp <- bselmodel(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", seed = 1)
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
    tmp <- bselmodel.mv(yi = yi, V = V_assink_factor, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", sample = 3000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_reg_exact <- scenario_fit("fit_bselmodel.mv_V_reg_exact", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink_factor, mods = ~ deltype, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", sample = 4000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_no_study_exact <- scenario_fit("fit_bselmodel.mv_V_no_study_exact", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink_factor, random = ~ 1 | study:esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", sample = 2000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_no_effect_exact <- scenario_fit("fit_bselmodel.mv_V_no_effect_exact", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink_factor, random = ~ 1 | study, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", sample = 3000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_fixed_exact <- scenario_fit("fit_bselmodel.mv_V_fixed_exact", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink_factor, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "exact", sample = 1500, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })

  ### Approximate model fits ----
  fit_bselmodel_approximate <- scenario_fit("fit_bselmodel_approximate", {
    tmp <- bselmodel(yi = yi, vi = vi, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_cluster_approximate <- scenario_fit("fit_bselmodel_cluster_approximate", {
    tmp <- bselmodel(yi = yi, vi = vi, cluster = study, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_cluster_reg_approximate <- scenario_fit("fit_bselmodel_cluster_reg_approximate", {
    tmp <- bselmodel(yi = yi, vi = vi, mods = ~ deltype, cluster = study, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", sample = 15000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel_fixed_approximate <- scenario_fit("fit_bselmodel_fixed_approximate", {
    tmp <- bselmodel(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", seed = 1)
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
    tmp <- bselmodel.mv(yi = yi, vi = vi, mods = ~ deltype, random = ~ 1 | study / esid, prior_heterogeneity = BayesTools::prior_random(study = BayesTools::random_block(parameterization = "centered")), measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", sample = 15000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_no_study_approximate <- scenario_fit("fit_bselmodel.mv_no_study_approximate", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study:esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_no_effect_approximate <- scenario_fit("fit_bselmodel.mv_no_effect_approximate", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | study, prior_heterogeneity = BayesTools::prior_random(parameterization = "centered"), measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", sample = 35000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_fixed_approximate <- scenario_fit("fit_bselmodel.mv_fixed_approximate", {
    tmp <- bselmodel.mv(yi = yi, vi = vi, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", seed = 1)
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
    tmp <- bselmodel.mv(yi = yi, V = V_assink, mods = ~ deltype, random = ~ 1 | study / esid, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", sample = 15000, parallel = TRUE, seed = 1)
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
    tmp <- bselmodel.mv(yi = yi, V = V_assink, random = ~ 1 | study, prior_heterogeneity = BayesTools::prior_random(parameterization = "centered"), measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", sample = 180000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_bselmodel.mv_V_fixed_approximate <- scenario_fit("fit_bselmodel.mv_V_fixed_approximate", {
    tmp <- bselmodel.mv(yi = yi, V = V_assink, prior_heterogeneity = NULL, measure = "SMD", data = dat.assink2016, selection_likelihood = "approximate", sample = 15000, parallel = TRUE, seed = 1)
    tmp <- add_marglik(tmp)
    tmp <- add_loo(tmp)
    return(tmp)
  })

  bselmodel_fit_names <- sort(ls(
    pattern = "^fit_bselmodel", envir = environment()
  ))
  bselmodel_fits <- mget(
    bselmodel_fit_names, envir = environment(), inherits = FALSE
  )
  bselmodel_fits <- bselmodel_fits[vapply(
    bselmodel_fits, inherits, logical(1L), what = "brma"
  )]
  scenario_text("bselmodel-fit-validity", do.call(rbind, lapply(
    names(bselmodel_fits),
    function(name) {

      fit <- bselmodel_fits[[name]]
      diagnostics <- attr(fit[["fit"]][["converged"]], "diagnostics")
      assessable  <- diagnostics[diagnostics[["state"]] == "assessable", ]
      loo_result  <- loo(fit)
      pareto_k    <- loo_result[["diagnostics"]][["pareto_k"]]
      data.frame(
        fit                       = name,
        likelihood                = fit[["selection_likelihood"]][["type"]],
        target                    = fit[["selection_likelihood"]][["target"]],
        converged                 = isTRUE(fit[["fit"]][["converged"]]),
        max_Rhat                  = max(assessable[["Rhat"]]),
        min_ESS                   = min(assessable[["ESS"]]),
        marginal_likelihood_finite = is.finite(logml(fit)),
        loo_finite                = all(is.finite(loo_result[["estimates"]][, 1L])),
        max_pareto_k              = max(pareto_k),
        pareto_k_above_0.7        = sum(pareto_k > 0.7),
        pareto_k_above_1          = sum(pareto_k > 1),
        row.names                 = NULL
      )
    }
  )))

  ### Model summaries ----
  scenario_text("bselmodel-summary-metafor", fit_bselmodel_metafor)

  scenario_text("bselmodel-summary-mv-V-exact",          summary(fit_bselmodel.mv_V_exact))
  scenario_text("bselmodel-summary-mv-V-no-study-exact", summary(fit_bselmodel.mv_V_no_study_exact))
  scenario_text("bselmodel-summary-mv-V-no-effect-exact", summary(fit_bselmodel.mv_V_no_effect_exact))
  scenario_text("bselmodel-summary-mv-V-fixed-exact",     summary(fit_bselmodel.mv_V_fixed_exact))
  scenario_text("bselmodel-summary-mv-V-reg-exact",       summary(fit_bselmodel.mv_V_reg_exact))
  scenario_text("bselmodel-summary-mv-vi-exact",          summary(fit_bselmodel.mv_exact))
  scenario_text("bselmodel-summary-mv-vi-no-study-exact", summary(fit_bselmodel.mv_no_study_exact))
  scenario_text("bselmodel-summary-mv-vi-no-effect-exact", summary(fit_bselmodel.mv_no_effect_exact))
  scenario_text("bselmodel-summary-mv-vi-fixed-exact",     summary(fit_bselmodel.mv_fixed_exact))
  scenario_text("bselmodel-summary-mv-vi-reg-exact",       summary(fit_bselmodel.mv_reg_exact))
  scenario_text("bselmodel-summary-simple-exact",          summary(fit_bselmodel_exact))
  scenario_text("bselmodel-summary-cluster-exact",         summary(fit_bselmodel_cluster_exact))
  scenario_text("bselmodel-summary-cluster-reg-exact",     summary(fit_bselmodel_cluster_reg_exact))
  scenario_text("bselmodel-summary-fixed-exact",           summary(fit_bselmodel_fixed_exact))

  scenario_text("bselmodel-summary-mv-V-approximate",          summary(fit_bselmodel.mv_V_approximate))
  scenario_text("bselmodel-summary-mv-V-no-study-approximate", summary(fit_bselmodel.mv_V_no_study_approximate))
  scenario_text("bselmodel-summary-mv-V-no-effect-approximate", summary(fit_bselmodel.mv_V_no_effect_approximate))
  scenario_text("bselmodel-summary-mv-V-fixed-approximate",     summary(fit_bselmodel.mv_V_fixed_approximate))
  scenario_text("bselmodel-summary-mv-V-reg-approximate",       summary(fit_bselmodel.mv_V_reg_approximate))
  scenario_text("bselmodel-summary-mv-vi-approximate",          summary(fit_bselmodel.mv_approximate))
  scenario_text("bselmodel-summary-mv-vi-no-study-approximate", summary(fit_bselmodel.mv_no_study_approximate))
  scenario_text("bselmodel-summary-mv-vi-no-effect-approximate", summary(fit_bselmodel.mv_no_effect_approximate))
  scenario_text("bselmodel-summary-mv-vi-fixed-approximate",     summary(fit_bselmodel.mv_fixed_approximate))
  scenario_text("bselmodel-summary-mv-vi-reg-approximate",       summary(fit_bselmodel.mv_reg_approximate))
  scenario_text("bselmodel-summary-simple-approximate",          summary(fit_bselmodel_approximate))
  scenario_text("bselmodel-summary-cluster-approximate",         summary(fit_bselmodel_cluster_approximate))
  scenario_text("bselmodel-summary-cluster-reg-approximate",     summary(fit_bselmodel_cluster_reg_approximate))
  scenario_text("bselmodel-summary-fixed-approximate",           summary(fit_bselmodel_fixed_approximate))

  # metafor::selmodel() has no rma.mv method, so only the shared univariate
  # selection-model target has a direct metafor comparison.
  scenario_text("bselmodel-metafor-comparison", data.frame(
    implementation        = c("metafor", "RoBMA exact", "RoBMA approximate"),
    mu                    = c(
      unname(fit_bselmodel_metafor[["beta"]][[1L]]),
      ex_r(fit_bselmodel_exact, "mu"),
      ex_r(fit_bselmodel_approximate, "mu")
    ),
    omega                 = c(
      unname(fit_bselmodel_metafor[["delta"]][[2L]]),
      ex_r(fit_bselmodel_exact, "omega[0.025,1]"),
      ex_r(fit_bselmodel_approximate, "omega[0.025,1]")
    ),
    total_random_sd       = c(
      sqrt(fit_bselmodel_metafor[["tau2"]]),
      ex_r(fit_bselmodel_exact, "tau"),
      ex_r(fit_bselmodel_approximate, "tau")
    ),
    row.names = NULL
  ))

  scenario_text("bselmodel-selection-targets", data.frame(
    model = rep(c("simple", "cluster", "mv_vi", "mv_correlated_V"), 2L),
    likelihood = rep(c("exact", "approximate"), each = 4L),
    target = vapply(list(
      fit_bselmodel_exact, fit_bselmodel_cluster_exact,
      fit_bselmodel.mv_exact, fit_bselmodel.mv_V_exact,
      fit_bselmodel_approximate, fit_bselmodel_cluster_approximate,
      fit_bselmodel.mv_approximate, fit_bselmodel.mv_V_approximate
    ), function(fit) fit[["selection_likelihood"]][["target"]], character(1L)),
    row.names = NULL
  ))

  ### Model-fit comparisons ----
  getloo_bselmodel <- function(fit) loo(fit)[["estimates"]]["looic", 1L]
  equivalent_bselmodel_fit <- function(
      fit_mv_V, fit_mv_V_no_study, fit_mv_V_no_effect, fit_mv_V_fixed,
      fit_mv_vi, fit_mv_vi_no_study, fit_mv_vi_no_effect, fit_mv_vi_fixed,
      fit_cluster, fit_simple, fit_fixed) {

    data.frame(
      structure = c("nested", "estimate", "study", "fixed"),
      logml_mv_V = c(logml(fit_mv_V), logml(fit_mv_V_no_study), logml(fit_mv_V_no_effect), logml(fit_mv_V_fixed)),
      logml_mv_vi = c(logml(fit_mv_vi), logml(fit_mv_vi_no_study), logml(fit_mv_vi_no_effect), logml(fit_mv_vi_fixed)),
      logml_bselmodel = c(logml(fit_cluster), logml(fit_simple), NA, logml(fit_fixed)),
      looic_mv_V = c(getloo_bselmodel(fit_mv_V), getloo_bselmodel(fit_mv_V_no_study), getloo_bselmodel(fit_mv_V_no_effect), getloo_bselmodel(fit_mv_V_fixed)),
      looic_mv_vi = c(getloo_bselmodel(fit_mv_vi), getloo_bselmodel(fit_mv_vi_no_study), getloo_bselmodel(fit_mv_vi_no_effect), getloo_bselmodel(fit_mv_vi_fixed)),
      looic_bselmodel = c(getloo_bselmodel(fit_cluster), getloo_bselmodel(fit_simple), NA, getloo_bselmodel(fit_fixed)),
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

  scenario_text("bselmodel-compare-mv-V-loo-exact", loo_model_weights(
    fit_bselmodel.mv_V_exact, fit_bselmodel.mv_V_no_study_exact, fit_bselmodel.mv_V_no_effect_exact, fit_bselmodel.mv_V_fixed_exact
  ))
  scenario_text("bselmodel-compare-mv-V-logml-exact", t(t(round(post_prob(
    fit_bselmodel.mv_V_exact, fit_bselmodel.mv_V_no_study_exact, fit_bselmodel.mv_V_no_effect_exact, fit_bselmodel.mv_V_fixed_exact
  ), 3))))
  scenario_text("bselmodel-compare-mv-vi-loo-exact", loo_model_weights(
    fit_bselmodel.mv_exact, fit_bselmodel.mv_no_study_exact, fit_bselmodel.mv_no_effect_exact, fit_bselmodel.mv_fixed_exact
  ))
  scenario_text("bselmodel-compare-mv-vi-logml-exact", t(t(round(post_prob(
    fit_bselmodel.mv_exact, fit_bselmodel.mv_no_study_exact, fit_bselmodel.mv_no_effect_exact, fit_bselmodel.mv_fixed_exact
  ), 3))))
  scenario_text("bselmodel-compare-specialized-loo-exact", loo_model_weights(
    fit_bselmodel_cluster_exact, fit_bselmodel_exact, fit_bselmodel_fixed_exact
  ))
  scenario_text("bselmodel-compare-specialized-logml-exact", t(t(round(post_prob(
    fit_bselmodel_cluster_exact, fit_bselmodel_exact, fit_bselmodel_fixed_exact
  ), 3))))

  scenario_text("bselmodel-compare-mv-V-loo-approximate", loo_model_weights(
    fit_bselmodel.mv_V_approximate, fit_bselmodel.mv_V_no_study_approximate, fit_bselmodel.mv_V_no_effect_approximate, fit_bselmodel.mv_V_fixed_approximate
  ))
  scenario_text("bselmodel-compare-mv-V-logml-approximate", t(t(round(post_prob(
    fit_bselmodel.mv_V_approximate, fit_bselmodel.mv_V_no_study_approximate, fit_bselmodel.mv_V_no_effect_approximate, fit_bselmodel.mv_V_fixed_approximate
  ), 3))))
  scenario_text("bselmodel-compare-mv-vi-loo-approximate", loo_model_weights(
    fit_bselmodel.mv_approximate, fit_bselmodel.mv_no_study_approximate, fit_bselmodel.mv_no_effect_approximate, fit_bselmodel.mv_fixed_approximate
  ))
  scenario_text("bselmodel-compare-mv-vi-logml-approximate", t(t(round(post_prob(
    fit_bselmodel.mv_approximate, fit_bselmodel.mv_no_study_approximate, fit_bselmodel.mv_no_effect_approximate, fit_bselmodel.mv_fixed_approximate
  ), 3))))
  scenario_text("bselmodel-compare-specialized-loo-approximate", loo_model_weights(
    fit_bselmodel_cluster_approximate, fit_bselmodel_approximate, fit_bselmodel_fixed_approximate
  ))
  scenario_text("bselmodel-compare-specialized-logml-approximate", t(t(round(post_prob(
    fit_bselmodel_cluster_approximate, fit_bselmodel_approximate, fit_bselmodel_fixed_approximate
  ), 3))))

  scenario_text("bselmodel-model-fit-reg-equivalent", data.frame(
    likelihood = c("exact", "approximate"),
    logml_mv_V = c(logml(fit_bselmodel.mv_V_reg_exact), logml(fit_bselmodel.mv_V_reg_approximate)),
    logml_mv_vi = c(logml(fit_bselmodel.mv_reg_exact), logml(fit_bselmodel.mv_reg_approximate)),
    logml_bselmodel = c(logml(fit_bselmodel_cluster_reg_exact), logml(fit_bselmodel_cluster_reg_approximate)),
    row.names = NULL
  ))

  ### Basic fit plots ----
  scenario_plot("bselmodel-posterior-rho", {
    par(mfrow = c(1, 2))
    plot(fit_bselmodel_cluster_exact, "rho", prior = TRUE, main = "exact")
    lines(fit_bselmodel.mv_exact, "var_prop(study)", col = "blue")
    plot(fit_bselmodel_cluster_approximate, "rho", prior = TRUE, main = "approximate")
    lines(fit_bselmodel.mv_approximate, "var_prop(study)", col = "blue")
  })
  scenario_plot("bselmodel-posterior-location", {
    par(mfrow = c(1, 2))
    plot(fit_bselmodel.mv_V_exact, "mu", prior = TRUE, xlim = c(-1, 1), main = "exact")
    lines(fit_bselmodel.mv_exact, "mu", col = "blue")
    plot(fit_bselmodel.mv_V_approximate, "mu", prior = TRUE, xlim = c(-1, 1), main = "approximate")
    lines(fit_bselmodel.mv_approximate, "mu", col = "blue")
  })
  scenario_plot("bselmodel-posterior-mod", {
    par(mfrow = c(1, 2))
    plot(fit_bselmodel.mv_V_reg_exact, "deltype", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 3), main = "exact")
    plot(fit_bselmodel.mv_V_reg_approximate, "deltype", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 3), main = "approximate")
  })
  scenario_plot("bselmodel-posterior-random-mv-V", {
    par(mfrow = c(2, 3))
    plot(fit_bselmodel.mv_V_exact, "sd_total", prior = TRUE)
    plot(fit_bselmodel.mv_V_exact, "study: sd", prior = TRUE)
    plot(fit_bselmodel.mv_V_exact, "esid_study: sd", prior = TRUE)
    plot(fit_bselmodel.mv_V_exact, "var_prop(esid_study)", prior = TRUE)
    plot(fit_bselmodel.mv_V_exact, "var_prop(study)", prior = TRUE)
  })
  scenario_plot("bselmodel-weightfunction", {
    par(mfrow = c(1, 2))
    plot_weightfunction(fit_bselmodel.mv_V_exact, main = "exact")
    plot_weightfunction(fit_bselmodel.mv_V_approximate, main = "approximate")
  })

  ### Hypotheses ----
  # A fixed regression-sized design keeps exact correlated-V ordinate checks
  # tractable while retaining their Monte Carlo diagnostics in the snapshots.
  bselmodel_density_control <- list(
    samples              = 100L,
    normalization_points = 20L
  )
  set.seed(1)
  BF_bselmodel_rho_exact <- scenario_time("BF_bselmodel_rho_exact", hypothesis(
    fit_bselmodel_cluster_exact,
    c("rho != 0 vs rho = 0", "rho != 1 vs rho = 1"),
    seed = 1,
    density_control = bselmodel_density_control
  ))
  BF_bselmodel_mv_rho_exact <- scenario_time("BF_bselmodel_mv_rho_exact", hypothesis(
    fit_bselmodel.mv_exact,
    c("var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1"),
    seed = 1,
    density_control = bselmodel_density_control
  ))
  BF_bselmodel_rho_approximate <- scenario_time("BF_bselmodel_rho_approximate", hypothesis(
    fit_bselmodel_cluster_approximate,
    c("rho != 0 vs rho = 0", "rho != 1 vs rho = 1"),
    seed = 1,
    density_control = bselmodel_density_control
  ))
  BF_bselmodel_mv_rho_approximate <- scenario_time("BF_bselmodel_mv_rho_approximate", hypothesis(
    fit_bselmodel.mv_approximate,
    c("var_prop(study) != 0 vs var_prop(study) = 0", "var_prop(study) != 1 vs var_prop(study) = 1"),
    seed = 1,
    density_control = bselmodel_density_control
  ))
  scenario_text("bselmodel-rho-bayes-factor-comparison", data.frame(
    likelihood = rep(c("exact", "approximate"), each = 2L),
    rho = rep(c(0, 1), 2L),
    qCMDE_bselmodel_BF = c(BF_bselmodel_rho_exact[["BF"]], BF_bselmodel_rho_approximate[["BF"]]),
    qCMDE_bselmodel_error_percent = c(BF_bselmodel_rho_exact[["BF_error"]], BF_bselmodel_rho_approximate[["BF_error"]]),
    qCMDE_mv_vi_BF = c(BF_bselmodel_mv_rho_exact[["BF"]], BF_bselmodel_mv_rho_approximate[["BF"]]),
    qCMDE_mv_vi_error_percent = c(BF_bselmodel_mv_rho_exact[["BF_error"]], BF_bselmodel_mv_rho_approximate[["BF_error"]]),
    marglik_bselmodel_BF = c(
      bf(fit_bselmodel_cluster_exact, fit_bselmodel_exact)[["bf"]], NA,
      bf(fit_bselmodel_cluster_approximate, fit_bselmodel_approximate)[["bf"]], NA
    ),
    marglik_mv_vi_BF = c(
      bf(fit_bselmodel.mv_exact, fit_bselmodel.mv_no_study_exact)[["bf"]],
      bf(fit_bselmodel.mv_exact, fit_bselmodel.mv_no_effect_exact)[["bf"]],
      bf(fit_bselmodel.mv_approximate, fit_bselmodel.mv_no_study_approximate)[["bf"]],
      bf(fit_bselmodel.mv_approximate, fit_bselmodel.mv_no_effect_approximate)[["bf"]]
    ),
    row.names = NULL
  ))

  set.seed(1)
  BF_bselmodel_random_allocation_exact <- scenario_time("BF_bselmodel_random_allocation_exact", hypothesis(
    fit_bselmodel.mv_V_exact,
    c(
      "var_prop(study) != 0 vs var_prop(study) = 0",
      "var_prop(study) != 1 vs var_prop(study) = 1"
    ),
    seed = 1,
    density_control = bselmodel_density_control
  ))
  set.seed(1)
  BF_bselmodel_random_sd_exact <- scenario_time("BF_bselmodel_random_sd_exact", hypothesis(
    fit_bselmodel.mv_V_exact,
    "sd_total = 0",
    seed = 1,
    density_control = bselmodel_density_control
  ))
  set.seed(1)
  BF_bselmodel_random_allocation_approximate <- scenario_time("BF_bselmodel_random_allocation_approximate", hypothesis(
    fit_bselmodel.mv_V_approximate,
    "var_prop(study) != 0 vs var_prop(study) = 0",
    seed = 1,
    density_control = bselmodel_density_control
  ))
  set.seed(1)
  BF_bselmodel_random_approximate_boundary_error <- testthat::expect_error(
    hypothesis(
      fit_bselmodel.mv_V_approximate,
      "var_prop(study) != 1 vs var_prop(study) = 1",
      seed = 1,
      density_control = bselmodel_density_control
    ),
    class = "RoBMA_density_ordinate_error"
  )
  testthat::expect_match(
    conditionMessage(BF_bselmodel_random_approximate_boundary_error),
    paste0(
      "qCMDE posterior ordinate for 'var_prop\\(study\\) = 1' was ",
      "rejected by diagnostics"
    )
  )
  set.seed(1)
  BF_bselmodel_random_approximate_sd_error <- testthat::expect_error(
    hypothesis(
      fit_bselmodel.mv_V_approximate,
      "sd_total = 0",
      seed = 1,
      density_control = bselmodel_density_control
    ),
    class = "RoBMA_density_ordinate_error"
  )
  testthat::expect_match(
    conditionMessage(BF_bselmodel_random_approximate_sd_error),
    paste0(
      "qCMDE posterior ordinate for 'sd_total = 0' was rejected by ",
      "diagnostics: posterior ordinate is zero or non-finite"
    )
  )
  exact_random_BF <- c(
    BF_bselmodel_random_allocation_exact[["BF"]],
    BF_bselmodel_random_sd_exact[["BF"]]
  )
  exact_random_error <- c(
    BF_bselmodel_random_allocation_exact[["BF_error"]],
    BF_bselmodel_random_sd_exact[["BF_error"]]
  )
  approximate_random_BF <- c(
    BF_bselmodel_random_allocation_approximate[["BF"]],
    NA_real_,
    NA_real_
  )
  approximate_random_error <- c(
    BF_bselmodel_random_allocation_approximate[["BF_error"]],
    NA_real_,
    NA_real_
  )
  scenario_text("bselmodel-random-bayes-factor-comparison", data.frame(
    likelihood = rep(c("exact", "approximate"), each = 3L),
    hypothesis = rep(c("rho != 0", "rho != 1", "sd != 0"), 2L),
    qCMDE_status = c(
      rep("available", 3L),
      "available", "rejected by diagnostics", "rejected by diagnostics"
    ),
    qCMDE_BF = c(exact_random_BF, approximate_random_BF),
    qCMDE_error_percent = c(exact_random_error, approximate_random_error),
    marglik_BF = c(
      bf(fit_bselmodel.mv_V_exact, fit_bselmodel.mv_V_no_study_exact)[["bf"]],
      bf(fit_bselmodel.mv_V_exact, fit_bselmodel.mv_V_no_effect_exact)[["bf"]],
      bf(fit_bselmodel.mv_V_exact, fit_bselmodel.mv_V_fixed_exact)[["bf"]],
      bf(fit_bselmodel.mv_V_approximate, fit_bselmodel.mv_V_no_study_approximate)[["bf"]],
      bf(fit_bselmodel.mv_V_approximate, fit_bselmodel.mv_V_no_effect_approximate)[["bf"]],
      bf(fit_bselmodel.mv_V_approximate, fit_bselmodel.mv_V_fixed_approximate)[["bf"]]
    ),
    row.names = NULL
  ))

  set.seed(1)
  scenario_text("bselmodel-mods-exact", scenario_time("BF_bselmodel_mods_exact", hypothesis(
    fit_bselmodel.mv_V_reg_exact,
    c(
      "deltype[general] = 0 vs deltype[general] != 0",
      "deltype[general] = 0 vs deltype[general] > 0",
      "deltype[general] > 0 vs deltype[general] < 0"
    ),
    seed = 1,
    density_control = bselmodel_density_control
  )))
  scenario_text("bselmodel-mods-approximate", scenario_time("BF_bselmodel_mods_approximate", hypothesis(
    fit_bselmodel.mv_V_reg_approximate,
    c(
      "deltype[general] = 0 vs deltype[general] != 0",
      "deltype[general] = 0 vs deltype[general] > 0",
      "deltype[general] > 0 vs deltype[general] < 0"
    ),
    seed = 1,
    density_control = bselmodel_density_control
  )))

  ### Pooled effects and predictions ----
  compare_bselmodel_pooled <- function(fit_RoBMA) {

    return(cbind.data.frame(
      "metafor" = t(data.frame(predict(fit_bselmodel_metafor))[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")]),
      "RoBMA"   = ex_p(fit_RoBMA)[, 1L]
    ))
  }
  scenario_text("bselmodel-pooled-effect-exact",       compare_bselmodel_pooled(fit_bselmodel_exact))
  scenario_text("bselmodel-pooled-effect-approximate", compare_bselmodel_pooled(fit_bselmodel_approximate))

  equivalent_bselmodel_pooled <- function(fit_cluster, fit_mv, fit_simple, fit_mv_no_study, fit_fixed, fit_mv_fixed) {

    cbind.data.frame(
      "cluster"      = ex_p(fit_cluster)[, 1L],
      "mv_nested"    = ex_p(fit_mv)[, 1L],
      "simple"       = ex_p(fit_simple)[, 1L],
      "mv_estimate"  = ex_p(fit_mv_no_study)[, 1L],
      "fixed"        = ex_p(fit_fixed)[, 1L],
      "mv_fixed"     = ex_p(fit_mv_fixed)[, 1L]
    )
  }
  scenario_text("bselmodel-pooled-equivalent-exact", equivalent_bselmodel_pooled(
    fit_bselmodel_cluster_exact, fit_bselmodel.mv_exact,
    fit_bselmodel_exact, fit_bselmodel.mv_no_study_exact,
    fit_bselmodel_fixed_exact, fit_bselmodel.mv_fixed_exact
  ))
  scenario_text("bselmodel-pooled-equivalent-approximate", equivalent_bselmodel_pooled(
    fit_bselmodel_cluster_approximate, fit_bselmodel.mv_approximate,
    fit_bselmodel_approximate, fit_bselmodel.mv_no_study_approximate,
    fit_bselmodel_fixed_approximate, fit_bselmodel.mv_fixed_approximate
  ))
  scenario_text("bselmodel-pooled-mv-V", cbind.data.frame(
    "exact"       = ex_p(fit_bselmodel.mv_V_exact)[, 1L],
    "approximate" = ex_p(fit_bselmodel.mv_V_approximate)[, 1L]
  ))

  compare_bselmodel_preds_reg <- function(fit_mv, fit_specialized = NULL, type = "terms") {

    cbind.data.frame(
      "bselmodel.mv" = unlist(data.frame(predict(fit_mv, type = type))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE),
      "bselmodel" = if (!is.null(fit_specialized)) {
        unlist(data.frame(predict(fit_specialized, type = type))[c(1, 50, 80), c("Mean", "CI_0.025", "CI_0.975")], use.names = FALSE)
      } else {
        rep(NA_real_, 9L)
      }
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
  scenario_text("bselmodel-marginal-means-exact", marginal_means(fit_bselmodel.mv_V_reg_exact))
  scenario_text("bselmodel-marginal-means-approximate", marginal_means(fit_bselmodel.mv_V_reg_approximate))
  scenario_plot("bselmodel-marginal-means-plot", {
    par(mfrow = c(1, 2))
    plot(marginal_means(fit_bselmodel.mv_V_reg_exact), "deltype", prior = TRUE, xlim = c(-2, 2), main = "exact")
    plot(marginal_means(fit_bselmodel.mv_V_reg_approximate), "deltype", prior = TRUE, xlim = c(-2, 2), main = "approximate")
  })

  ### Heterogeneity ----
  scenario_text("bselmodel-summary-heterogeneity-mv-vi-exact", summary_heterogeneity(fit_bselmodel.mv_exact))
  scenario_text("bselmodel-summary-heterogeneity-cluster-exact", summary_heterogeneity(fit_bselmodel_cluster_exact))
  scenario_text("bselmodel-summary-heterogeneity-mv-V-exact", summary_heterogeneity(fit_bselmodel.mv_V_exact))
  scenario_text("bselmodel-summary-heterogeneity-study-exact", summary_heterogeneity(fit_bselmodel.mv_V_no_effect_exact))
  scenario_text("bselmodel-summary-heterogeneity-fixed-exact", summary_heterogeneity(fit_bselmodel.mv_V_fixed_exact))
  scenario_text("bselmodel-summary-heterogeneity-mv-vi-approximate", summary_heterogeneity(fit_bselmodel.mv_approximate))
  scenario_text("bselmodel-summary-heterogeneity-cluster-approximate", summary_heterogeneity(fit_bselmodel_cluster_approximate))
  scenario_text("bselmodel-summary-heterogeneity-mv-V-approximate", summary_heterogeneity(fit_bselmodel.mv_V_approximate))
  scenario_text("bselmodel-summary-heterogeneity-study-approximate", summary_heterogeneity(fit_bselmodel.mv_V_no_effect_approximate))
  scenario_text("bselmodel-summary-heterogeneity-fixed-approximate", summary_heterogeneity(fit_bselmodel.mv_V_fixed_approximate))

  ### Random effects ----
  ranef_bselmodel.mv_V_exact <- scenario_time("ranef_bselmodel.mv_V_exact", ranef(fit_bselmodel.mv_V_exact))
  ranef_bselmodel.mv_exact <- scenario_time("ranef_bselmodel.mv_exact", ranef(fit_bselmodel.mv_exact))
  ranef_bselmodel_cluster_exact <- scenario_time("ranef_bselmodel_cluster_exact", ranef(fit_bselmodel_cluster_exact))
  ranef_bselmodel.mv_V_approximate <- scenario_time("ranef_bselmodel.mv_V_approximate", ranef(fit_bselmodel.mv_V_approximate))
  ranef_bselmodel.mv_approximate <- scenario_time("ranef_bselmodel.mv_approximate", ranef(fit_bselmodel.mv_approximate))
  ranef_bselmodel_cluster_approximate <- scenario_time("ranef_bselmodel_cluster_approximate", ranef(fit_bselmodel_cluster_approximate))
  scenario_text("bselmodel-random-effects-exact", head(as.data.frame(ranef_bselmodel.mv_V_exact), 12L))
  scenario_text("bselmodel-random-effects-approximate", head(as.data.frame(ranef_bselmodel.mv_V_approximate), 12L))

  plot_bselmodel_ranef_equivalence <- function(ranef_mv, ranef_cluster) {

    par(mfrow = c(1, 2))
    scenario_agreement_plot(as.data.frame(ranef_mv$study)[["Mean"]], as.data.frame(ranef_cluster$cluster)[["Mean"]], main = "study")
    scenario_agreement_plot(as.data.frame(ranef_mv$esid_study)[["Mean"]], as.data.frame(ranef_cluster$estimate)[["Mean"]], main = "esid_study")
    return(invisible(NULL))
  }
  scenario_plot("bselmodel-ranef-mv-vi-exact", plot_bselmodel_ranef_equivalence(ranef_bselmodel.mv_exact, ranef_bselmodel_cluster_exact))
  scenario_plot("bselmodel-ranef-mv-vi-approximate", plot_bselmodel_ranef_equivalence(ranef_bselmodel.mv_approximate, ranef_bselmodel_cluster_approximate))

  # Release fits that are no longer needed before the memory-intensive
  # deletion diagnostics. Scenario caches retain every fit on disk.
  diagnostic_fit_names <- c(
    "fit_bselmodel_cluster_exact",
    "fit_bselmodel_cluster_approximate",
    "fit_bselmodel.mv_exact",
    "fit_bselmodel.mv_approximate",
    "fit_bselmodel.mv_V_exact",
    "fit_bselmodel.mv_V_approximate",
    "fit_bselmodel.mv_V_reg_exact",
    "fit_bselmodel.mv_V_reg_approximate"
  )
  rm(
    list  = setdiff(bselmodel_fit_names, diagnostic_fit_names),
    envir = environment()
  )
  rm(bselmodel_fits)
  invisible(gc())

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
  scenario_plot("bselmodel-marginal-diagnostics-mv-vi-exact", plot_bselmodel_diagnostic_equivalence(fit_bselmodel_cluster_exact, fit_bselmodel.mv_exact))
  scenario_plot("bselmodel-marginal-diagnostics-mv-vi-approximate", plot_bselmodel_diagnostic_equivalence(fit_bselmodel_cluster_approximate, fit_bselmodel.mv_approximate))
  scenario_text("bselmodel-diagnostics-mv-V", head(data.frame(
    exact_residual       = residuals(fit_bselmodel.mv_V_exact),
    approximate_residual = residuals(fit_bselmodel.mv_V_approximate),
    exact_student        = rstudent(fit_bselmodel.mv_V_exact)[["z"]],
    approximate_student  = rstudent(fit_bselmodel.mv_V_approximate)[["z"]]
  ), 12L))
  scenario_text(
    "bselmodel-selection-approximation-diagnostics-mv-V",
    selection_approximation_diagnostics(
      fit_bselmodel.mv_V_approximate,
      max_posterior_samples = 128L,
      latent_samples        = 256L,
      seed                  = 1L
    )
  )

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
