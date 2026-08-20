if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("bcg")
# testthat::test_file("tests/scenarios/test-bcg.R")

### Description
# use the bcg dataset to fully test the feature suit for standard meta-analysis and meta-regressions

testthat::test_that("BCG Simple Fits", {
  set.seed(1)

  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)

  fit_simple_metafor <- scenario_fit("fit_simple_metafor", {
    metafor::rma(yi = yi, vi = vi, data = dat, method = "REML")
  })

  fit_simple <- scenario_fit("fit_simple", {
    tmp <- brma(yi = yi, vi = vi, data = dat, measure = "RR", seed = 1)
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })

  ### basic fit summary ----
  set.seed(1)
  fit_simple_metafor
  scenario_text("fit_simple_summary", {summary(fit_simple)})
  scenario_text("fit_simple_summary_heterogeneity", {summary_heterogeneity(fit_simple)})
  scenario_text("fit_simple_print_prior", {print_prior(fit_simple)})
  scenario_text("fit_simple_interpret", {interpret(fit_simple)})

  ### basic fit plots ----
  set.seed(1)
  scenario_plot("fit_simple_posterior_mu", {
    plot(fit_simple, "mu", ylim = c(0, 2.5), prior = TRUE)
    lines(fit_simple, "mu", density_method = "IWMDE", lty = 2)
    lines(fit_simple, "mu", density_method = "qCMDE", lty = 2)
  })

  scenario_plot("fit_simple_posterior_tau", {
    plot(fit_simple, "tau", ylim = c(0, 3), prior = TRUE)
    lines(fit_simple, "tau", density_method = "IWMDE", lty = 2)
    lines(fit_simple, "tau", density_method = "qCMDE", lty = 2)
  })

  # metafor::forest(fit_simple_metafor)
  scenario_plot("fit_simple_forest", {metafor::forest(as_metafor_forest(fit_simple))})

  ### point hypothesis tests ----
  set.seed(1)
  mu0_seq <- c(-1, -0.5, 0, 0.5)
  fit_simple_null_list <- list()
  for(i in seq_along(mu0_seq)){
    fit_simple_null_list[[i]] <- scenario_fit(paste0("fit_simple_null_", mu0_seq[i]), {
      tmp <- brma(yi = yi, vi = vi, data = dat, measure = "RR",
                  prior_effect = prior("spike", list(mu0_seq[i])))
      tmp <- add_loo(tmp)
      tmp <- add_marglik(tmp)
      return(tmp)
    })
  }
  names(fit_simple_null_list) <- paste0("mu=",mu0_seq)

  # get BFs
  BFs_marglik <- sapply(fit_simple_null_list, function(fit0) bf(fit_simple, fit0))

  # compute via density methods
  set.seed(1)
  BFs_IWMDE  <- sapply(mu0_seq, function(mu0) hypothesis(fit_simple, hypothesis = paste0("mu=",mu0), density_method = "IWMDE"))
  BFs_qCMDE  <- sapply(mu0_seq, function(mu0) hypothesis(fit_simple, hypothesis = paste0("mu=",mu0), density_method = "qCMDE"))
  BFs_normal <- sapply(mu0_seq, function(mu0) hypothesis(fit_simple, hypothesis = paste0("mu=",mu0), density_method = "normal"))

  # compare
  scenario_plot("mu_BF_comparison", {
    plot(mu0_seq, log(unlist(BFs_marglik[1,])), type = "b", ylab = "logBF")
    lines(mu0_seq + 0.02, log(unlist(BFs_IWMDE[3,])),  lty = 2)
    lines(mu0_seq + 0.04, log(unlist(BFs_qCMDE[3,])),  lty = 3)
    lines(mu0_seq + 0.06, log(unlist(BFs_normal[3,])), lty = 4)
  })


  # tau
  set.seed(1)
  tau0_seq <- c(0, 0.25, 0.5)
  fit_simple_null_list_tau <- list()
  for(i in seq_along(tau0_seq)){
    fit_simple_null_list_tau[[i]] <- scenario_fit(paste0("fit_simple_null_tau_", tau0_seq[i]), {
      tmp <- brma(yi = yi, vi = vi, data = dat, measure = "RR",
                  prior_heterogeneity = prior("spike", list(tau0_seq[i])))
      tmp <- add_loo(tmp)
      tmp <- add_marglik(tmp)
      return(tmp)
    })
  }
  names(fit_simple_null_list_tau) <- paste0("tau=",tau0_seq)

  # get BFs
  BFs_tau_marglik <- sapply(fit_simple_null_list_tau, function(fit0) bf(fit_simple, fit0))

  # compute via density methods
  set.seed(1)
  BFs_tau_IWMDE  <- sapply(tau0_seq, function(tau0) hypothesis(fit_simple, hypothesis = paste0("tau=",tau0), density_method = "IWMDE"))
  BFs_tau_qCMDE  <- sapply(tau0_seq, function(tau0) hypothesis(fit_simple, hypothesis = paste0("tau=",tau0), density_method = "qCMDE"))
  BFs_tau_normal <- sapply(tau0_seq, function(tau0) hypothesis(fit_simple, hypothesis = paste0("tau=",tau0), density_method = "normal"))

  # compare
  scenario_plot("tau_BF_comparison", {
    plot(tau0_seq, log(unlist(BFs_tau_marglik[1,])), type = "l", ylab = "logBF")
    lines(tau0_seq + 0.02, log(unlist(BFs_tau_IWMDE[3,])), lty = 2)
    lines(tau0_seq + 0.04, log(unlist(BFs_tau_qCMDE[3,])), lty = 3)
    lines(tau0_seq + 0.06, log(unlist(BFs_tau_normal[3,])), lty = 4)
  })

  ### directional hypothesis tests ----
  set.seed(1)
  fit_simple_minus <- scenario_fit("fit_simple_minus", {
    tmp <- brma(yi = yi, vi = vi, data = dat, measure = "RR", prior_effect = prior("normal", list(0, 1), list(-Inf, -0.5)))
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_simple_plus <- scenario_fit("fit_simple_plus", {
    tmp <- brma(yi = yi, vi = vi, data = dat, measure = "RR", prior_effect = prior("normal", list(0, 1), list(-0.5, Inf)))
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_simple_null <- scenario_fit("fit_simple_null", {
    tmp <- brma(yi = yi, vi = vi, data = dat, measure = "RR", prior_effect = prior("spike", list(0)))
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })

  BF_minus0    <- bf(fit_simple_minus, fit_simple_null)
  BF_plus0     <- bf(fit_simple_plus,  fit_simple_null)
  BF_minusplus <- bf(fit_simple_minus, fit_simple_plus)

  scenario_text("fit_simple_BF_marglik", {rbind(BF_minus0, BF_plus0, BF_minusplus)})
  scenario_text("fit_simple_BF_IWMDE",   {hypothesis(fit_simple, hypothesis = c("mu < -0.5 vs mu = 0", "mu > -0.5 vs mu = 0", "mu < -0.5 vs mu > -0.5"), density_method = "IWMDE")})
  scenario_text("fit_simple_BF_qCMDE",   {hypothesis(fit_simple, hypothesis = c("mu < -0.5 vs mu = 0", "mu > -0.5 vs mu = 0", "mu < -0.5 vs mu > -0.5"), density_method = "qCMDE")})
  scenario_text("fit_simple_BF_normal",  {hypothesis(fit_simple, hypothesis = c("mu < -0.5 vs mu = 0", "mu > -0.5 vs mu = 0", "mu < -0.5 vs mu > -0.5"), density_method = "normal")}) # unsurprisingly much more off


  ### diagnostic plots ----
  set.seed(1)
  scenario_plot("fit_simple_radial", {
    par(mfrow = c(1, 2))
    radial(fit_simple)
    metafor::radial(fit_simple_metafor)
  })
  scenario_plot("fit_simple_funnel",  {
    par(mfrow = c(1, 2))
    funnel(fit_simple)
    metafor::funnel(fit_simple_metafor, addtau2 = TRUE)
  })
  scenario_plot("fit_simple_qqnorm",  {
    par(mfrow = c(1, 2))
    qqnorm(fit_simple, ylim = c(-3, 3), xlim = c(-2, 2))
    qqnorm(fit_simple_metafor, ylim = c(-3, 3), xlim = c(-2, 2))
  })

  ### basic fit diagnostics ----
  set.seed(1)
  influence(fit_simple_metafor)
  scenario_text("fit_simple_influence", {influence(fit_simple)})

  ### likelihood weights ----
  # Integer likelihood weights are equivalent to repeating the corresponding
  # likelihood contribution.
  likelihood_weights             <- rep(1L, nrow(dat))
  likelihood_weights[c(4, 8, 12)] <- c(3L, 10L, 4L)
  dat_duplicated                 <- dat[rep(seq_len(nrow(dat)), likelihood_weights), ]

  fit_simple_weighted <- scenario_fit("fit_simple_weighted", {
    tmp <- brma(yi = yi, vi = vi, weights = likelihood_weights, data = dat, measure = "RR", seed = 3)
    tmp <- suppressWarnings(add_loo(tmp))
    return(tmp)
  })
  fit_simple_duplicated <- scenario_fit("fit_simple_duplicated", {
    brma(yi = yi, vi = vi, data = dat_duplicated, measure = "RR", seed = 4)
  })
  # metafor's custom weights are not likelihood weights, so retain both its
  # weighted and duplicated-data fits for comparison.
  fit_simple_weighted_metafor <- scenario_fit("fit_simple_weighted_metafor", {
    metafor::rma(yi = yi, vi = vi, weights = likelihood_weights, data = dat, method = "REML")
  })
  fit_simple_duplicated_metafor <- scenario_fit("fit_simple_duplicated_metafor", {
    metafor::rma(yi = yi, vi = vi, data = dat_duplicated, method = "REML")
  })

  simple_parameters <- c("mu", "tau")
  scenario_text("fit_simple_weighted_comparison", {data.frame(
    RoBMA              = ex(fit_simple, simple_parameters),
    RoBMA_weighted     = ex(fit_simple_weighted, simple_parameters),
    RoBMA_duplicated   = ex(fit_simple_duplicated, simple_parameters),
    metafor            = ex(fit_simple_metafor, simple_parameters),
    metafor_weighted   = ex(fit_simple_weighted_metafor, simple_parameters),
    metafor_duplicated = ex(fit_simple_duplicated_metafor, simple_parameters)
  )})

  fit_simple_influence                  <- suppressWarnings(influence(fit_simple))[["inf"]]
  fit_simple_weighted_influence         <- suppressWarnings(influence(fit_simple_weighted))[["inf"]]
  fit_simple_influence_metafor          <- metafor::influence.rma.uni(fit_simple_metafor)[["inf"]]
  fit_simple_weighted_influence_metafor <- metafor::influence.rma.uni(fit_simple_weighted_metafor)[["inf"]]
  scenario_text("fit_simple_weighted_dffits", {data.frame(
    study            = seq_len(nrow(dat)),
    weight           = likelihood_weights,
    RoBMA            = fit_simple_influence[["dffits"]],
    RoBMA_weighted   = fit_simple_weighted_influence[["dffits"]],
    metafor          = fit_simple_influence_metafor[["dffits"]],
    metafor_weighted = fit_simple_weighted_influence_metafor[["dffits"]]
  )})
  scenario_text("fit_simple_weighted_cooks", {data.frame(
    study            = seq_len(nrow(dat)),
    weight           = likelihood_weights,
    RoBMA            = fit_simple_influence[["cook.d"]],
    RoBMA_weighted   = fit_simple_weighted_influence[["cook.d"]],
    metafor          = fit_simple_influence_metafor[["cook.d"]],
    metafor_weighted = fit_simple_weighted_influence_metafor[["cook.d"]]
  )})
  scenario_text("fit_simple_weighted_hatvalues", {data.frame(
    study            = seq_len(nrow(dat)),
    weight           = likelihood_weights,
    RoBMA            = fit_simple_influence[["hat"]],
    RoBMA_weighted   = fit_simple_weighted_influence[["hat"]],
    metafor          = fit_simple_influence_metafor[["hat"]],
    metafor_weighted = fit_simple_weighted_influence_metafor[["hat"]]
  )})


  ### null model fitting works ----
  set.seed(1)
  fit_simple_fe <- scenario_fit("fit_simple_fe", {
    tmp <- brma(yi = yi, vi = vi, data = dat, measure = "RR", seed = 1, prior_heterogeneity = prior("spike", list(0)))
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_simple_fe0 <- scenario_fit("fit_simple_fe0", {
    tmp <- brma(yi = yi, vi = vi, data = dat, measure = "RR", seed = 1, prior_effect = prior("spike", list(0)), prior_heterogeneity = prior("spike", list(0)))
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })

  # summary
  scenario_text("fit_simple_fe", {summary(fit_simple_fe)})
  scenario_text("fit_simple_fe0", {summary(fit_simple_fe0)})

  # tests
  BF_fe0 <- bf(fit_simple_fe, fit_simple_fe0)

  scenario_text("fit_simple0_BF_marglik", {BF_fe0})
  scenario_text("fit_simple0_BF_IWMDE",   {hypothesis(fit_simple_fe, hypothesis = c("mu = 0"), density_method = "IWMDE")})
  scenario_text("fit_simple0_BF_qCMDE",   {hypothesis(fit_simple_fe, hypothesis = c("mu = 0"), density_method = "qCMDE")})

})

testthat::test_that("BCG Meta-Regression", {

  set.seed(1)
  data(dat.bcg, package = "metadat")

  ### categorical regression ----
  dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)

  fit_reg1_metafor <- scenario_fit("fit_reg1_metafor", {
    metafor::rma(yi = yi, vi = vi, mods = ~ alloc, data = dat, method = "REML")
  })

  fit_reg1 <- scenario_fit("fit_reg1", {
    tmp <- brma(yi = yi, vi = vi, mods = ~ alloc, data = dat, measure = "RR", seed = 2)
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })

  fit_reg1_emm  <- marginal_means(fit_reg1)
  fit_reg1_emm0 <- marginal_means(fit_reg1, bf = TRUE)

  ### basic fit summary ----
  fit_reg1_metafor
  scenario_text("fit_reg1_summary", {summary(fit_reg1)})
  scenario_text("fit_reg1_summary_heterogeneity", {summary_heterogeneity(fit_reg1)})
  scenario_text("fit_reg1_emm", {fit_reg1_emm})
  scenario_text("fit_reg1_emm0", {fit_reg1_emm0})
  scenario_text("fit_reg1_print_prior", {print_prior(fit_reg1)})
  scenario_text("fit_reg1_interpret", {interpret(fit_reg1)})

  ### basic fit plots ----
  set.seed(1)
  scenario_plot("fit_reg1_posterior_mu", {
    plot(fit_reg1, "mu", ylim = c(0, 2), prior = TRUE)
    lines(fit_reg1, "mu", density_method = "IWMDE", lty = 2)
    lines(fit_reg1, "mu", density_method = "qCMDE", lty = 2)
  })

  scenario_plot("fit_reg1_posterior_tau", {
    plot(fit_reg1, "tau", ylim = c(0, 3), prior = TRUE)
    lines(fit_reg1, "tau", density_method = "IWMDE", lty = 2)
    lines(fit_reg1, "tau", density_method = "qCMDE", lty = 2)
  })

  scenario_plot("fit_reg1_posterior_alloc", {
    plot(fit_reg1, "alloc", ylim = c(0, 1.5), prior = TRUE)
    lines(fit_reg1, "alloc", density_method = "IWMDE", lty = 2)
    lines(fit_reg1, "alloc", density_method = "qCMDE", lty = 2)
  })

  scenario_plot("fit_reg1_regplot", {regplot(fit_reg1, mod = "alloc")})
  scenario_plot("fit_reg1_emmplot", {
    plot(fit_reg1_emm, parameter = "alloc", prior = TRUE, xlim = c(-3, 3))
    lines(fit_reg1_emm, "alloc", density_method = "IWMDE", lty = 2)
    lines(fit_reg1_emm, "alloc", density_method = "qCMDE", lty = 2)
  })

  # with transformation
  scenario_plot("fit_reg1_posterior_alloc_exp", {
    plot(fit_reg1, "alloc", ylim = c(0, 2), prior = TRUE, transform = "EXP", xlim = c(0, 4))
    lines(fit_reg1, "alloc", density_method = "IWMDE", lty = 2, transform = "EXP")
  })

  scenario_plot("fit_reg1_emmplot_exp", {
    plot(fit_reg1_emm, "alloc", ylim = c(0, 5), prior = TRUE, transform = "EXP", xlim = c(0, 2))
    lines(fit_reg1_emm, "alloc", density_method = "IWMDE", lty = 2, transform = "EXP")
  })


  ### directional hypothesis tests ----
  set.seed(1)
  # cannot test multiple factor level simultaneously, so we check null model consistency
  fit_reg1_null <- scenario_fit("fit_reg1_null", {
    tmp <- brma(yi = yi, vi = vi, mods = ~ alloc, data = dat, measure = "RR", prior_mods = prior("spike", list(0)))
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_reg1_null_via_simple <- scenario_fit("fit_reg1_null_via_simple", {
    tmp <- brma(yi = yi, vi = vi, data = dat, measure = "RR")
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  BF_nullnull  <- bf(fit_reg1_null, fit_reg1_null_via_simple)

  # temp_draws <- as_draws_df(fit_reg1)
  # mean(temp_draws$`mu_alloc[1]` < 0.0) / (1 - mean(temp_draws$`mu_alloc[1]` < 0.0))
  set.seed(1)
  scenario_text("fit_reg1_BF_marglik", {BF_nullnull})
  scenario_text("fit_reg1_BF_default", {hypothesis(fit_reg1, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0", "alloc[random] < 0 vs alloc[random] > 0"))})
  scenario_text("fit_reg1_BF_IWMDE",   {hypothesis(fit_reg1, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0", "alloc[random] < 0 vs alloc[random] > 0"), density_method = "IWMDE")})
  scenario_text("fit_reg1_BF_qCMDE",   {hypothesis(fit_reg1, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0", "alloc[random] < 0 vs alloc[random] > 0"), density_method = "qCMDE")})

  scenario_text("fit_reg1_BF_cross_level", {rbind(
    # these are computed by re-arranging the levels into a single posterior row
    hypothesis(fit_reg1, hypothesis = "alloc[random] < alloc[systematic] vs alloc[random] = alloc[systematic]"),
    hypothesis(fit_reg1, hypothesis = "alloc[random] - alloc[systematic] < 0 vs alloc[random] - alloc[systematic] = 0"),
    hypothesis(fit_reg1, hypothesis = "alloc[random] < alloc[systematic] vs alloc[random] = alloc[systematic]",         density_method = "IWMDE"),
    hypothesis(fit_reg1, hypothesis = "alloc[random] - alloc[systematic] < 0 vs alloc[random] - alloc[systematic] = 0", density_method = "IWMDE"),
    hypothesis(fit_reg1, hypothesis = "alloc[random] < alloc[systematic] vs alloc[random] = alloc[systematic]",         density_method = "qCMDE")
  )})
  scenario_text("fit_reg1_BF_cross_level_default", {rbind(
    hypothesis(fit_reg1, hypothesis = "alloc[alternate] > alloc[random] vs alloc[alternate] = alloc[random]"),
    hypothesis(fit_reg1, hypothesis = "0                > alloc[random] vs 0 = alloc[random]"),
    hypothesis(fit_reg1, hypothesis = "alloc[alternate] > alloc[random] vs alloc[alternate] = alloc[random]",   density_method = "qCMDE"),
    hypothesis(fit_reg1, hypothesis = "0                > alloc[random] vs 0 = alloc[random]",                  density_method = "qCMDE")
  )})

  ### directional hypothesis tests for marginal means ----

  # these two match because its the default level
  set.seed(1)
  scenario_text("fit_reg1_mm_BF1", {rbind(
    hypothesis(fit_reg1_emm, hypothesis = c("alloc[alternate] < 0 vs alloc[alternate] = 0")),
    hypothesis(fit_reg1,     hypothesis = c("intercept < 0 vs intercept = 0"))
  )})
  # the mm is larger because its additive effect
  scenario_text("fit_reg1_mm_BF2", {rbind(
    hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0")),
    hypothesis(fit_reg1,     hypothesis = c("alloc[random] < 0 vs alloc[random] = 0"))
  )})
  scenario_text("fit_reg1_mm_BF3", {rbind(
    hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0")),
    hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0"), density_method = "IWMDE"),
    hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0"), density_method = "qCMDE")
  )})
  scenario_text("fit_reg1_mm_BF4", {
    hypothesis(fit_reg1_emm, hypothesis = "alloc[random] < alloc[alternate] vs alloc[random] > alloc[alternate]")
  })

  scenario_text("fit_reg1_mm_BF5", {rbind(
    hypothesis(fit_reg1_emm, hypothesis = "alloc[random] < alloc[alternate] vs alloc[random] = alloc[alternate]"),
    hypothesis(fit_reg1_emm, hypothesis = "alloc[random] - alloc[alternate] < 0 vs alloc[random] - alloc[alternate] = 0"),
    hypothesis(fit_reg1_emm, hypothesis = "alloc[random] < alloc[alternate] vs alloc[random] = alloc[alternate]", density_method = "IWMDE"),
    hypothesis(fit_reg1_emm, hypothesis = "alloc[random] < alloc[alternate] vs alloc[random] = alloc[alternate]", density_method = "qCMDE")
  )})

  # FUTURE: this would be nice to have
  # hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < alloc[alternate] < alloc[systematic] vs alloc[random] = alloc[alternate] = alloc[systematic] "))


  ### basic fit diagnostics ----
  set.seed(1)
  influence(fit_reg1_metafor)
  scenario_text("fit_reg1_influence", {suppressWarnings(influence(fit_reg1))})

  scenario_plot("fit_reg1_funnel",  {
    par(mfrow = c(1, 2))
    suppressWarnings(funnel(fit_reg1, ylim = c(1, 0)))
    metafor::funnel(fit_reg1_metafor, ylim = c(1, 0), type = "rstudent")
  })


  ### continuous regression ----
  set.seed(1)
  dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)

  fit_reg2_metafor <- scenario_fit("fit_reg2_metafor", {
    metafor::rma(yi = yi, vi = vi, mods = ~ ablat, data = dat, method = "REML")
  })

  fit_reg2 <- scenario_fit("fit_reg2", {
    tmp <- brma(yi = yi, vi = vi, mods = ~ ablat, data = dat, measure = "RR", seed = 2)
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })

  ### basic fit summary ----
  fit_reg2_metafor
  scenario_text("fit_reg2_summary", {summary(fit_reg2)})
  scenario_text("fit_reg2_summary_std", {summary(fit_reg2, standardized_coefficients = TRUE)})
  scenario_text("fit_reg2_summary_heterogeneity", {summary_heterogeneity(fit_reg2)})

  ### basic fit plots ----
  set.seed(1)
  scenario_plot("fit_reg2_posterior_ablat", {
    plot(fit_reg2, "ablat",  prior = TRUE)
    lines(fit_reg2, "ablat", density_method = "IWMDE", lty = 2)
    lines(fit_reg2, "ablat", density_method = "qCMDE", lty = 2)
  })
  scenario_plot("fit_reg2_regplot", {regplot(fit_reg2, mod = "ablat")})
  scenario_plot("fit_reg2_regplot_si",  {regplot(fit_reg2, mod = "ablat", si = TRUE)})
  scenario_plot("fit_reg2_regplot_si2", {regplot(fit_reg2, mod = "ablat", si = TRUE, sei = 0.1)})

  ### directional hypothesis tests ----
  set.seed(1)
  fit_reg2_minus <- scenario_fit("fit_reg2_minus", {
    tmp <- brma(yi = yi, vi = vi, mods = ~ ablat, data = dat, measure = "RR", prior_mods = prior("normal", list(0, 0.5), list(-Inf, -0.25)))
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_reg2_plus <- scenario_fit("fit_reg2_plus", {
    tmp <- brma(yi = yi, vi = vi, mods = ~ ablat, data = dat, measure = "RR", prior_mods = prior("normal", list(0, 0.5), list(-0.25, Inf)))
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_reg2_null <- scenario_fit("fit_reg2_null", {
    tmp <- brma(yi = yi, vi = vi, mods = ~ ablat, data = dat, measure = "RR", prior_mods = prior("spike", list(0)))
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })

  BF_minus0    <- bf(fit_reg2_minus, fit_reg2_null)
  BF_plus0     <- bf(fit_reg2_plus,  fit_reg2_null)
  BF_minusplus <- bf(fit_reg2_minus, fit_reg2_plus)

  # Important to realize that the comparison is on the standardized coefficients scale
  set.seed(1)
  scenario_text("fit_reg2_BF_marglik", {rbind(BF_minus0, BF_plus0, BF_minusplus)})
  scenario_text("fit_reg2_BF_IWMDE",   {hypothesis(fit_reg2, hypothesis = c("ablat < -0.25 vs ablat = 0", "ablat > -0.25 vs ablat = 0", "ablat < -0.25 vs ablat > -0.25"), density_method = "IWMDE",  standardized_coefficients = TRUE)})
  scenario_text("fit_reg2_BF_qCMDE",   {hypothesis(fit_reg2, hypothesis = c("ablat < -0.25 vs ablat = 0", "ablat > -0.25 vs ablat = 0", "ablat < -0.25 vs ablat > -0.25"), density_method = "qCMDE",  standardized_coefficients = TRUE)})
  scenario_text("fit_reg2_BF_normal",  {hypothesis(fit_reg2, hypothesis = c("ablat < -0.25 vs ablat = 0", "ablat > -0.25 vs ablat = 0", "ablat < -0.25 vs ablat > -0.25"), density_method = "normal", standardized_coefficients = TRUE)})

  # Non-standardized need to rescale the comparison point
  set.seed(1)
  0.25 / sd(dat$ablat)
  scenario_text("fit_reg2_BF_IWMDE_direct",   {hypothesis(fit_reg2, hypothesis = c("ablat < -0.01730933 vs ablat = 0", "ablat > -0.01730933 vs ablat = 0", "ablat < -0.01730933 vs ablat > -0.01730933"), density_method = "IWMDE")})
  scenario_text("fit_reg2_BF_qCMDE_direct",   {hypothesis(fit_reg2, hypothesis = c("ablat < -0.01730933 vs ablat = 0", "ablat > -0.01730933 vs ablat = 0", "ablat < -0.01730933 vs ablat > -0.01730933"), density_method = "qCMDE")})
  scenario_text("fit_reg2_BF_normal_direct",  {hypothesis(fit_reg2, hypothesis = c("ablat < -0.01730933 vs ablat = 0", "ablat > -0.01730933 vs ablat = 0", "ablat < -0.01730933 vs ablat > -0.01730933"), density_method = "normal")})


  ### regression with an interaction ----
  fit_reg3_metafor <- scenario_fit("fit_reg3_metafor", {
    metafor::rma(yi = yi, vi = vi, mods = ~ alloc * ablat, data = dat, method = "REML")
  })

  fit_reg3 <- scenario_fit("fit_reg3", {
    tmp <- brma(yi = yi, vi = vi, mods = ~ alloc * ablat, data = dat, measure = "RR", seed = 2)
    tmp <- suppressWarnings(add_loo(tmp))
    tmp <- add_marglik(tmp)
    return(tmp)
  })

  fit_reg3_summary <- summary(fit_reg3)
  scenario_text("fit_reg3_summary", {fit_reg3_summary})
  reg3_parameters <- rownames(fit_reg3_summary[["estimates_mods"]])
  scenario_text("fit_reg3_coefficient_comparison", {data.frame(
    term      = reg3_parameters,
    RoBMA     = ex(fit_reg3, reg3_parameters),
    metafor   = ex(fit_reg3_metafor, reg3_parameters),
    row.names = NULL
  )})

  fit_reg3_emm <- marginal_means(fit_reg3)
  scenario_text("fit_reg3_marginal_means", {fit_reg3_emm})

  ### interaction hypothesis tests ----
  # Compare the interaction estimates across allocation groups.
  set.seed(1)
  fit_reg3_hypotheses       <- hypothesis(fit_reg3, hypothesis = c("alloc:ablat[random] < 0", "alloc:ablat[random] = 0", "alloc:ablat[random] < alloc:ablat[systematic]"))
  fit_reg3_hypotheses_IWMDE <- hypothesis(fit_reg3, hypothesis = c("alloc:ablat[random] < 0", "alloc:ablat[random] = 0", "alloc:ablat[random] < alloc:ablat[systematic]"), density_method = "IWMDE")
  scenario_text("fit_reg3_interaction_hypotheses", cbind(
    fit_reg3_hypotheses,
    fit_reg3_hypotheses_IWMDE[,3:4]
  ))

  ### numerical diagnostic comparisons ----
  fit_reg3_rstudent          <- suppressWarnings(rstudent(fit_reg3))
  fit_reg3_dfbetas           <- suppressWarnings(dfbetas(fit_reg3))
  fit_reg3_influence_metafor <- metafor::influence.rma.uni(fit_reg3_metafor)
  fit_reg3_dfbetas_metafor   <- metafor::dfbetas.rma.uni(fit_reg3_metafor)

  # RoBMA standardizes continuous moderators before forming interactions.
  dat_scaled       <- dat
  dat_scaled$ablat <- as.numeric(scale(dat_scaled$ablat))
  fit_reg3_vif_metafor <- metafor::rma(yi = yi, vi = vi, mods = ~ alloc * ablat, data = dat_scaled, method = "REML")
  fit_reg3_vif_metafor <- metafor::vif(fit_reg3_vif_metafor, btt = list(alloc = 2:3, ablat = 4, interaction = 5:6))[["vif"]]
  fit_reg3_vif_metafor <- vapply(fit_reg3_vif_metafor, function(x) x[["sif"]], numeric(1))
  fit_reg3_vif_robma   <- vif(fit_reg3, posterior_correlation = FALSE)[["vif"]]
  scenario_text("fit_reg3_vif_comparison", {data.frame(
    term    = fit_reg3_vif_robma[["term"]],
    RoBMA   = fit_reg3_vif_robma[["GVIF^(1/(2*df))"]],
    metafor = unname(fit_reg3_vif_metafor)
  )})

  # These are visual comparisons, not identities: RoBMA's studentized and
  # influence diagnostics use posterior/PSIS rather than classical deletion.
  # The sparse allocation-by-ablat groups also make classical deletion unstable,
  # whereas RoBMA's priors regularize the deleted interaction estimates.
  fit_reg3_diagnostics_robma <- list(
    "Residuals"           = residuals(fit_reg3),
    "Rstandard (z)"       = rstandard(fit_reg3)[["z"]],
    "Rstudent (z)"        = fit_reg3_rstudent[["z"]],
    "Hat values"          = hatvalues(fit_reg3),
    "DFBETAS"             = as.numeric(as.matrix(fit_reg3_dfbetas)),
    "DFFITS"              = suppressWarnings(dffits(fit_reg3)),
    "Cook's distance"     = suppressWarnings(cooks.distance(fit_reg3)),
    "COVRATIO"            = suppressWarnings(covratio(fit_reg3))
  )
  fit_reg3_diagnostics_metafor <- list(
    "Residuals"           = metafor::residuals.rma(fit_reg3_metafor),
    "Rstandard (z)"       = metafor::rstandard.rma.uni(fit_reg3_metafor)[["z"]],
    "Rstudent (z)"        = metafor::rstudent.rma.uni(fit_reg3_metafor)[["z"]],
    "Hat values"          = metafor::hatvalues.rma.uni(fit_reg3_metafor),
    "DFBETAS"             = unlist(unclass(fit_reg3_dfbetas_metafor)[seq_len(ncol(fit_reg3_dfbetas))], use.names = FALSE),
    "DFFITS"              = fit_reg3_influence_metafor[["inf"]][["dffits"]],
    "Cook's distance"     = fit_reg3_influence_metafor[["inf"]][["cook.d"]],
    "COVRATIO"            = fit_reg3_influence_metafor[["inf"]][["cov.r"]]
  )
  scenario_plot("fit_reg3_diagnostic_comparison", {
    layout(matrix(
      c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 8),
      nrow = 3, byrow = TRUE
    ))
    par(mar = c(3.2, 3.2, 2, 0.5), mgp = c(1.9, 0.6, 0))
    for (diagnostic in names(fit_reg3_diagnostics_robma)) {
      metafor_value <- as.numeric(fit_reg3_diagnostics_metafor[[diagnostic]])
      robma_value   <- as.numeric(fit_reg3_diagnostics_robma[[diagnostic]])
      scenario_agreement_plot(metafor_value, robma_value, diagnostic)
      if (identical(diagnostic, names(fit_reg3_diagnostics_robma)[[1]])) {
        legend(
          "topleft", legend = "band: ±0.1 SD(x)",
          fill = "grey90", border = NA, bty = "n", cex = 0.7
        )
      }
    }
  })

  ### direct predictions ----
  newdata <- expand.grid(
    alloc          = levels(factor(dat$alloc)),
    ablat          = c(15, 35, 55),
    KEEP.OUT.ATTRS = FALSE
  )
  newmods <- stats::model.matrix(~ alloc * ablat, data = newdata)[, -1, drop = FALSE]
  scenario_text("fit_reg3_prediction_comparison", {data.frame(
    newdata,
    RoBMA   = data.frame(predict(fit_reg3, newdata = newdata, type = "terms"))[,"Mean"],
    metafor = predict(fit_reg3_metafor, newmods = newmods)[["pred"]]
  )})

  ### marginal means ----
  marginal_draws    <- fit_reg3_emm[["inference"]][["averaged"]][["mu_alloc__xXx__ablat"]]
  marginal_newdata  <- attr(marginal_draws, "data")
  marginal_newdata$ablat <- mean(dat$ablat) + marginal_newdata$ablat * stats::sd(dat$ablat)
  marginal_newmods       <- stats::model.matrix(~ alloc * ablat, data = marginal_newdata)[, -1, drop = FALSE]
  scenario_text("fit_reg3_marginal_means_comparison", {data.frame(
    marginal_newdata,
    RoBMA   = vapply(marginal_draws, mean, numeric(1)),
    metafor = predict(fit_reg3_metafor, newmods = marginal_newmods)[["pred"]]
  )})

})

