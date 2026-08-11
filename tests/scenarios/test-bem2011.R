if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
REGENERATE_SCENARIO_FILES <- FALSE
SHOW_SCENARIO_OUTPUT      <- FALSE
scenario_start("bem2011")
# testthat::test_file("tests/scenarios/test-bem2011.R")

### Description
# use the bem dataset to fully test the feature suit for publication bias-adjusted meta-analysis

testthat::test_that("Bem simple models", {
  set.seed(1)

  data(Bem2011, package = "RoBMA")

  fit_simple_metafor <- metafor::rma(yi = d, sei = se, data = Bem2011, method = "ML")
  fit_3PSM_metafor   <- metafor::selmodel(fit_simple_metafor, type = "step", steps = c(0.025))
  fit_PET_metafor    <- metafor::rma(yi = d, sei = se, mods = ~ se, data = Bem2011, method = "ML")

  fit_simple <- scenario_fit("fit_simple", {
    tmp <- brma(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1)
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_3PSM <- scenario_fit("fit_3PSM", {
    tmp <- bselmodel(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1)
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })
  fit_PET <- scenario_fit("fit_PET", {
    tmp <- bPET(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1,
                sample = 20000, adapt = 5000, burnin = 5000)
    tmp <- add_loo(tmp)
    tmp <- add_marglik(tmp)
    return(tmp)
  })

  ### simple summary ----
  scenario_text("fit_simple_summary", {summary(fit_simple)})
  scenario_text("fit_3PSM_summary", {summary(fit_3PSM)})
  scenario_text("fit_PET_summary", {summary(fit_PET)})

  # simple model comparison
  scenario_text("fit_bf_comparison", {post_prob(fit_simple, fit_3PSM, fit_PET)})
  scenario_text("fit_loo_comparison", {loo_model_weights(fit_simple, fit_3PSM, fit_PET)})

  ### basic fit plots ----
  set.seed(1)
  scenario_plot("fit_posterior_mu", {
    plot(fit_simple, "mu", ylim = c(0, 12.5), xlim = c(-0.5, 0.5), prior = TRUE)
    lines(fit_simple, "mu", density_method = "IWMDE", lty = 2)

    lines(fit_3PSM, "mu", col = "blue")
    lines(fit_3PSM, "mu", density_method = "IWMDE", lty = 2, col = "blue")

    lines(fit_PET, "mu", col = "red")
    lines(fit_PET, "mu", density_method = "qCMDE", lty = 2, col = "red", density_control = list(samples = 2000))
  })

  ### bias specific plots ----
  set.seed(1)
  scenario_plot("fit_weight", {plot_weightfunction(fit_3PSM)})
  scenario_plot("fit_PET", {plot_pet_peese(fit_PET, ylim = c(-1, 5), xlim = c(0, 1))})

  ### fit specific plots ----
  set.seed(1)
  z_fit_simple <- as_zplot(fit_simple)
  z_fit_3PSM   <- as_zplot(fit_3PSM)
  z_fit_PET    <- as_zplot(fit_PET)

  scenario_plot("zplot", {
    hist(z_fit_3PSM, from = -4)

    lines(z_fit_simple)
    lines(z_fit_3PSM, col = "blue")
    lines(z_fit_PET, col = "red")
  })

  set.seed(1)
  scenario_plot("funnel_simple", {funnel(fit_simple)})
  scenario_plot("funnel_3PSM", {funnel(fit_3PSM)})
  scenario_plot("funnel_PET", {funnel(fit_PET)})

  set.seed(1)
  scenario_plot("qqnorm_simple", {qqnorm(fit_simple)})
  scenario_plot("qqnorm_3PSM", {qqnorm(fit_3PSM)})
  scenario_plot("qqnorm_PET", {qqnorm(fit_PET)})

  ### point hypothesis tests for publication bias adjustment ----
  set.seed(1)
  mu0_seq <- c(-0.25, 0, 0.25)
  fit_3PSM_null_list <- list()
  fit_PET_null_list  <- list()
  for(i in seq_along(mu0_seq)){
    fit_3PSM_null_list[[i]] <- scenario_fit(paste0("fit_3PSM_null_", mu0_seq[i]), {
      tmp <- bselmodel(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1,
                  prior_effect = prior("spike", list(mu0_seq[i])))
      tmp <- add_marglik(tmp)
      return(tmp)
    })
    fit_PET_null_list[[i]] <- scenario_fit(paste0("fit_PET_null_", mu0_seq[i]), {
      tmp <- bPET(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1,
                       prior_effect = prior("spike", list(mu0_seq[i])))
      tmp <- add_marglik(tmp)
      return(tmp)
    })
  }

  names(fit_3PSM_null_list) <- paste0("mu=",mu0_seq)
  names(fit_PET_null_list)  <- paste0("mu=",mu0_seq)

  # get BFs
  BFs_3PSM_marglik <- sapply(fit_3PSM_null_list, function(fit0) bf(fit_3PSM, fit0))
  BFs_PET_marglik  <- sapply(fit_PET_null_list,  function(fit0) bf(fit_PET,  fit0))

  # compute via density methods
  set.seed(1)
  BFs_3PSM_IWMDE  <- sapply(mu0_seq, function(mu0) hypothesis(fit_3PSM, hypothesis = paste0("mu=",mu0), density_method = "IWMDE"))
  BFs_3PSM_qCMDE  <- sapply(mu0_seq, function(mu0) hypothesis(fit_3PSM, hypothesis = paste0("mu=",mu0), density_method = "qCMDE"))
  BFs_3PSM_normal <- sapply(mu0_seq, function(mu0) hypothesis(fit_3PSM, hypothesis = paste0("mu=",mu0), density_method = "normal"))

  BFs_PET_IWMDE  <- sapply(mu0_seq, function(mu0) hypothesis(fit_PET, hypothesis = paste0("mu=",mu0), density_method = "IWMDE"))
  BFs_PET_qCMDE  <- sapply(mu0_seq, function(mu0) hypothesis(fit_PET, hypothesis = paste0("mu=",mu0), density_method = "qCMDE"))
  BFs_PET_normal <- sapply(mu0_seq, function(mu0) hypothesis(fit_PET, hypothesis = paste0("mu=",mu0), density_method = "normal"))

  # compare
  scenario_plot("mu_BF_3PSM_comparison", {
    # the reported error% is large for the first two points
    plot(mu0_seq, log(unlist(BFs_3PSM_marglik[1,])), type = "b", ylab = "logBF", ylim = c(-1, 10))
    lines(mu0_seq + 0.01, log(unlist(BFs_3PSM_IWMDE[3,])),  lty = 2)
    lines(mu0_seq + 0.02, log(unlist(BFs_3PSM_qCMDE[3,])),  lty = 3)
    lines(mu0_seq + 0.03, log(unlist(BFs_3PSM_normal[3,])), lty = 4)
  })
  scenario_plot("mu_BF_PET_comparison", {
    plot(mu0_seq, log(unlist(BFs_PET_marglik[1,])), type = "b", ylab = "logBF")
    lines(mu0_seq + 0.01, log(unlist(BFs_PET_IWMDE[3,])),  lty = 2)
    lines(mu0_seq + 0.02, log(unlist(BFs_PET_qCMDE[3,])),  lty = 3)
    lines(mu0_seq + 0.03, log(unlist(BFs_PET_normal[3,])), lty = 4)
  })

})

testthat::test_that("Bem BMA models", {
  set.seed(1)

  data(Bem2011, package = "RoBMA")

  fit_BMA <- scenario_fit("fit_BMA", {
    tmp <- BMA(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_RoBMA <- scenario_fit("fit_RoBMA", {
    tmp <- RoBMA(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1,
                 adapt = 5000, burnin = 5000, sample = 20000)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_BMA_con <- scenario_fit("fit_BMA_con", {
    tmp <- BMA(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1,
               prior_effect_null = FALSE)
    tmp <- add_loo(tmp)
    return(tmp)
  })
  fit_RoBMA_con <- scenario_fit("fit_RoBMA_con", {
    tmp <- RoBMA(yi = d, sei = se, data = Bem2011, measure = "SMD", seed = 1,
                 prior_effect_null = FALSE,
                 adapt = 5000, burnin = 5000, sample = 20000)
    tmp <- add_loo(tmp)
    return(tmp)
  })

  ### simple summary ----
  scenario_text("fit_BMA_summary",   {summary(fit_BMA)})
  scenario_text("fit_RoBMA_summary", {summary(fit_RoBMA)})
  scenario_text("fit_BMA_summary_con",   {summary(fit_BMA_con)})
  scenario_text("fit_RoBMA_summary_con", {summary(fit_RoBMA_con)})
  scenario_text("fit_BMA_summary_con2",   {summary(fit_BMA, conditional = TRUE)})
  scenario_text("fit_RoBMA_summary_con2", {summary(fit_RoBMA, conditional = TRUE)})

  scenario_text("fit_BMA_summary_models",   {summary_models(fit_BMA)})
  scenario_text("fit_RoBMA_summary_models", {summary_models(fit_RoBMA)})


  ### basic fit plots ----
  set.seed(1)
  scenario_plot("fit_average_posterior_mu", {
    par(mar = c(4, 4, 1, 4))
    plot(fit_BMA, "mu", ylim = c(0, 12.5), xlim = c(-0.5, 0.5), ylim2 = c(0, 1), prior = TRUE)
    lines(fit_BMA, "mu", density_method = "qCMDE", lty = 2)

    lines(fit_RoBMA, "mu", col = "blue")
    lines(fit_RoBMA, "mu", density_method = "qCMDE", lty = 2, col = "blue", density_control = list(samples = 2000))
  })

  set.seed(1)
  scenario_plot("fit_conditional_posterior_mu", {
    plot(fit_BMA, "mu", ylim = c(0, 12.5), xlim = c(-0.5, 0.5), prior = TRUE, conditional = TRUE)
    lines(fit_BMA, "mu", density_method = "IWMDE", lty = 2, conditional = TRUE)
    lines(fit_BMA_con, "mu", density_method = "IWMDE", lty = 3)

    lines(fit_RoBMA, "mu", col = "blue", conditional = TRUE)
    # fails on ESS
    # lines(fit_RoBMA, "mu", density_method = "qCMDE", lty = 2, col = "blue", density_control = list(samples = Inf), conditional = TRUE)
    lines(fit_RoBMA_con, "mu", density_method = "qCMDE", lty = 3, col = "blue", density_control = list(samples = 2000))
  })

  ### bias specific plots ----
  set.seed(1)
  scenario_plot("fit_RoBMA_weight", {plot_weightfunction(fit_RoBMA)})
  scenario_plot("fit_RoBMA_PET",    {plot_pet_peese(fit_RoBMA, ylim = c(-1, 5), xlim = c(0, 1))})

  ### fit specific plots ----
  set.seed(1)
  z_fit_BMA   <- as_zplot(fit_BMA)
  z_fit_RoBMA <- as_zplot(fit_RoBMA)

  scenario_plot("zplot_averaged", {
    hist(z_fit_RoBMA, from = -4)

    lines(z_fit_BMA)
    lines(z_fit_RoBMA, col = "blue")
  })

  set.seed(1)
  scenario_plot("funnel_BMA",   {funnel(fit_BMA)})
  scenario_plot("funnel_RoBMA", {funnel(fit_RoBMA)})

  set.seed(1)
  scenario_plot("qqnorm_BMA",   {qqnorm(fit_BMA)})
  scenario_plot("qqnorm_RoBMA", {qqnorm(fit_RoBMA)})

  ### hypotheses ----

  # full fit product space BF = 0.664 (fit_RoBMA)
  scenario_text("density_test_at_0", rbind(
    hypothesis(fit_RoBMA,     hypothesis = "mu = 0", density_method = "qCMDE", density_control = list(samples = 2000), conditional = TRUE),
    hypothesis(fit_RoBMA_con, hypothesis = "mu = 0", density_method = "qCMDE", density_control = list(samples = 2000))
  ))

  scenario_text("density_ineq", rbind(
    suppressWarnings(hypothesis(fit_BMA, hypothesis = "mu > 0.2 vs mu < 0.2")), # this one should be different because rellies on model-averaged fit
    hypothesis(fit_BMA,     hypothesis = "mu > 0.2 vs mu < 0.2", conditional = TRUE),
    hypothesis(fit_BMA_con, hypothesis = "mu > 0.2 vs mu < 0.2")
  ))

  # check loo comparison works
  scenario_text("fit_loo_comparison_RoBMA", {loo_model_weights(fit_BMA, fit_BMA_con, fit_RoBMA, fit_RoBMA_con)})

})
