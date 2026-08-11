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
  scenario_text("fit_simple_summary", {print(summary(fit_simple))})
  scenario_text("fit_3PSM_summary", {print(summary(fit_3PSM))})
  scenario_text("fit_PET_summary", {print(summary(fit_PET))})

  # simple model comparison
  scenario_text("fit_bf_comparison", {print(post_prob(fit_simple, fit_3PSM, fit_PET))})
  scenario_text("fit_loo_comparison", {print(loo_model_weights(fit_simple, fit_3PSM, fit_PET))})

  ### basic fit plots ----
  set.seed(1)
  scenario_plot("fit_posterior_mu", {
    plot(fit_simple, "mu", ylim = c(0, 12.5), xlim = c(-0.5, 0.5), prior = TRUE)
    lines(fit_simple, "mu", density_method = "IWMDE", lty = 2)

    lines(fit_3PSM, "mu", col = "blue")
    lines(fit_3PSM, "mu", density_method = "IWMDE", lty = 2, col = "blue")

    lines(fit_PET, "mu", col = "red")
    lines(fit_PET, "mu", density_method = "qCMDE", lty = 2, col = "red", density_control = list(max_samples = 2000))
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
  scenario_plot("funnel_simple", {
    funnel(fit_simple)
  })
  scenario_plot("funnel_3PSM", {
    funnel(fit_3PSM)
  })
  scenario_plot("funnel_PET", {
    funnel(fit_PET)
  })

  set.seed(1)
  scenario_plot("qqnorm_simple", {
    qqnorm(fit_simple)
  })
  scenario_plot("qqnorm_3PSM", {
    qqnorm(fit_3PSM)
  })
  scenario_plot("qqnorm_PET", {
    qqnorm(fit_PET)
  })

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
