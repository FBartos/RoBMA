source("tests/scenarios/helper-scenarios.R")
REGENERATE_SCENARIO_FILES <- FALSE
scenario_start("bcg")

### Description
# use the bcg dataset to fully test the feature suit for standard meta-analysis and meta-regressions

testthat::test_that("BCG", {
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
  fit_simple_metafor
  scenario_text("fit_simple_summary", {print(summary(fit_simple))})
  scenario_text("fit_simple_summary_heterogeneity", {print(summary_heterogeneity(fit_simple))})
  scenario_text("fit_simple_print_prior", {print_prior(fit_simple)})
  scenario_text("fit_simple_interpret", {print(interpret(fit_simple))})

  ### basic fit plots ----
  scenario_plot("fit_simple_posterior_mu", {
    # TODO: there seems to be more density for the IWMDE/qCMDE figures (the lines is on average higher??)
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

  ### hypothesis tests ----
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
  BFs_IWMDE  <- sapply(mu0_seq, function(mu0) hypothesis(fit_simple, hypothesis = paste0("mu=",mu0), density_method = "IWMDE"))
  BFs_qCMDE  <- sapply(mu0_seq, function(mu0) hypothesis(fit_simple, hypothesis = paste0("mu=",mu0), density_method = "qCMDE"))
  BFs_normal <- sapply(mu0_seq, function(mu0) hypothesis(fit_simple, hypothesis = paste0("mu=",mu0), density_method = "normal"))

  # compare
  scenario_plot("mu_BF_comparison", {
    plot(mu0_seq, log(unlist(BFs_marglik[1,])), type = "l", ylab = "logBF")
    lines(mu0_seq + 0.02, log(unlist(BFs_IWMDE[3,])), lty = 2)
    lines(mu0_seq + 0.04, log(unlist(BFs_qCMDE[3,])), lty = 3)
    lines(mu0_seq + 0.06, log(unlist(BFs_normal[3,])), lty = 4)
  })


  # tau
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


  ### diagnostic plots ----
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
  dev.off()

  ### basic fit diagnostics ----
  influence(fit_simple_metafor)
  scenario_text("fit_simple_influence", {print(influence(fit_simple))})


  ### meta-regressions ----

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
  scenario_text("fit_reg1_summary", {print(summary(fit_reg1))})
  scenario_text("fit_reg1_summary_heterogeneity", {print(summary_heterogeneity(fit_reg1))})
  scenario_text("fit_reg1_emm", {print(fit_reg1_emm)})
  scenario_text("fit_reg1_emm0", {print(fit_reg1_emm0)})
  scenario_text("fit_reg1_print_prior", {print_prior(fit_reg1)})
  scenario_text("fit_reg1_interpret", {print(interpret(fit_reg1))})

  ### basic fit plots ----
  scenario_plot("fit_reg1_posterior_mu", {
    # TODO: how is it possible that the MCSE is so high? this should be an extremely simple case!
    plot(fit_reg1, "mu", ylim = c(0, 2), prior = TRUE)
    lines(fit_reg1, "mu", density_method = "IWMDE", lty = 2)
    lines(fit_reg1, "mu", density_method = "qCMDE", lty = 2)
  })

  scenario_plot("fit_reg1_posterior_tau", {
    # TODO: how is it possible that the MCSE is so high? this should be an extremely simple case!
    plot(fit_reg1, "tau", ylim = c(0, 3), prior = TRUE)
    lines(fit_reg1, "tau", density_method = "IWMDE", lty = 2)
    lines(fit_reg1, "tau", density_method = "qCMDE", lty = 2)
  })

  scenario_plot("fit_reg1_posterior_alloc", {
    # TODO: how is it possible that the MCSE is so high? this should be an extremely simple case!
    plot(fit_reg1, "alloc", ylim = c(0, 1.5), prior = TRUE)
    lines(fit_reg1, "alloc", density_method = "IWMDE", lty = 2)
    lines(fit_reg1, "alloc", density_method = "qCMDE", lty = 2)
  })

  scenario_plot("fit_reg1_regplot", {regplot(fit_reg1, mod = "alloc")})
  scenario_plot("fit_reg1_emmplot", {plot(fit_reg1_emm, parameter = "alloc", prior = TRUE, xlim = c(-3, 3))})

  # TODO: extend with BF comparison from different null models once density method for regression is resolved


  ### null model works ----


})

