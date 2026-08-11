if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
REGENERATE_SCENARIO_FILES <- FALSE
SHOW_SCENARIO_OUTPUT      <- FALSE
scenario_start("bcg")
# testthat::test_file("tests/scenarios/bcg.R")

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
  scenario_text("fit_simple_summary", {print(summary(fit_simple))})
  scenario_text("fit_simple_summary_heterogeneity", {print(summary_heterogeneity(fit_simple))})
  scenario_text("fit_simple_print_prior", {print_prior(fit_simple)})
  scenario_text("fit_simple_interpret", {print(interpret(fit_simple))})

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

  scenario_text("fit_simple_BF_marglik", {print(rbind(BF_minus0, BF_plus0, BF_minusplus))})
  scenario_text("fit_simple_BF_IWMDE",   {print(hypothesis(fit_simple, hypothesis = c("mu < -0.5 vs mu = 0", "mu > -0.5 vs mu = 0", "mu < -0.5 vs mu > -0.5"), density_method = "IWMDE"))})
  scenario_text("fit_simple_BF_qCMDE",   {print(hypothesis(fit_simple, hypothesis = c("mu < -0.5 vs mu = 0", "mu > -0.5 vs mu = 0", "mu < -0.5 vs mu > -0.5"), density_method = "qCMDE"))})
  scenario_text("fit_simple_BF_normal",  {print(hypothesis(fit_simple, hypothesis = c("mu < -0.5 vs mu = 0", "mu > -0.5 vs mu = 0", "mu < -0.5 vs mu > -0.5"), density_method = "normal"))}) # unsurprisingly much more off


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
  scenario_text("fit_simple_influence", {print(influence(fit_simple))})


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
  scenario_text("fit_simple_fe", {print(summary(fit_simple_fe))})
  scenario_text("fit_simple_fe0", {print(summary(fit_simple_fe0))})

  # tests
  BF_fe0 <- bf(fit_simple_fe, fit_simple_fe0)

  scenario_text("fit_simple0_BF_marglik", {print(BF_fe0)})
  scenario_text("fit_simple0_BF_IWMDE",   {print(hypothesis(fit_simple_fe, hypothesis = c("mu = 0"), density_method = "IWMDE"))})
  scenario_text("fit_simple0_BF_qCMDE",   {print(hypothesis(fit_simple_fe, hypothesis = c("mu = 0"), density_method = "qCMDE"))})

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
  scenario_text("fit_reg1_summary", {print(summary(fit_reg1))})
  scenario_text("fit_reg1_summary_heterogeneity", {print(summary_heterogeneity(fit_reg1))})
  scenario_text("fit_reg1_emm", {print(fit_reg1_emm)})
  scenario_text("fit_reg1_emm0", {print(fit_reg1_emm0)})
  scenario_text("fit_reg1_print_prior", {print_prior(fit_reg1)})
  scenario_text("fit_reg1_interpret", {print(interpret(fit_reg1))})

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
  scenario_text("fit_reg1_BF_marglik", {print(BF_nullnull)})
  scenario_text("fit_reg1_BF_default", {print(hypothesis(fit_reg1, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0", "alloc[random] < 0 vs alloc[random] > 0")))})
  scenario_text("fit_reg1_BF_IWMDE",   {print(hypothesis(fit_reg1, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0", "alloc[random] < 0 vs alloc[random] > 0"), density_method = "IWMDE"))})
  scenario_text("fit_reg1_BF_qCMDE",   {print(hypothesis(fit_reg1, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0", "alloc[random] < 0 vs alloc[random] > 0"), density_method = "qCMDE"))})

  # FUTURE:
  # Cross-level point equality is a pending grammar/atom-semantics decision;
  # see .agents/instructions-decisions.md.
  # hypothesis(fit_reg1, hypothesis = c("alloc[random] < alloc[systematic] vs alloc[random] = alloc[systematic]"))
  # hypothesis(fit_reg1, hypothesis = c("alloc[random] - alloc[systematic] < 0 vs alloc[random] - alloc[systematic] = 0"))

  ### directional hypothesis tests for marginal means ----

  # these two match because its the default level
  set.seed(1)
  scenario_text("fit_reg1_mm_BF1", {print(rbind(
    hypothesis(fit_reg1_emm, hypothesis = c("alloc[alternate] < 0 vs alloc[alternate] = 0")),
    hypothesis(fit_reg1,     hypothesis = c("intercept < 0 vs intercept = 0"))
  ))})
  # the mm is larger because its additive effect
  scenario_text("fit_reg1_mm_BF2", {print(rbind(
    hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0")),
    hypothesis(fit_reg1,     hypothesis = c("alloc[random] < 0 vs alloc[random] = 0"))
  ))})
  scenario_text("fit_reg1_mm_BF3", {print(rbind(
    hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0")),
    hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0"), density_method = "IWMDE"),
    hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < 0 vs alloc[random] = 0"), density_method = "qCMDE")
  ))})
  scenario_text("fit_reg1_mm_BF4", {print(
    hypothesis(fit_reg1_emm, hypothesis = "alloc[random] < alloc[alternate] vs alloc[random] > alloc[alternate]")
  )})

  # FUTURE:
  # Symbolic right-hand-side equality is a pending grammar/atom-semantics
  # decision; see .agents/instructions-decisions.md.
  # hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < alloc[alternate] vs alloc[random] = alloc[alternate]"))

  # FUTURE: this would be nice to have
  # hypothesis(fit_reg1_emm, hypothesis = c("alloc[random] < alloc[alternate] < alloc[systematic] vs alloc[random] = alloc[alternate] = alloc[systematic] "))


  ### basic fit diagnostics ----
  set.seed(1)
  influence(fit_reg1_metafor)
  scenario_text("fit_reg1_influence", {print(suppressWarnings(influence(fit_reg1)))})

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
  scenario_text("fit_reg2_summary", {print(summary(fit_reg2))})
  scenario_text("fit_reg2_summary_std", {print(summary(fit_reg2, standardized_coefficients = TRUE))})
  scenario_text("fit_reg2_summary_heterogeneity", {print(summary_heterogeneity(fit_reg2))})

  ### basic fit plots ----
  set.seed(1)
  scenario_plot("fit_reg2_posterior_ablat", {
    plot(fit_reg2, "ablat",  prior = TRUE)
    lines(fit_reg2, "ablat", density_method = "IWMDE", lty = 2)
    lines(fit_reg2, "ablat", density_method = "qCMDE", lty = 2)
  })
  scenario_plot("fit_reg2_regplot", {regplot(fit_reg2, mod = "ablat")})

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
  scenario_text("fit_reg2_BF_marglik", {print(rbind(BF_minus0, BF_plus0, BF_minusplus))})
  scenario_text("fit_reg2_BF_IWMDE",   {print(hypothesis(fit_reg2, hypothesis = c("ablat < -0.25 vs ablat = 0", "ablat > -0.25 vs ablat = 0", "ablat < -0.25 vs ablat > -0.25"), density_method = "IWMDE",  standardized_coefficients = TRUE))})
  scenario_text("fit_reg2_BF_qCMDE",   {print(hypothesis(fit_reg2, hypothesis = c("ablat < -0.25 vs ablat = 0", "ablat > -0.25 vs ablat = 0", "ablat < -0.25 vs ablat > -0.25"), density_method = "qCMDE",  standardized_coefficients = TRUE))})
  scenario_text("fit_reg2_BF_normal",  {print(hypothesis(fit_reg2, hypothesis = c("ablat < -0.25 vs ablat = 0", "ablat > -0.25 vs ablat = 0", "ablat < -0.25 vs ablat > -0.25"), density_method = "normal", standardized_coefficients = TRUE))})

  # Non-standardized need to rescale the comparison point
  set.seed(1)
  0.25 / sd(dat$ablat)
  scenario_text("fit_reg2_BF_IWMDE_direct",   {print(hypothesis(fit_reg2, hypothesis = c("ablat < -0.01730933 vs ablat = 0", "ablat > -0.01730933 vs ablat = 0", "ablat < -0.01730933 vs ablat > -0.01730933"), density_method = "IWMDE"))})
  scenario_text("fit_reg2_BF_qCMDE_direct",   {print(hypothesis(fit_reg2, hypothesis = c("ablat < -0.01730933 vs ablat = 0", "ablat > -0.01730933 vs ablat = 0", "ablat < -0.01730933 vs ablat > -0.01730933"), density_method = "qCMDE"))})
  scenario_text("fit_reg2_BF_normal_direct",  {print(hypothesis(fit_reg2, hypothesis = c("ablat < -0.01730933 vs ablat = 0", "ablat > -0.01730933 vs ablat = 0", "ablat < -0.01730933 vs ablat > -0.01730933"), density_method = "normal"))})

})

